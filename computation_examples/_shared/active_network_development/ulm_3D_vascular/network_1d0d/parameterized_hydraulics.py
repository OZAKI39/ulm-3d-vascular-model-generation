"""Geometry parameters -> analytical edge R -> complete steady SI network.

This module never loads a CFD field. Saved H0 *geometry* is a portable input;
the legacy H0 state solution is only a regression oracle in the tests.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from types import MappingProxyType
import hashlib
import json

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .audit import sha256
from .hydraulic_resistance import linear_radius_resistance, MU_PA_S, ROI_TARGET_Q_M3_S
from .idealized_h0 import H0Domain, h0_boundary, port_signs
from .network_solver import solve_network, operating_point_scale

MODEL_NAME = 'PARAMETERIZED_GEOMETRY_HYDRAULIC_0D_V1'
PORT_ORDER = ('INLET', 'O1', 'O2', 'O3')
OUTLET_ORDER = PORT_ORDER[1:]


def readonly(values, dtype=None):
    """An owned immutable snapshot, including protection against setflags(write=True)."""
    a = np.ascontiguousarray(values, dtype=dtype)
    return np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)


def _labels(n, edges):
    u, v = np.asarray(edges, dtype=int).T
    graph = coo_matrix((np.ones(2*len(u)), (np.r_[u,v], np.r_[v,u])), shape=(n,n))
    return connected_components(graph, directed=False)


def load_parameterized_domain(geometry_npz, ports_json):
    """Load the existing SI H0 graph, not H0_solution_si.npz or any 3D result.

    Unlike legacy load_h0(), this loader does not depend on archived absolute
    paths to FEM files. It rechecks geometry/topology and recomputes port signs.
    """
    mapping = json.loads(Path(ports_json).read_text())
    if (mapping.get('status') != 'PASS' or mapping.get('extra_saved_roi_ports') != 0
            or mapping.get('expected_fem_ports') != 4):
        raise ValueError('Unverified four-port geometry mapping')
    # Retain geometry only; extension/FEM/pressure metadata is not a model input.
    ports = [{k:p[k] for k in ('name','network_node_id','outward_normal')}
             for p in mapping['ports']]
    if sorted(p['name'] for p in ports) != sorted(PORT_ORDER):
        raise ValueError('Exactly INLET, O1, O2, O3 required')
    ports.sort(key=lambda p: PORT_ORDER.index(p['name']))
    with np.load(geometry_npz, allow_pickle=False) as g:
        ids, xyz, r, e = [g[k].copy() for k in ('ids','xyz_m','radius_m','edges')]
        inside = g['roi_internal_edge_mask'].copy()
        source, terminals = h0_boundary(ids,e)
        if source != int(g['source_index']) or not np.array_equal(terminals,g['terminal_indices']):
            raise ValueError('Saved source/terminal indices disagree with geometry')
        signs = port_signs(ids,xyz,e,inside,ports)
        pn = np.array([s['node_index'] for s in signs])
        if not np.array_equal(pn,g['port_indices']):
            raise ValueError('Saved port order/indices disagree with verified port mapping')
        _, exterior = _labels(len(ids),e[~inside])
        return H0Domain(ids,xyz,r,e,g['original_edge'].copy(),g['fractions'].copy(),
                        inside,ports,pn,signs,source,terminals,exterior)


def downstream_node_mask(domain, port_name):
    """Exterior connected component, excluding the fixed real-cut node.

    Connectivity is undirected: reversing SWC edge storage cannot change this
    mask. This safely reuses external_radius_variant's physical rule without
    changing that protected legacy function.
    """
    if port_name == 'O3':
        raise ValueError('NOT_AVAILABLE: O3 has no saved downstream geometry')
    if port_name not in ('O1','O2'):
        raise ValueError('Downstream scaling supported only for O1/O2')
    pn = next(s['node_index'] for s in port_signs(domain.ids,domain.xyz_m,domain.edges,
                                               domain.internal_roi,domain.ports) if s['port']==port_name)
    _, labels = _labels(len(domain.ids),domain.edges[~domain.internal_roi])
    mask = labels == labels[pn]
    mask[pn] = False
    if not mask.any() or mask[domain.source] or np.any(mask[np.unique(domain.edges[domain.internal_roi])]):
        raise ValueError('Invalid downstream component: absent, reaches source or crosses ROI')
    if not np.any(mask[domain.terminals]):
        raise ValueError('Downstream component has no reference terminal')
    return readonly(mask)


def downstream_edge_mask(domain, port_name):
    mask = downstream_node_mask(domain,port_name).copy()
    pn = next(p['node_index'] for p in domain.signs if p['port']==port_name)
    mask[pn] = True
    return readonly(~domain.internal_roi & np.all(mask[domain.edges],axis=1))


@dataclass(frozen=True)
class GeometryHydraulicCache:
    node_ids: np.ndarray
    xyz_m: np.ndarray
    baseline_radius_m: np.ndarray
    edge_nodes: np.ndarray
    edge_length_m: np.ndarray
    baseline_r0_m: np.ndarray
    baseline_r1_m: np.ndarray
    baseline_resistance_pa_s_m3: np.ndarray
    internal_roi_edge_mask: np.ndarray
    downstream_node_masks: object
    downstream_edge_masks: object
    source_index: int
    terminal_indices: np.ndarray
    port_indices: np.ndarray
    port_edge_indices: np.ndarray
    port_coefficients: np.ndarray
    input_hashes: object
    geometry_sha256: str

    @classmethod
    def from_domain(cls, domain, input_hashes=None):
        ids, xyz, radii, e = domain.ids,domain.xyz_m,domain.radius_m,domain.edges
        n = len(ids)
        if (len(np.unique(ids))!=n or xyz.shape!=(n,3) or radii.shape!=(n,)
                or not np.all(np.isfinite(xyz)) or not np.all(np.isfinite(radii)) or np.any(radii<=0)):
            raise ValueError('Invalid finite SI geometry/radius or duplicate node ID')
        if (e.ndim!=2 or e.shape[1]!=2 or not np.issubdtype(e.dtype,np.integer)
                or np.any(e<0) or np.any(e>=n) or np.any(e[:,0]==e[:,1])):
            raise ValueError('Invalid edge topology')
        if domain.internal_roi.dtype!=bool or domain.internal_roi.shape!=(len(e),):
            raise ValueError('ROI mask must be one boolean per edge')
        if _labels(n,e)[0]!=1:
            raise ValueError('Source/reference component check failed: network must be connected')
        source, terminals = h0_boundary(ids,e)
        if source!=domain.source or not np.array_equal(terminals,domain.terminals):
            raise ValueError('Source/reference definitions differ from H0')
        ports=sorted(domain.ports,key=lambda p: PORT_ORDER.index(p['name']))
        if tuple(p['name'] for p in ports)!=PORT_ORDER:
            raise ValueError('Duplicate/missing port identity')
        for p in ports:
            normal=np.asarray(p['outward_normal'])
            if normal.shape!=(3,) or not np.isfinite(normal).all() or not np.isclose(np.linalg.norm(normal),1.,rtol=1e-9,atol=1e-12):
                raise ValueError('Port outward normal must be finite unit vector')
        signs=port_signs(ids,xyz,e,domain.internal_roi,ports)
        node_masks={p:downstream_node_mask(domain,p) for p in ('O1','O2')}
        edge_masks={p:downstream_edge_mask(domain,p) for p in ('O1','O2')}
        if np.any(node_masks['O1'] & node_masks['O2']):
            raise ValueError('Overlapping downstream components')
        length=np.linalg.norm(xyz[e[:,1]]-xyz[e[:,0]],axis=1)
        if np.any(length<=0): raise ValueError('Solver edge length must be positive')
        r0,r1=radii[e[:,0]],radii[e[:,1]]
        resistance=linear_radius_resistance(length,r0,r1,MU_PA_S)
        if np.any(resistance<=0):raise ValueError('Positive finite baseline resistance required')
        digest=hashlib.sha256()
        for a in (ids,xyz,radii,e,domain.internal_roi):
            digest.update(str((a.shape,a.dtype.str)).encode());digest.update(np.ascontiguousarray(a).tobytes())
        digest.update(json.dumps(signs,sort_keys=True).encode())
        return cls(readonly(ids),readonly(xyz),readonly(radii),readonly(e),readonly(length),
                   readonly(r0),readonly(r1),readonly(resistance),readonly(domain.internal_roi),
                   MappingProxyType(node_masks),MappingProxyType(edge_masks),source,readonly(terminals),
                   readonly([s['node_index'] for s in signs],int),readonly([s['roi_segment'] for s in signs],int),
                   readonly([s['coefficient'] for s in signs],int),MappingProxyType(dict(input_hashes or {})),digest.hexdigest())


def load_geometry_cache(geometry_npz, ports_json):
    domain=load_parameterized_domain(geometry_npz,ports_json)
    return GeometryHydraulicCache.from_domain(domain,{
        'geometry_npz_sha256':sha256(geometry_npz),'port_mapping_sha256':sha256(ports_json)})


@dataclass(frozen=True)
class ParameterizedHydraulicSpec:
    s_O1: float = 1.0
    s_O2: float = 1.0
    mu_pa_s: float = MU_PA_S
    distal_reference_pa: float = 0.0
    o3_mode: str = 'legacy_reference'
    terminal_resistance_O3_pa_s_m3: float | None = None

    def validate(self):
        for name in ('s_O1','s_O2','mu_pa_s'):
            x=getattr(self,name)
            if not np.isfinite(x) or x<=0: raise ValueError(name+' must be finite positive')
        if not np.isfinite(self.distal_reference_pa): raise ValueError('Nonfinite reference pressure')
        if self.o3_mode not in ('legacy_reference','terminal_resistance'): raise ValueError('Unsupported O3 mode')
        if self.o3_mode=='legacy_reference' and self.terminal_resistance_O3_pa_s_m3 is not None:
            raise ValueError('Legacy O3 must not specify a terminal resistance')
        if self.o3_mode=='terminal_resistance':
            r=self.terminal_resistance_O3_pa_s_m3
            if r is None or not np.isfinite(r) or r<=0:
                raise ValueError('O3 terminal resistance must be finite positive')


@dataclass(frozen=True)
class ParameterizedHydraulicResult:
    spec: ParameterizedHydraulicSpec
    target_roi_flow_m3_s: float
    pressure_pa: np.ndarray
    edge_flow_m3_s: np.ndarray
    node_outflow_m3_s: np.ndarray
    effective_radius_m: np.ndarray
    edge_resistance_pa_s_m3: np.ndarray
    edge_nodes: np.ndarray
    port_pressure_pa: np.ndarray
    port_flow_m3_s: np.ndarray
    signed_outlet_fractions: np.ndarray
    positive_outlet_fractions: np.ndarray | None
    mass_audit: dict
    reverse_flow_audit: dict
    solver_audit: dict
    radius_provenance: dict
    closure_provenance: dict
    geometry_sha256: str
    input_hashes: dict

    def summary(self):
        return dict(model_name=MODEL_NAME,parameters=asdict(self.spec),
            roi_target_flow_m3_s=self.target_roi_flow_m3_s,
            inlet_pressure_pa=float(self.port_pressure_pa[0]),
            real_cut_pressure_pa=dict(zip(PORT_ORDER,self.port_pressure_pa.tolist())),
            port_flow_m3_s=dict(zip(PORT_ORDER,self.port_flow_m3_s.tolist())),
            outlet_flow_fraction=dict(zip(OUTLET_ORDER,self.signed_outlet_fractions.tolist())),
            normalized_positive_outlet_fraction=None if self.positive_outlet_fractions is None else dict(zip(OUTLET_ORDER,self.positive_outlet_fractions.tolist())),
            mass_audit=self.mass_audit,reverse_flow_audit=self.reverse_flow_audit,
            solver_audit=self.solver_audit,radius_provenance=self.radius_provenance,
            closure_provenance=self.closure_provenance,source_geometry_sha256=self.geometry_sha256,
            input_hashes=self.input_hashes)


def solve_parameterized_operating_point(cache, spec, target_roi_flow_m3_s=ROI_TARGET_Q_M3_S):
    """Independently construct R, assemble/solve the whole network, then scale.

    Solve relative to a common Pd first. Adding Pd after scaling is a gauge
    shift; multiplying Pd by the operating-point factor would be incorrect.
    """
    spec.validate()
    radii=cache.baseline_radius_m.copy()
    for name in ('O1','O2'): radii[cache.downstream_node_masks[name]]*=getattr(spec,'s_'+name)
    e=cache.edge_nodes
    resistance=linear_radius_resistance(cache.edge_length_m,radii[e[:,0]],radii[e[:,1]],spec.mu_pa_s)
    n=len(radii); terminals=cache.terminal_indices
    bc={int(i):0. for i in terminals};bc[cache.source_index]=1.
    closure=dict(mode='legacy_reference',source='DIRICHLET_REFERENCE_TERMINAL',
                 reference_pressure_pa=spec.distal_reference_pa,virtual_nodes=[],virtual_edges=[])
    if spec.o3_mode=='terminal_resistance':
        o3=int(cache.port_indices[3])
        if o3 not in terminals: raise ValueError('O3 must be the existing reference terminal')
        e=np.vstack([e,[o3,n]])
        resistance=np.r_[resistance,spec.terminal_resistance_O3_pa_s_m3]
        terminals=np.r_[terminals[terminals!=o3],n]
        del bc[o3];bc[n]=0.
        closure=dict(mode='terminal_resistance',source='MODEL_CLOSURE_PARAMETER',
            terminal_resistance_pa_s_m3=spec.terminal_resistance_O3_pa_s_m3,
            reference_pressure_pa=spec.distal_reference_pa,original_node_count=n,
            virtual_nodes=[dict(augmented_index=n,name='DISTAL_O3',original_swc_node_id=None)],
            virtual_edges=[dict(augmented_index=len(e)-1,node_indices=[o3,n],original_swc_edge_id=None)])
        n+=1
    unit=solve_network(n,e,resistance,bc)
    unit_ports=unit.flow[cache.port_edge_indices]*cache.port_coefficients
    scale=operating_point_scale(unit_ports[0],target_roi_flow_m3_s)
    p=unit.pressure*scale+spec.distal_reference_pa
    q=unit.flow*scale;nodal=unit.node_outflow*scale;portq=unit_ports*scale
    fixed=np.array(sorted(bc));free=np.setdiff1d(np.arange(n),fixed)
    sourceq=float(nodal[cache.source_index])
    if not np.isfinite(sourceq) or sourceq<=0: raise ValueError('Invalid source flow orientation')
    boundary=float(abs(nodal[fixed].sum())/sourceq)
    internal=float(np.max(np.abs(nodal[free]))/sourceq) if len(free) else 0.
    roi=float(abs(portq[0]-portq[1:].sum())/portq[0])
    target_error=float(abs(portq[0]-target_roi_flow_m3_s)/target_roi_flow_m3_s)
    residual=max(boundary,internal,roi,target_error)
    if not np.isfinite(residual) or residual>1e-9: raise ValueError('Parameterized network mass balance failed')
    fractions=portq[1:]/portq[0]
    # Never silently clip backflow. A positive-only summary is unavailable when
    # any outlet reverses; the full signed state and reverse audit remain intact.
    positive=portq[1:]/portq[1:].sum() if np.all(portq[1:]>=0) and portq[1:].sum()>0 else None
    reverse=dict(port_backflow={k:bool(v<0) for k,v in zip(PORT_ORDER,portq)},
                 negative_algebraic_edge_count=int(np.sum(q<0)),
                 edge_sign_meaning='stored u->v only; negative edge Q is not automatically physiological backflow',
                 clipping_applied=False)
    mass=dict(status='PASS',max_relative_residual=residual,global_relative_residual=boundary,
              max_internal_relative_residual=internal,roi_relative_residual=roi,
              roi_target_relative_error=target_error,source_inflow_m3_s=sourceq,
              terminal_outflow_m3_s=float(-nodal[terminals].sum()))
    prov={name:dict(source='EXISTING_GEOMETRY_EFFECTIVE_RADIUS_PARAMETER',
                   radius_scale_dimensionless=getattr(spec,'s_'+name),
                   scaled_nodes=int(cache.downstream_node_masks[name].sum()),
                   affected_edges=int(cache.downstream_edge_masks[name].sum()),
                   real_cut_radius_unchanged=True,first_edge_rule='fixed cut radius to scaled downstream radius; analytical taper')
          for name in ('O1','O2')}
    return ParameterizedHydraulicResult(spec,float(target_roi_flow_m3_s),readonly(p),readonly(q),readonly(nodal),
        readonly(radii),readonly(resistance),readonly(e),readonly(p[cache.port_indices]),readonly(portq),
        readonly(fractions),None if positive is None else readonly(positive),mass,reverse,
        dict(unit.audit,operating_point_scale_pa=scale,fixed_pressure_node_indices=fixed.tolist(),
             pressure_reference_pa=spec.distal_reference_pa),prov,closure,cache.geometry_sha256,dict(cache.input_hashes))
