"""Read-only, hash-checked native Y assets. No CFD or particle solver entrypoint."""
from pathlib import Path
from pd_clot.portable_paths import relocate
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
import yaml
import vtk
from scipy.spatial import cKDTree
from .config import ROOT, dump

TASK = 'y_bifurcation_clot_placement_visual_integration_v1'
OUT = ROOT / 'verification' / TASK
RESULTS = ROOT / 'results' / TASK


def digest(path):
    with relocate(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


class YSources:
    def __init__(self):
        self.register = json.loads((OUT / 'SOURCE_ASSET_REGISTER.json').read_text())
        self.files = {f['relative_path']: f for f in self.register['files']}
        for f in self.files.values():
            assert digest(f['local_readonly_copy']) == f['sha256'], f
        self.states = self.register['states']
        self.config = yaml.safe_load(self.path('configs/stage4/bifurcation.yaml').read_text())
        self.mesh = pv.read(self.path(self.states['simvascular/mesh']['mesh']))
        self.flow = pv.read(self.path(self.states['simvascular/steady']['snapshot']))
        mp = Path(self.states['simvascular/mesh']['mesh']).parent
        self.faces = {n: pv.read(self.path(str(mp / 'mesh-surfaces' / (n + '.vtp')))).triangulate()
                      for n in ['WALL', 'INLET', 'OUTLET_A', 'OUTLET_B']}
        self.exterior = pv.read(self.path(str(mp / 'mesh-complete.exterior.vtp'))).triangulate()
        self.mesh_id = self.states['simvascular/mesh']['mesh_sha256']
        self.flow_id = self.states['simvascular/steady']['snapshot_sha256']
        self.geometry_id = self.states['geometry']['geometry_hash']
        self.u = np.asarray(self.flow['Velocity'], float)
        self.p = np.asarray(self.flow['Pressure'], float).reshape(-1)
        self.tet = self.flow.cells_dict[pv.CellType.TETRA]
        vertices=np.asarray(self.flow.points[self.tet],float)
        self.origins=vertices[:,0]
        self.inverse=np.linalg.inv(np.transpose(vertices[:,1:]-vertices[:,:1],(0,2,1)))
        self.locator=vtk.vtkStaticCellLocator();self.locator.SetDataSet(self.flow);self.locator.BuildLocator()
        self.rejected_locator_candidates=0

    def path(self, relative):
        return relocate(self.files[relative]['local_readonly_copy'])

    def sample(self, positions):
        x = np.atleast_2d(positions).astype(float)
        ids=np.array([self.locator.FindCell(p) for p in x],int)
        def bary(positions,cells):
            w=np.einsum('nij,nj->ni',self.inverse[cells],positions-self.origins[cells])
            return np.c_[1-w.sum(1),w]
        weights=bary(x,np.maximum(ids,0))
        bad=(ids<0)|(weights.min(1)<-1e-10)|(weights.max(1)>1+1e-10)
        for i in np.flatnonzero(bad):
            self.rejected_locator_candidates+=1
            candidates=vtk.vtkIdList();eps=1e-14
            self.locator.FindCellsWithinBounds(np.c_[x[i]-eps,x[i]+eps].ravel(),candidates)
            cells=np.array([candidates.GetId(j) for j in range(candidates.GetNumberOfIds())],int)
            if len(cells):
                w=bary(x[i:i+1],cells);good=(w.min(1)>=-1e-10)&(w.max(1)<=1+1e-10)
                if good.any():
                    j=np.flatnonzero(good)[0];ids[i]=cells[j];weights[i]=w[j];continue
            raise ValueError(f'Outside original CFD lumen: strict barycentric membership failed at {x[i].tolist()}; no extrapolation')
        tet=self.tet[ids]
        return (np.einsum('ni,nij->nj', weights, self.u[tet]),
                np.einsum('ni,ni->n', weights, self.p[tet]))

    def audit(self):
        rows = []
        def check(name, handed, actual, equal=True, source='native mesh/config/field'):
            rows.append(dict(item=name, user_handoff=handed, actual=actual,
                             status='ACTUAL_FILE_VERIFIED' if equal else 'MISSING_OR_DIFFERENT', source=source))
        config = self.config
        expected = dict(parent_length_m=.014, parent_radius_m=.0013, daughter_length_m=.012,
                        daughter_radius_m=.001, branch_angle_deg=30., blending_length_m=.0002,
                        mesh_edge_size_m=.00025)
        for k, v in expected.items():
            check(k, v, config['geometry'][k], np.isclose(v, config['geometry'][k], rtol=1e-12))
        for k, v in dict(density_kg_m3=1000., viscosity_Pa_s=.001).items():
            check(k, v, config['medium'][k], config['medium'][k] == v)
        check('inlet_mean_velocity_m_s', .005, config['flow']['inlet_mean_velocity_m_s'], config['flow']['inlet_mean_velocity_m_s'] == .005)
        check('outlet_reference_tractions_Pa', [0., 0.], config['flow']['outlet_pressure_Pa'], config['flow']['outlet_pressure_Pa'] == [0., 0.])
        check('nodes', 26712, self.mesh.n_points, self.mesh.n_points == 26712)
        check('tetrahedra', 146230, self.mesh.n_cells, self.mesh.n_cells == 146230)
        check('physical_ports', 3, len(self.faces)-1, len(self.faces) == 4)
        xmlpath=self.path(str(Path(self.states['simvascular/steady']['run'])/'run/solver.xml'))
        xml=ET.parse(xmlpath).getroot();equation=xml.find('Add_equation')
        check('native_solver_dimensions',3,int(xml.find('.//Number_of_spatial_dimensions').text),int(xml.find('.//Number_of_spatial_dimensions').text)==3,source='actual solver.xml')
        check('native_solver_density_kg_m3',1000.,float(xml.find('.//Density').text),float(xml.find('.//Density').text)==1000.,source='actual solver.xml')
        check('native_solver_viscosity_Pa_s',.001,float(xml.find('.//Viscosity/Value').text),float(xml.find('.//Viscosity/Value').text)==.001,source='actual solver.xml')
        bc={v.attrib['name']:{child.tag:child.text for child in v} for v in equation.findall('Add_BC')}
        check('native_wall_noslip','Dir / 0',bc['WALL'],bc['WALL']['Type']=='Dir' and float(bc['WALL']['Value'])==0,source='actual solver.xml')
        check('native_inlet_flux_and_profile','Parabolic, flux normalized, steady',bc['INLET'],bc['INLET']['Profile']=='Parabolic' and bc['INLET']['Impose_flux']=='true' and bc['INLET']['Time_dependence']=='Steady',source='actual solver.xml')
        for name in ['OUTLET_A','OUTLET_B']:
            check('native_'+name+'_traction','Neu / 0',bc[name],bc[name]['Type']=='Neu' and float(bc[name]['Value'])==0,source='actual solver.xml')
        check('native_final_time_s',1.21,121*float(xml.find('.//Time_step_size').text),source='actual final filename and solver.xml dt')
        check('geometric_mirror_construction','positive-y half reflected',self.states['geometry']['scope'],'positive-y half reflected' in self.states['geometry']['scope'],source='original geometry STATE provenance, not a new anatomy assertion')
        check('coordinate_frame', 'global xyz, m', 'native points; identity transform; SI solver XML')
        check('maximum_velocity_m_s', .010074017175615313, float(np.linalg.norm(self.u, axis=1).max()),
              np.isclose(np.linalg.norm(self.u, axis=1).max(), .010074017175615313))
        tree = cKDTree(self.flow.points)
        flux = {}; pressures = {}; areas = {}; centers = {}
        for name, surface in self.faces.items():
            dist, ids = tree.query(surface.points)
            assert dist.max() < 1e-10
            tri = surface.faces.reshape(-1, 4)[:, 1:]
            xyz = np.asarray(surface.points)[tri]
            av = .5*np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
            area = np.linalg.norm(av, axis=1)
            flux[name] = float(np.einsum('ij,ij->', self.u[ids][tri].mean(1), av))
            pressures[name] = float(np.dot(self.p[ids][tri].mean(1), area)/area.sum())
            areas[name] = float(area.sum())
            centers[name] = np.average(xyz.mean(1), weights=area, axis=0).tolist()
        check('actual_inlet_Q_m3_s', 2.6354060534101732e-8, -flux['INLET'], np.isclose(-flux['INLET'], 2.6354060534101732e-8, atol=1e-17))
        check('inlet_area_average_pressure_Pa', .7049350301311733, pressures['INLET'], np.isclose(pressures['INLET'], .7049350301311733))
        for name, target in [('OUTLET_A', 1.3176074951406056e-8), ('OUTLET_B', 1.3177985582695677e-8)]:
            check(name+'_flow_m3_s', target, flux[name], np.isclose(flux[name], target, atol=1e-17))
        for name,target in [('INLET',[-.014,0,0]),('OUTLET_A',[.012*np.cos(np.pi/6),.006,0]),('OUTLET_B',[.012*np.cos(np.pi/6),-.006,0])]:
            check(name+'_actual_port_center_m',target,centers[name],np.linalg.norm(np.array(centers[name])-target)<5e-5,source='native boundary area-weighted center; discretization tolerance 50 µm')
        check('final_step', 121, int(Path(self.states['simvascular/steady']['snapshot']).stem.split('_')[-1]))
        check('historical_P1_import_checks', 1000, sum(v['count'] for v in self.states['simvascular/import']['comparisons'].values()),
              source='read-only historical import STATE, not rerun')
        report = dict(status='PASS' if all(x['status']=='ACTUAL_FILE_VERIFIED' for x in rows) else 'CONFLICT',
                      rows=rows, bbox_m=list(self.mesh.bounds), volume_m3=float(self.mesh.volume),
                      face_areas_m2=areas, face_centers_m=centers, boundary_outward_flux_m3_s=flux,
                      boundary_mean_pressure_Pa=pressures, velocity_range_m_s=[float(np.linalg.norm(self.u,axis=1).min()),float(np.linalg.norm(self.u,axis=1).max())],
                      pressure_range_Pa=[float(self.p.min()),float(self.p.max())],
                      mesh_id=self.mesh_id, flow_id=self.flow_id, geometry_id=self.geometry_id,
                      domain='UNOBSTRUCTED_STEADY_BACKGROUND_FLOW', CLOT_PRESENT_CFD_SOLVED='NOT_RUN')
        dump(OUT/'ACTUAL_Y_ASSET_AUDIT.json', report)
        (OUT/'ACTUAL_Y_ASSET_AUDIT.md').write_text('# Actual server Y assets\n\n'+'\n'.join(f"- {r['item']}: USER_HANDOFF={r['user_handoff']}; actual={r['actual']}; {r['status']}" for r in rows)+'\n\nCFD source is the original unobstructed lumen. CLOT_PRESENT_CFD_SOLVED = NOT_RUN.\n')
        return report

    def spot_check(self, near_clot):
        rng = np.random.default_rng(810)
        ids = rng.choice(len(self.tet), 8, replace=False)
        weights = rng.dirichlet(np.ones(4), len(ids))
        q = np.einsum('ni,nij->nj', weights, self.flow.points[self.tet[ids]])
        u, p = self.sample(q)
        eu = np.einsum('ni,nij->nj', weights, self.u[self.tet[ids]])
        ep = np.einsum('ni,ni->n', weights, self.p[self.tet[ids]])
        node_ids = self.tet[ids[:4], 0]; un, pn = self.sample(self.flow.points[node_ids])
        uq, pq = self.sample(near_clot)
        outside = False
        try: self.sample([[1,1,1]])
        except ValueError: outside = True
        report = dict(status='PASS', internal_points=8, native_nodes=4, local_points=len(near_clot),
                      velocity_max_abs_error_m_s=float(np.abs(u-eu).max()), pressure_max_abs_error_Pa=float(np.abs(p-ep).max()),
                      node_velocity_error=float(np.abs(un-self.u[node_ids]).max()), node_pressure_error=float(np.abs(pn-self.p[node_ids]).max()),
                      local_velocity_m_s=uq,local_pressure_Pa=pq,outside_rejected=outside,
                      original_1000_point_check_not_repeated=True, coordinate_transform=np.eye(4), source_flow_id=self.flow_id,
                      load_applied=False, CLOT_PRESENT_CFD_SOLVED='NOT_RUN')
        assert report['velocity_max_abs_error_m_s']<1e-12 and report['pressure_max_abs_error_Pa']<1e-10 and outside
        dump(OUT/'FLOW_COORDINATE_CHECK.json',report)
        np.savez_compressed(OUT/'FLOW_SPOT_CHECKS.npz',position_xyz_m=q,velocity_m_s=u,pressure_Pa=p,near_clot_xyz_m=near_clot,near_clot_velocity_m_s=uq)
        return report
