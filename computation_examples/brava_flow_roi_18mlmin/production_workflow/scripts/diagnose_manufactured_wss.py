"""Raw WSS interpolation diagnostic; keeps the original calibration gates unchanged."""
from pathlib import Path
import argparse
import json
import numpy as np
import pyvista as pv

ROOT=Path(__file__).resolve().parents[1]


def main(name):
    case=ROOT/'manufactured_pipe'/name
    config=json.loads((ROOT/'straight_pipe/manufactured'/name/'INPUT.json').read_text())
    proto=json.loads((ROOT/'manufactured_pipe/PROTOCOL.json').read_text())
    radius,length,mean,mu=config['radius_m'],config['length_m'],proto['mean_velocity_m_s'],proto['mu_Pa_s']
    with np.load(case/'mesh/BOUNDARY_TOPOLOGY.npz') as d:
        p,t=d['points_m'],d['tetra'];mask=d['facet_tags']==1
        wall,owner=d['boundary_triangles'][mask],d['adjacent_tetra'][mask]
    with np.load(case/'cfd/analysis/FIXED_PROBES.npz') as d:
        gradient=d['gradient_all_tetra'];targets=d['points_m'][d['groups']==0].copy()
    exact_nodal=np.zeros_like(p);exact_nodal[:,2]=2*mean*(1-np.sum(p[:,:2]**2,axis=1)/radius**2)
    exact_nodal[np.unique(wall)]=0
    dx=p[t[:,1:]]-p[t[:,:1]]
    gi=np.linalg.solve(dx,exact_nodal[t[:,1:]]-exact_nodal[t[:,:1]]).transpose(0,2,1)
    xyz=p[wall];normals=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);area=np.linalg.norm(normals,axis=1)/2
    normals/=2*area[:,None]
    def tangential(g):
        stress=mu*(g[owner]+g[owner].transpose(0,2,1))
        traction=np.einsum('nij,nj->ni',stress,normals)
        return traction-np.einsum('ni,ni->n',traction,normals)[:,None]*normals
    actual,interpolated=tangential(gradient),tangential(gi)
    w,wi=np.linalg.norm(actual,axis=1),np.linalg.norm(interpolated,axis=1)
    exact=4*mu*mean/radius
    targets[:,:2]*=radius/np.linalg.norm(targets[:,:2],axis=1)[:,None]
    poly=pv.PolyData(p,np.c_[np.full(len(wall),3),wall].ravel())
    selected=poly.find_closest_cell(targets)
    central=(xyz[:,:,2].mean(axis=1)>=.4*length)&(xyz[:,:,2].mean(axis=1)<=.6*length)
    def summary(indices,weights):
        def rms(value):return float(np.sqrt(np.dot(weights,value**2)/weights.sum())/exact)
        return dict(CFD_magnitude_relative_L2=rms(w[indices]-exact),
                    analytic_nodal_interpolant_magnitude_relative_L2=rms(wi[indices]-exact),
                    CFD_minus_interpolant_tangential_vector_relative_L2=rms(np.linalg.norm(actual[indices]-interpolated[indices],axis=1)))
    result=dict(status='DIAGNOSTIC_ONLY_NO_GATE_CHANGE',case=name,
                fixed_wall_targets=summary(selected,np.ones(len(selected))),
                central_all_facets_area_weighted=summary(central,area[central]),
                interpretation='Raw adjacent-tetra symmetric-gradient traction. All-facet metric is diagnostic, not a replacement for the frozen nearest-facet target WSS gate. Norms do not add. Wall-target traces at surface edges can select one of adjacent discontinuous P1 gradients; no smoothing is applied.')
    a=json.loads((case/'cfd/analysis/CALIBRATION.json').read_text())
    assert abs(result['fixed_wall_targets']['CFD_magnitude_relative_L2']-a['raw_WSS_relative_L2'])<1e-12
    (case/'cfd/analysis/WSS_INTERPOLATION_DIAGNOSTIC.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('case');args=parser.parse_args();main(args.case)
