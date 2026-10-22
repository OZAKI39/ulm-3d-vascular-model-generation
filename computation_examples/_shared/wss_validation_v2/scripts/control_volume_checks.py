"""Check measured cap/section flux differences against P1 volume divergence.

Clipped control volumes are diagnostic subsets of actual CFD, not new solutions.
"""
import argparse,csv
import numpy as np
import pyvista as pv
from case_common import *
from vessel_regions import PORTS
from analyze_vessel import csvout

p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);args=p.parse_args()
case=args.case.resolve();m=np.load(case/'SV_MESH/mesh_arrays.npz');f=np.load(case/'frozen_flow/flow_arrays_si.npz')
x=m['points_m'];t=m['tetra'];u=f['velocity_m_s'];qtarget=json.loads((case/'policy.json').read_text())['Q_target_m3_s']
sections={r['section']:r for r in csv.DictReader((case/'reports/internal_sections.csv').open())}
flows={r['boundary']:r for r in csv.DictReader((case/'reports/boundary_flows.csv').open())}
rows=[]
for port in PORTS:
    name=port['name'];origin=np.array(port['real_cut_xyz_um'])*1e-6;n=np.array(port['outward_normal']);L=port['fem_cap']['real_cut_to_cap_axial_um']*1e-6
    lower=origin+.25*L*n;upper=origin+.5*L*n
    axial=(x-origin)@n;radial=np.linalg.norm(x-origin-axial[:,None]*n,axis=1);limit=3*port['radius_um']*1e-6
    select=(axial[t].max(axis=1)>=.25*L)&(radial[t].min(axis=1)<limit)
    tt=t[select];ids,inv=np.unique(tt,return_inverse=True);local_t=inv.reshape(-1,4)
    grid=pv.UnstructuredGrid(np.column_stack([np.full(len(tt),4),local_t]).ravel(),np.full(len(tt),pv.CellType.TETRA,np.uint8),x[ids])
    grid.cell_data['Divergence_s_inv']=np.trace(wss.p1_gradients(x,tt,u),axis1=1,axis2=2)
    for endpoint in ['half','cap']:
        slab=grid.clip(normal=n,origin=lower,invert=False)
        if endpoint=='half':slab=slab.clip(normal=n,origin=upper,invert=True)
        assert slab.n_cells>0
        linked=slab.connectivity();point=origin+(.375 if endpoint=='half' else .625)*L*n
        inside=int(linked.find_containing_cell(point));assert inside>=0,(name,endpoint,'center outside clipped control volume')
        region=linked.cell_data['RegionId'][inside];slab=linked.extract_cells(linked.cell_data['RegionId']==region)
        s=(slab.points-origin)@n;r=np.linalg.norm(slab.points-origin-s[:,None]*n,axis=1)
        assert r.max()<limit and s.min()>=.25*L-1e-12
        if endpoint=='half':assert s.max()<=.5*L+1e-12
        sized=slab.compute_cell_sizes(length=False,area=False,volume=True)
        vol=np.asarray(sized.cell_data['Volume']);assert np.all(vol>=0)
        div=np.asarray(sized.cell_data['Divergence_s_inv']);integral=float(np.dot(vol,div))
        qlo=float(sections[name+'_extension_0.25']['Q_m3_s'])
        role='INLET' if name=='INLET' else 'OUTLET_0'+name[1]
        qhi=float(sections[name+'_extension_0.5']['Q_m3_s']) if endpoint=='half' else float(flows[role]['Q_outward_m3_s'])
        difference=qhi-qlo;residual=integral-difference
        rows.append(dict(case=case.name,region=name+'_extension_0.25_to_'+endpoint,unit_flux='m3/s',volume_um3=float(vol.sum()*1e18),clipped_cells=len(vol),upper_or_cap_Q_m3_s=qhi,lower_section_Q_m3_s=qlo,flux_difference_m3_s=difference,volume_integral_divergence_m3_s=integral,divergence_theorem_residual_m3_s=residual,residual_relative_to_inlet=abs(residual)/qtarget,flux_difference_pct_of_inlet=100*difference/qtarget,complete_lumen_radial_guard=True,weight='exact_P1_gradient_times_clipped_cell_volume;P1_section_and_cap_flux',reference='same_actual_CFD;diagnostic_identity_not_pointwise_mass_conservation'))
csvout(case/'reports/control_volume_flux.csv',rows)
allrows=[]
for path in sorted(V.glob('stage[34]/*/reports/control_volume_flux.csv')):allrows.extend(list(csv.DictReader(path.open())))
csvout(V/'data/control_volume_flux.csv',allrows)
maximum=max(r['residual_relative_to_inlet'] for r in rows)
print(case.name,'maximum divergence-theorem discrepancy / inlet',maximum)
assert maximum<1e-7,'Diagnostic clipping/flux identity needs review; not accepted as a flow conservation error'
