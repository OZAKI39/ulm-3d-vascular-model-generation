"""Independent preservation checks and a small exact-geometry cross-check."""
from pathlib import Path
import json,hashlib,sys,xml.etree.ElementTree as ET
import numpy as np
from audit_math import align_saved,wall_coefficient
from run_audit import init,wall_sample,ENV,dump
ROOT=Path(__file__).resolve().parents[1]

def main():
    baseline=json.loads((ROOT/'data/readonly_baseline.json').read_text());changed=[]
    for p,h in baseline.items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h:changed.append(p)
    assert not changed,changed
    cfg=json.loads((ROOT/'data/local_config.json').read_text());init(cfg)
    from particle_3d.wall_gap import wall_gap
    from particle_3d.particle_shapes import Sphere
    from particle_3d.nearfield_regularization import NearFieldRegularizationV1
    from particle_3d.hydrodynamic_resistance import PhysicalNearField
    inv=json.loads((ROOT/'data/trajectory_inventory.json').read_text());chosen=[]
    for r in inv['records']:
        if r['sample_count']==0:continue
        m=json.loads(Path(r['local_path']).with_suffix('.json').read_text())
        if m.get('nearfield_states',{}).get('CONTINUUM_HANDOFF_CONTACT',0)>3:chosen.append(r)
        if len(chosen)==2:break
    rows=[]
    for r in chosen:
        with np.load(r['local_path']) as f:s=f['samples'].copy()
        x=align_saved(s);idx=np.unique(np.r_[np.linspace(0,len(s)-1,16).astype(int),np.argmin(s[:,14])])
        distance,n,tri=wall_sample(x[idx]);a=r['radius_m']
        for j,i in enumerate(idx):
            sphere=Sphere(x[i],a);g=wall_gap(sphere,ENV['wall'])
            z=NearFieldRegularizationV1().evaluate(PhysicalNearField(1,None,g.gap_m,g.normal_inward,g.roundoff_m),{1:sphere},ENV['mu'])[0]
            h=distance[j]-a
            assert abs(h-g.gap_m)<2e-17
            assert np.linalg.norm(n[j]-g.normal_inward)<1e-10
            np.testing.assert_allclose(wall_coefficient(a,h,ENV['mu']),z,rtol=1e-7,atol=1e-30)
            rows.append(dict(dataset=r['dataset'],id=r['id'],state=int(s[i,17]),index=int(i),
                             gap_error_m=float(abs(h-g.gap_m)),normal_error=float(np.linalg.norm(n[j]-g.normal_inward))))
    newcase=Path(cfg['new_flow']).parents[1];xml=ET.parse(newcase/'run/solver.xml')
    rho=float(xml.findtext('.//Add_equation/Density'));mu=float(xml.findtext('.//Viscosity/Value'))
    assert rho==ENV['rho'] and mu==ENV['mu']
    new_manifest=json.loads((newcase/'frozen_flow/manifest.json').read_text())
    physics=json.loads((newcase/'reports/physics_validation.json').read_text())
    dump(ROOT/'data/independent_verification.json',dict(status='PASS',protected_file_count=len(baseline),changed_files=changed,
        geometry_parity_samples=len(rows),geometry_rows=rows,new_solver_density=rho,new_solver_viscosity=mu,
        new_native_manifest=new_manifest,new_physics_validation=physics,formal_solver_modified=False,trajectory_recomputed=False))
    print('READONLY_AND_GEOMETRY_PASS',len(baseline),len(rows))

if __name__=='__main__':main()
