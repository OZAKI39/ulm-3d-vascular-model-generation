"""Executable A-F acceptance report; small controls, no BraVa computations."""
import argparse,copy,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pd_clot.runner import simulate
from pd_clot.geometry import make_cloud
from pd_clot.mechanics import evaluate_cloud,prepare_shape
from pd_clot.damage import instantaneous,fragments,particle_damage
from pd_clot.output import write_json,write_vtk,write_collection


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reuse-completed',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parents[1];output=root/'verification';output.mkdir(exist_ok=True)
    c=json.loads((root/'configs/straight_pipe.json').read_text())
    c['clot']['cells']=[6,4,4];c['clot']['origin_m']=[-.000375,-.00025,-.0011]
    c['simulation'].update(number_of_macro_steps=3,representative_frequency_Hz=1000.,warmup_cycles=0)
    c['streaming'].update(bubble_center_m=[0,0,-.0005],include_pipe_traction=False)
    checks={};runs={}
    def run(name,conf):
        path=output/name
        if args.reuse_completed and (path/'SUMMARY.json').exists():
            if json.loads((path/'CONFIG.json').read_text())!=conf:raise ValueError('Stored control configuration differs')
            z=np.load(path/'states.npz')
            r=dict(summary=json.loads((path/'SUMMARY.json').read_text()),history=json.loads((path/'history.json').read_text()),
                cloud=make_cloud(conf['clot']),x=z['x'][-1],integrity=z['integrity'][-1],states=[dict(damage=z['damage'][-1])])
        else:r=simulate(conf,path,quiet=True)
        runs[name]=r;return r
    a=copy.deepcopy(c);a['streaming']['traction_scale_Pa']=0;r=run('A_no_load',a)
    checks['A_no_load']=r['summary']['final']['maximum_displacement_um']<1e-7 and r['integrity'].min()==1 and r['summary']['final']['number_of_fragments']==1
    a=copy.deepcopy(c);a['streaming'].update(traction_scale_Pa=1.,uniform=True);r=run('B_small_uniform',a)
    checks['B_small_uniform']=r['summary']['final']['maximum_displacement_um']>0 and r['integrity'].min()==1
    a=copy.deepcopy(c);a['streaming']['traction_scale_Pa']=160.;a['damage'].update(Q_ref=.004,C_damage=.000002);r=run('C_localized',a)
    g=r['cloud'];D=r['states'][-1]['damage'];upper=g.X[:,2]>g.X[:,2].mean()
    checks['C_localization']=D[upper].mean()>D[g.fixed].mean() and D.max()>0
    checks['D_history_monotone']=bool(np.all(np.diff([z['mean_particle_damage'] for z in r['history']])>=0))
    stronger=copy.deepcopy(a);stronger['streaming']['traction_scale_Pa']=180.;rs=run('C_stronger_control',stronger)
    checks['stronger_load_control']=rs['history'][-1]['mean_particle_damage']>=r['history'][-1]['mean_particle_damage']
    a['damage']['C_damage']=0;r=run('E_zero_rate',a)
    checks['E_zero_rate']=bool(np.all(r['integrity']==1))
    a['damage']['enabled']=False;r=run('damage_disabled',a)
    checks['damage_disabled']=bool(np.all(r['integrity']==1))
    # F is a deliberately prescribed cleavage displacement; no claim of flow-caused fragmentation.
    g=make_cloud(c['clot']);x=g.X.copy();x[g.X[:,2]>g.X[:,2].mean(),2]+=5*g.spacing
    stretch=np.linalg.norm(x[g.pairs[:,1]]-x[g.pairs[:,0]],axis=1)/g.length
    integrity=instantaneous(np.ones(len(g.pairs)),stretch,1.01,1.02)
    shape=prepare_shape(g,integrity,c['safety']);mech=evaluate_cloud(g,x,integrity,c['material'],c['safety'],shape)
    labels,degree,frag=fragments(g,integrity,0);D,Dmax=particle_damage(g,integrity)
    fdir=output/'F_artificial_cleavage'
    fsummary=dict(kind='ARTIFICIAL_KINEMATIC_SOFTWARE_TEST',opening_m=5*g.spacing,s1=1.01,s2=1.02,**frag)
    if args.reuse_completed and (fdir/'SUMMARY.json').exists():
        if json.loads((fdir/'SUMMARY.json').read_text())!=fsummary:raise ValueError('Cleavage evidence differs')
    else:
        fdir.mkdir(exist_ok=False)
        write_vtk(fdir/'vtk',0,g,x,np.zeros_like(x),integrity,dict(particle_damage=D,max_bond_damage=Dmax,
            number_of_active_bonds=degree,strain_energy_density=mech[3],fragment_id=labels),np.zeros_like(g.face_position))
        write_collection(fdir,[(0,0)])
        write_json(fdir/'SUMMARY.json',fsummary)
    checks['F_fragmentation']=frag['number_of_fragments']==2
    checks['all_fixed_bases_unchanged']=all(r['history'][-1]['fixed_base_error_m']==0 for r in runs.values())
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,
        results={name:r['summary']['final'] for name,r in runs.items()},artificial_fragmentation=frag,
        limitations='A-F software verification. Artificial cleavage is imposed kinematics; not ultrasound fragmentation.')
    write_json(output/'ACCEPTANCE.json',result);print(json.dumps(result,indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
