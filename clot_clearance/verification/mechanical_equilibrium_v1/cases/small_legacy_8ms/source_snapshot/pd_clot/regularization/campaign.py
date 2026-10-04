"""Reproducible, non-overwriting regularization studies; no appearance fitting."""
import argparse,copy,hashlib,json
from pathlib import Path
from .coupon import calibrate
from .runner import run
from .audit import audit_run

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'configs/streaming_regularized_demo_coarse.json'


def run_study(study, output_root, sweep_parameter=None):
    out=Path(output_root);out.mkdir(parents=True,exist_ok=True)
    base=json.loads(BASE.read_text());base_cal=ROOT/'verification/regularization/calibration_coarse/CALIBRATION.json'

    def case(name,c,calibration=None):
        dest=out/name
        if (dest/'SUMMARY.json').exists():
            if json.loads((dest/'CONFIG.json').read_text())!=c:raise ValueError(f'Existing config differs: {dest}')
            audit_run(dest)
            print('REUSE_COMPLETED',str(dest),flush=True);return
        if dest.exists():raise FileExistsError(f'Unfinished run preserved; use a new output root: {dest}')
        if calibration is None:
            cp=out/'calibration'/name/'CALIBRATION.json'
            if not cp.exists():calibrate(c,c['clot']['particle_spacing_m'],c['regularized_damage']['Gc_demo_J_m2'],cp.parent)
            calibration=cp
        print('START_CASE',name,flush=True)
        result=run(c,dest,calibration)
        audit_run(dest)
        print('END_CASE',name,result['elapsed_wall_s'],flush=True)

    if study in ['demo','all']:case('coarse',copy.deepcopy(base),base_cal)
    if study in ['cycle','all']:
        for delta in [500,250]:
            c=copy.deepcopy(base);c['damage']['DeltaN']=delta;c['simulation']['number_of_macro_steps']=25000//delta
            case(f'cycle_{delta}',c,base_cal)
    if study in ['mesh','all']:
        for name,scale in [('medium',1.5),('fine',2.)]:
            c=copy.deepcopy(base);c['clot']['particle_spacing_m']/=scale;c['clot']['cells']=[int(v*scale) for v in base['clot']['cells']]
            case(name,c)
    if study in ['timestep','all']:
        c=copy.deepcopy(base);c['simulation']['dt_s']/=2
        # Full exposure avoids silently choosing a window without a failure event.
        case('dt_half',c,base_cal)
    if study in ['sweep','all']:
        parameters=[('Gc_demo_J_m2',[.005,.02]),('m_E',[.8,1.2]),('C_E',[.0005,.002]),('streaming_velocity_scale_m_s',[.005,.02])]
        for name,values in parameters:
            if sweep_parameter is not None and sweep_parameter!=name:continue
            for value in values:
                c=copy.deepcopy(base)
                if name=='streaming_velocity_scale_m_s':c['streaming'][name]=value;cal=base_cal
                else:c['regularized_damage'][name]=value;cal=None
                case(f'sweep_{name}_{value:g}',c,cal)
    if study in ['sweep_surface','all']:
        for Gc in [.005,.02]:
            if sweep_parameter is not None and sweep_parameter!=str(Gc):continue
            for U in [.005,.02]:
                c=copy.deepcopy(base);c['regularized_damage']['Gc_demo_J_m2']=Gc;c['streaming']['streaming_velocity_scale_m_s']=U
                cp=out/'calibration'/f'sweep_Gc_demo_J_m2_{Gc:g}'/'CALIBRATION.json'
                if not cp.exists():raise FileNotFoundError('Run the Gc sensitivity calibration first')
                case(f'surface_Gc_{Gc:g}_U_{U:g}',c,cp)


def main():
    p=argparse.ArgumentParser();p.add_argument('--study',choices=['demo','cycle','mesh','timestep','sweep','sweep_surface','all'],required=True)
    p.add_argument('--output-root',type=Path,default=ROOT/'results/streaming_regularized_demo')
    p.add_argument('--sweep-parameter',help='Optional independent sweep subset; no material fitting')
    a=p.parse_args();run_study(a.study,a.output_root,a.sweep_parameter)


if __name__=='__main__':main()
