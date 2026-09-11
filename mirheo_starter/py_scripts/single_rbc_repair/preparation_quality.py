"""CPU preparation-window evaluation with the predeclared unchanged criteria."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from py_scripts.fluid_physics.common import now,sha256_file,write_json
from py_scripts.single_rbc_benchmark.analysis import read_csv
from py_scripts.single_rbc_benchmark.physics import read_off,geometry
from .geometry_checks import checked_mesh


def inspect(derived, criteria_file, output):
    start=time.perf_counter();d=Path(derived).resolve();spec=json.loads((d/'spec.json').read_text());c=spec['config']
    criteria=json.loads(Path(criteria_file).read_text())['common_preparation_criteria']
    ref,faces=read_off(spec['mesh']);scale=geometry(ref,faces)['a'];ref-=ref.mean(0)
    raw=read_csv(d/'vertices.csv');frames=[]
    for step in sorted(set(raw['step'])):
        rows=raw[raw['step']==step];v,g=checked_mesh(rows['vertex'],np.column_stack([rows[k] for k in ('x','y','z')]),faces,c['geometry']['periodic_length'])
        frames.append(dict(step=int(step),elapsed_time_star=int(step)*spec['dt'],centered=v-v.mean(0),theta=g['theta_deg']))
    width=criteria['window_time_star'];windows=[]
    final=frames[-1]['elapsed_time_star']
    for index in range(int((final+1e-9)//width)):
        lo=index*width;hi=lo+width;selected=[f for f in frames if lo-1e-9<=f['elapsed_time_star']<hi-1e-9]
        if not selected:continue
        mean=np.mean([f['centered'] for f in selected],axis=0)
        angle=float(np.degrees(.5*np.angle(np.mean(np.exp(2j*np.radians([f['theta'] for f in selected]))))))
        windows.append(dict(begin=lo,end=hi,samples=len(selected),reference_shape_rms_over_a=float(np.sqrt(np.mean(np.sum((mean-ref)**2,axis=1)))/scale),mean_coordinates=mean,mean_angle_deg=angle))
    transitions=[]
    for before,after in zip(windows,windows[1:]):
        residual=float(np.sqrt(np.mean(np.sum((after['mean_coordinates']-before['mean_coordinates'])**2,axis=1)))/scale/width)
        angle=abs((after['mean_angle_deg']-before['mean_angle_deg']+90)%180-90)/width
        transitions.append(dict(before_window=[before['begin'],before['end']],after_window=[after['begin'],after['end']],residual_shape_rms_over_a_per_time=residual,orientation_change_degrees_per_time=angle,
            residual_pass=residual<=criteria['residual_shape_rms_over_a_per_time'],orientation_pass=angle<=criteria['orientation_change_degrees_per_time']))
    required=criteria['minimum_consecutive_windows'];last=transitions[-required:]
    ready=len(last)==required and all(x['residual_pass'] and x['orientation_pass'] for x in last)
    for w in windows:w.pop('mean_coordinates')
    result=dict(recorded_at=now(),source_directory=str(d),criteria=criteria,
        definition='Center translation removed per frame; nonoverlapping 2-unit windows average corresponding vertex coordinates, with no rotation or time fitting. Residual is RMS between consecutive window means divided by reference a and window duration. Angles use the same XZ geometry definition and axial circular means.',
        reference_shape_scope='Distance to the shared reference is diagnostic only, not a newly measured HemoCell-versus-Mirheo acceptance test.',
        final_saved_elapsed_time_star=final,final_reference_shape_rms_over_a=float(np.sqrt(np.mean(np.sum((frames[-1]['centered']-ref)**2,axis=1)))/scale),
        windows=windows,transitions=transitions,last_required_transitions_pass=ready,
        status='RESIDUAL_WINDOWS_PASS_OTHER_GATES_STILL_REQUIRED' if ready else 'PREPARATION_RESIDUAL_CRITERIA_FAILED_OR_INCOMPLETE',
        new_HemoCell_match_measured=False,qualified_speedup=None,
        source_sha256={str(p):sha256_file(p) for p in (d/'spec.json',d/'vertices.csv',Path(criteria_file).resolve(),Path(spec['mesh']))},cpu_analysis_s=time.perf_counter()-start)
    write_json(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--derived',required=True);p.add_argument('--criteria',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    r=inspect(a.derived,a.criteria,a.output);print(json.dumps({k:v for k,v in r.items() if k not in ('source_sha256','windows','transitions')},ensure_ascii=False,indent=2))
