"""Standalone CPU figures of actual authorized observations; no solver imports."""
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from py_scripts.fluid_physics.common import now,sha256_file,write_json
from py_scripts.single_rbc_benchmark.analysis import read_csv
from py_scripts.single_rbc_benchmark.physics import read_off
from .workflow import load_config,paths


def main():
    start=time.perf_counter();c=load_config();b,runs,out=paths(c);folder=out/'authorized_figures';folder.mkdir(exist_ok=False)
    inputs=[];fig,axes=plt.subplots(2,2,figsize=(11,6),layout='constrained')
    for task,label in [('A4_continuous_native_observation','A4 dt=0.001, planned t=5'),('A5_continuous_preparation_30','A5 dt=0.001, failed at step 6407'),('A6_continuous_preparation_half_dt','A6 dt=0.0005, completed t=30')]:
        p=b/'runtime_reviews'/task/'review/frame_metrics.csv';inputs.append(p);a=read_csv(p)
        spec=runs/'gpu'/task/'spec.json';inputs.append(spec);dt=json.loads(spec.read_text())['dt'];t=a['step']*dt
        for ax,key,mult,title in [(axes[0,0],'D',1,'D'),(axes[0,1],'theta_deg',1,'XZ angle [deg]'),(axes[1,0],'area_relative_drift',100,'|A/A0 - 1| [%]'),(axes[1,1],'max_wlc_extension_fraction',1,'Maximum WLC extension fraction')]:
            ax.plot(t,a[key]*mult,label=label,lw=1.2);ax.set(xlabel='Zero-shear elapsed t*',ylabel=title);ax.grid(alpha=.2)
    axes[1,0].axhline(2,color='red',ls='--',lw=1);axes[1,1].axhline(1,color='red',ls='--',lw=1)
    axes[0,0].legend(fontsize=8);fig.suptitle('Stored preparation states only; Gamma = 0 for every curve')
    first=folder/'stored_preparation.svg';fig.savefig(first);plt.close(fig)
    ref,faces=read_off(b/'common_reference.off');inputs.append(b/'common_reference.off')
    p=runs/'gpu/A5_continuous_preparation_30/bounce_relaxation_outer_local_6407_membrane.csv';inputs.append(p);raw=read_csv(p);order=np.argsort(raw['id']);assert np.array_equal(raw['id'][order],np.arange(len(ref)))
    old=np.column_stack([raw[k] for k in ('old_x','old_y','old_z')])[order];new=np.column_stack([raw[k] for k in ('x','y','z')])[order]
    edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,points,title,limits in [(axes[0],old,'Step 6407: old state; WLC max = 7.486',(-9,9)),(axes[1],new,'Step 6407: before bounce; WLC max = 88.425',(-70,70))]:
        ax.add_collection(LineCollection(ref[edges][:,:,[0,2]],colors='#6e8995',linewidths=.4,alpha=.5,label='Reference'))
        ax.add_collection(LineCollection(points[edges][:,:,[0,2]],colors='#af3532',linewidths=.45,alpha=.7,label='Actual native mesh'))
        ax.set(xlim=limits,ylim=limits,xlabel='rank-local x',ylabel='rank-local z',title=title);ax.set_aspect('equal');ax.grid(alpha=.2)
        ax.axhline(12,ls='--',c='#555',lw=.7);ax.axhline(-12,ls='--',c='#555',lw=.7)
    axes[0].legend(fontsize=8);fig.suptitle('Actual pre-bounce coordinates; no folding or rotation to hide the failure')
    second=folder/'failure_geometry.svg';fig.savefig(second);plt.close(fig)
    record=dict(recorded_at=now(),scope='CPU figures from saved states; not new simulation or qualified endpoint comparison',
                source_sha256={str(p):sha256_file(p) for p in inputs},figures={str(p):sha256_file(p) for p in (first,second)},cpu_render_s=time.perf_counter()-start)
    write_json(b/'authorized_figures.json',record);print(json.dumps(record,indent=2))


if __name__=='__main__':main()
