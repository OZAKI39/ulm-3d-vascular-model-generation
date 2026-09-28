"""Differentiate internal plateaus, screen projection and terminal-frame hold."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parents[1];RUN=HERE.parent

def longest(mask,weights):
    best=now=0.
    for on,w in zip(mask,weights):
        now=now+float(w) if on else 0.;best=max(best,now)
    return best

def main():
    rows=json.loads((RUN/'data/metrics.json').read_text());m=json.loads((RUN/'data/render_manifest.json').read_text())
    ages=np.linspace(0,m['age_end_s'],m['video_frames']);layout=m['animation_layout']
    d=np.array(layout['camera_position'])-layout['camera_focal_point'];d/=np.linalg.norm(d)
    right=np.cross([0,0,1],d);right/=np.linalg.norm(right);up=np.cross(d,right)
    scale=1080/(2*layout['parallel_scale']);output=[]
    for r in rows:
        p=RUN/'tracks'/f"mb_{r['particle_id']:06d}"/'trajectory.npz'
        with p.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==r['trajectory_file_sha256']
        a=np.load(p)['samples'];dt=np.diff(a[:,0]);ds=np.linalg.norm(np.diff(a[:,1:4],axis=0),axis=1)
        speed=ds/dt;stationary=longest(speed<=1e-9,dt)*1e3;exact=longest(ds==0,dt)*1e3
        xyz=np.column_stack([np.interp(ages,a[:,0],a[:,j]) for j in [1,2,3]])*1e6
        projected=xyz@np.column_stack([right,up])*scale
        dpix=np.linalg.norm(np.diff(projected,axis=0),axis=1)
        active=ages[1:]<=a[-1,0]
        low=(dpix<.5)&active
        runframes=int(round(longest(low,np.ones(len(low)))))
        retained=max(0.,m['age_end_s']-a[-1,0])/m['age_end_s']*m['video_seconds'] if r['status']=='SUPPORTED_STATIONARY' else 0.
        output.append(dict(particle_id=r['particle_id'],status=r['status'],max_exact_center_plateau_ms=exact,
            max_displacement_speed_le_1e_9_m_s_ms=stationary,maximum_subpixel_active_intervals=runframes,
            subpixel_max_video_duration_s=runframes/m['video_fps'],post_termination_marker_hold_video_s=retained))
    with (HERE/'data/apparent_pauses.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    finished=[r for r in output if r['status']=='COMPLETED'];stops=[r for r in output if r['status']=='SUPPORTED_STATIONARY']
    result=dict(all_pass=True,physical_frame_spacing_ms=float(np.diff(ages)[0]*1e3),
        diagnostic_zero_speed_m_s=1e-9,diagnostic_subpixel_distance_px=.5,
        completed_tracks_with_near_zero_plateau_at_least_0p5ms=sum(r['max_displacement_speed_le_1e_9_m_s_ms']>=.5 for r in finished),
        completed_max_near_zero_plateau_ms=max(r['max_displacement_speed_le_1e_9_m_s_ms'] for r in finished),
        completed_max_exact_plateau_ms=max(r['max_exact_center_plateau_ms'] for r in finished),
        completed_with_at_least_3_subpixel_intervals=sum(r['maximum_subpixel_active_intervals']>=3 for r in finished),
        completed_max_subpixel_video_s=max(r['subpixel_max_video_duration_s'] for r in finished),
        stopped_terminal_marker_hold_video_s=[min(r['post_termination_marker_hold_video_s'] for r in stops),max(r['post_termination_marker_hold_video_s'] for r in stops)],
        scope='Subpixel diagnostic is specific to the main movie camera and threshold, not a physical stop classification. No frame interpolation beyond the original endpoint is claimed as physics.')
    (HERE/'data/apparent_pauses_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
