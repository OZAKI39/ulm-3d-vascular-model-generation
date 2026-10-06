"""Freeze original-data-only coordinates and decisions BEFORE reference payloads."""
from common import *
import datetime

def main():
    target=HERE/'EVENT_FIELD_PROTOCOL.json'
    if target.exists(): raise RuntimeError('PROTOCOL_ALREADY_FROZEN')
    events={}
    release={1189:3.9634867119504023,1329:10.115530289085944,642:1.046333099758476}
    inputs={}
    for pid in PIDS:
        rp=OLD/f'result_{pid}.json'; op=OLD/f'original_{pid}.npz'
        r=read(rp); z=np.load(op); m=np.load(OLD/f'measured_{pid}.npz')
        a=r['accounting'][0]['from_state']; b=r['accounting'][0]['to_state']
        n=unit(a['normal']); tangent=unit(np.asarray(b['x_m'])-a['x_m'])
        if pid!=822: tangent=unit(tangent-n*np.dot(n,tangent))
        else: n=unit(np.eye(3)[np.argmin(abs(tangent))]-tangent*np.dot(tangent,np.eye(3)[np.argmin(abs(tangent))]))
        origin=a['nearest_point_m'] if pid!=822 else a['x_m']
        radius=r['radius_m']; t=z['physical_time_s']
        fast=float(t[np.flatnonzero((t>=release[pid]) & (z['wall_gap_m']>=radius))[0]]) if pid!=822 else b['time_s']
        phases=[float(t[0]),release.get(pid,a['time_s']),fast,float(t[-1])]
        # All generated lengths are determined by original geometry only.
        events[pid]=dict(window_s=r['window_s'],actual_range_s=r['actual_range_s'],
          radius_m=radius,original_track=str(op.relative_to(BASE)),raw_track=f'population_v1/analysis/realization_0000_contactfix/tracks/mb_{pid:06d}.npz',
          accepted_count=len(t),phase_boundaries_s=phases,phase_names=['slow','release','fast'] if pid!=822 else ['birth_to_A','A_to_B','after_B'],
          phase_rule='Original first FACE-to-EDGE release, then first accepted surface gap >= radius; no reference-field dependent selection.',
          A={k:a[k] for k in ['time_s','x_m','tetra','triangle','normal','nearest_point_m','gap_um','u_speed_mm_s']},
          B={k:b[k] for k in ['time_s','x_m','tetra','triangle','normal','nearest_point_m','gap_um','u_speed_mm_s']},
          plane=dict(origin_m=origin,e1=tangent.tolist(),e2=n.tolist(),normal=np.cross(tangent,n).tolist(),
             bounds_rule='Bounds of original A tetra vertices projected to this basis, padded by 10 percent; 161 by 161 fixed grid.',resolution=[161,161]),
          normal_profiles=dict(origins_m=[a['nearest_point_m'],b['nearest_point_m']],normals=[a['normal'],b['normal']],
              distance_m=np.linspace(0,200e-6,401).tolist(),validity='Stop at first non-original FACE/projection or domain departure; retain invalid rows; wall endpoint exact no-slip trace.',
              radius_semantics='d = center-wall distance; h = d-radius only for this recorded sphere radius'),
          CV=dict(seed_cells=[a['tetra']] if pid!=822 else 'all original cells of A-to-B segment',face_neighbor_layers=[0,1,2],
              surface='Oriented exterior triangular faces of each original tetra union; include ALL wall and lateral faces, never interpret as complete channel cross sections.'))
        for p in [rp,op,OLD/f'measured_{pid}.npz',OLD/f'ledger_{pid}.json']:
            inputs[str(p.relative_to(BASE))]=sha(p)
    put('PROTOCOL_SOURCE_SHA256.json',inputs)
    put('EVENT_FIELD_PROTOCOL.json',dict(status='FROZEN_BEFORE_REFERENCE_VELOCITY_PAYLOAD',
      frozen_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),events=events,
      fields=dict(priority=NAMES,alternative_if_ineligible=['M2_medium_frame180','M2_fine_frame180'],max_reference_payloads=3,
        admission='Full physical boundary identity, exact solver/BC/physics, native actual steady field and hashes, node/cell mapping. Do not read additional candidates to improve agreement.'),
      ROI='All original window crossed cells plus TWO face-neighbor layers. No extension planned; visualization only.',
      path_sampling='Enumerate every actual segment in each selected mesh geometrically. Union all face fractions, then same physical midpoints and two Gauss points per subsegment for all fields. Original accepted positions retained separately. Midpoint G is constant on each union span; two-point Gaussian integrates squared vector differences exactly. Scalar speed/direction quadrature refinement compared with four-point rule.',
      weights='Primary arc length per exact union span; secondary original dt multiplied by span fraction. Separate original slow/release/fast. Missing samples retain original denominators.',
      normal_definition='Same original finite triangle closest-point geometry at each coordinate; FACE uses inward original facet normal, EDGE/VERTEX use original finite-feature separation normal and record candidate facet ties/ambiguity. ID822 directional interpretation uses fixed path frame.',
      CV_integration='Nested centroid subdivision of each complete original CV face, levels 0,1,2,3,4,5,6 (1 to 4096 equal-area triangles). GPU P1 at centroids; retain all per-face integrals, successive differences and numerical floor. Original mesh exact face mean separately. If error unresolved do not adjudicate flux trend. Wall face no-slip exact trace.',
      rise822='Same original endpoint levels + 10/90 percent of original scalar increase. Enumerate all crossings and monotonicity; missing/multiple roots explicitly recorded.',
      units=dict(position='m / plots um or mm',velocity='m/s / plots mm/s',G='s^-1',flux='m^3/s'),
      numerical_controls='FP64, strict nonnegative barycentric interior. Boundary trace only on verified exact original boundary. Roundoff allowance used for root reconciliation, never projection/extrapolation. Direction N/A below operation-scaled FP64 resolution; no arbitrary epsilon normalization.',
      judgments='Report actual absolute and event-amplitude-normalized changes; repeated gain/direction/length supports only tested resolutions. Mechanism-altering normal sign/direction or derivative differences -> FIELD_RESOLUTION_SENSITIVE. No universal percentage pass, GCI, Richardson or truth claim.',
      statuses=['INPUT_IDENTITY_FAIL','LOCAL_NUMERICAL_DEFECT_CONFIRMED','FIELD_RESOLUTION_SENSITIVE','FIELD_FEATURE_SUPPORTED_ON_TESTED_RESOLUTIONS','NOT_EVALUABLE_REFERENCE_MISSING','NOT_EVALUABLE_GEOMETRY_OR_BC_CONFOUNDED','NOT_EVALUABLE_NOT_STEADY','PHYSICAL_VALIDATION_NOT_PERFORMED'],
      resources=dict(diagnostic_cumulative_seconds=5400,optional_CFD_seconds=14400,total_seconds=21600,new_CFD_default=0,
          max_tree_RSS_GiB=10,max_nonreclaimable_GiB=34,GPU_max_allocated_GiB=8,cpu_workers=2,BLAS_threads=1,
          actual_preflight='RTX4090 total24564/free15477 MiB, cpu.max768000/100000, cpuset0-15, mem.max45156925440B, current45140172800B largely reclaimable file cache, disk44GB free; recheck each supervised stage.'),
      stop='No new meshes, trajectories, alternative solver, boundary change or further mesh search after controlled comparison. All jobs supervised; exit on limits. Optional CFD only after all section8 conditions, not currently triggered.',
      physical_validation='PHYSICAL_VALIDATION_NOT_PERFORMED',production_readonly=True))
    print('PROTOCOL FROZEN',sha(target))

if __name__=='__main__': main()
