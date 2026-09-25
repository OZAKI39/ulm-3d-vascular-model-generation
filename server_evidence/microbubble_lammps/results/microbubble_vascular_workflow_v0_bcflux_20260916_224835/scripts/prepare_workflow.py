#!/usr/bin/env python3
"""Derived geometry bounds and immutable size-stream input. No flow reconstruction."""
from pathlib import Path
import sys,json,heapq,hashlib
import numpy as np
import yaml
from flow_geometry import Geometry,FrozenSampler
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S/'inputs'))
from sonovue_sampler import SonoVueDistribution
g=Geometry(S/'geometry/GEOMETRY_ARRAYS.npz');field=FrozenSampler(S/'fields/FROZEN_FLOW_FIELD_V0.h5')
raw=(S/'geometry/CLOSED_REFERENCE.stl').read_bytes();records=np.frombuffer(raw,offset=84,dtype=np.dtype('V50'))
for region in range(3):
    subset=records[g.g['classes']==region]
    (S/f'geometry/REGION_{region}.stl').write_bytes(raw[:80]+np.array([len(subset)],'<u4').tobytes()+subset.tobytes())
section=dict(np.load(S/'geometry/INJECTION_SECTION.npz'));tri=section['triangles'];normal=section['normal']
np.savetxt(S/'geometry/injection_triangles.txt',tri.reshape(-1,9),fmt='%.17g')
manifest=json.loads((S/'geometry/BOUNDARY_MANIFEST.json').read_text())
with (S/'geometry/ports.txt').open('w') as f:
    for p in manifest['ports']:f.write(p['name']+' '+' '.join(format(x,'.17g') for x in p['center_m']+p['outward_unit_normal'])+' '+p['stl']+'\n')
# Distance to a fixed wall set is 1-Lipschitz. Every point in a triangle is
# within max vertex distance of its centroid. Heap upper bounds certify
# the global maximum, not just a local optimizer or half-width heuristic.
heap=[];counter=0;best=0.;witness=None
def add(t):
    global counter,best,witness
    c=t.mean(axis=0);d=g.distance(c)[0];status,u=field.query([c])
    if status[0]==0 and -u[0]@normal>0 and d>best:best=d;witness=c.copy()
    bound=d+np.linalg.norm(t-c,axis=1).max()+2e-14
    heapq.heappush(heap,(-bound,counter,t));counter+=1
for t in tri:add(t)
while -heap[0][0]-best>1e-11:
    _,_,t=heapq.heappop(heap);a,b,c=t;ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
    for sub in [[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]:add(np.array(sub))
    assert counter<1000000,'GEOMETRY_CERTIFICATION_BUDGET'
upper=-heap[0][0];clearance=field.dx
proof=dict(status='PASS',method='continuous triangle branch-and-bound using 1-Lipschitz wall distance',wall='WALL_ONLY.stl',cap_excluded=True,lower_max_distance_m=best,upper_max_distance_m=upper,witness_m=witness.tolist(),certification_tolerance_m=1e-11,triangles_evaluated=counter,clearance_m=clearance,admissible_radius_lower_m=best-clearance,admissible_radius_upper_m=upper-clearance,uncertain_band_policy='stop explicitly, never classify by finite failed random search')
(S/'validation/INLET_SIZE_CAPACITY.json').write_text(json.dumps(proof,indent=2)+'\n')
dist=SonoVueDistribution(S/'inputs/FROZEN_SONOVUE_HISTOGRAM.csv');seed=42;u=np.random.default_rng(seed).random(20000);diam=dist.inverse_cdf(u)
np.savetxt(S/'inputs/size_stream.txt',np.c_[u,diam*.5e-6],fmt='%.17g',header='NOT EXPERIMENTAL CONCENTRATION; frozen inverse CDF; uniform_u radius_m')
(S/'inputs/SIZE_STREAM_PROVENANCE.json').write_text(json.dumps(dict(seed=seed,draws_precomputed=len(u),consumption='rank0 only; index advances only on committed transaction; never resample a pending candidate',sampler_sha256=hashlib.sha256((S/'inputs/sonovue_sampler.py').read_bytes()).hexdigest(),histogram_sha256=dist.histogram_sha256,diameter_to_radius_m='diameter_um * 0.5e-6',stream_sha256=hashlib.sha256((S/'inputs/size_stream.txt').read_bytes()).hexdigest()),indent=2)+'\n')
Q=json.loads((S/'validation/AUTHORITATIVE_INLET_FLOW_AUDIT.json').read_text())['authoritative_flow_rate'];speed=Q/float(section['area'])
typical_a=float(dist.inverse_cdf(.5))*.5e-6;clear_time=2*(typical_a+clearance)/speed
xyz=g.g['points'];length=float(np.linalg.norm(xyz.max(0)-xyz.min(0)));transit=length/speed
contract=dict(disclaimer='NOT EXPERIMENTAL CONCENTRATION',flux_mode='TEST_ONLY_DIRECT_NUMBER_FLUX',production_bubble_concentration_status='UNSPECIFIED',authoritative_Q_in_m3_s=Q,injection_position_policy='VALID_SAMPLER_SUPPORT_ONLY',position_distribution_status='WORKFLOW_V0_APPROXIMATION',default_position_mode='VALID_AREA_UNIFORM',alternative_position_mode='VALID_SUPPORT_FLUX_WEIGHTED',typical_radius_m=typical_a,nominal_Q_over_section_area_m_s=speed,typical_inlet_clearance_time_s=clear_time,geometry_extent_time_scale_s=transit,flux_choice_basis='LOW 0.12 / clearing time; MEDIUM 0.65 / clearing time; STRESS 8 / clearing time; bounding-box extent/speed is only a nominal transit estimate, not actual residence time',cases={
 'LOW':dict(number_flux_per_s=.12/clear_time,max_time_s=.25,control_basis='SOURCE_POPULATION'),
 'MEDIUM':dict(number_flux_per_s=.65/clear_time,max_time_s=.15,control_basis='ADMITTED_POPULATION'),
 'CAPACITY_STRESS':dict(number_flux_per_s=8/clear_time,max_time_s=.08,control_basis='ADMITTED_POPULATION'),
 'FLUX_WEIGHTED':dict(number_flux_per_s=.12/clear_time,max_time_s=.12,control_basis='SOURCE_POPULATION',position_mode='VALID_SUPPORT_FLUX_WEIGHTED')},dt_max_s=1e-4,dt_min_s=1e-10,wall_margin_m=1e-10,injection_wall_clearance_m=clearance,injection_pair_clearance_m=2e-9,pending_max_count=20,pending_max_age_s=.04,max_source_draws_per_step=1000,position_attempts_per_candidate=4096,max_retry_count=24,max_steps=100000,nearwall_gap_over_radius_threshold=.2,nearwall_min_continuous_residence_s=.001)
contract['cases']['CAPACITY_STRESS_LIMIT']=dict(number_flux_per_s=16/clear_time,max_time_s=.08,control_basis='ADMITTED_POPULATION')
contract['stress_followup_rationale']='First 8/clearing-time stress reaches wall stall with 17 pending before the unchanged 20 pending gate; 16/clearing-time stress explicitly exercises capacity gate without changing geometry or safety.'
contract['rolling_window_duration_s']=.05
(S/'configs/test_flux_cases.yaml').write_text(yaml.safe_dump(contract,sort_keys=False))
remote='/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_'+S.name
for name,case in contract['cases'].items():
    cfg=dict(stage=remote,disclaimer='NOT_EXPERIMENTAL_CONCENTRATION',flux_mode=contract['flux_mode'],authoritative_Q=Q,number_flux=case['number_flux_per_s'],control_basis=case['control_basis'],position_mode=case.get('position_mode',contract['default_position_mode']),max_time=case['max_time_s'],dt_max=contract['dt_max_s'],dt_min=contract['dt_min_s'],wall_margin=contract['wall_margin_m'],injection_clearance=contract['injection_wall_clearance_m'],pair_clearance=contract['injection_pair_clearance_m'],size_lower=proof['admissible_radius_lower_m'],size_upper=proof['admissible_radius_upper_m'],steps_cap=contract['max_steps'],max_retry=contract['max_retry_count'],draw_budget=contract['max_source_draws_per_step'],position_attempts=contract['position_attempts_per_candidate'],max_pending=contract['pending_max_count'],pending_max_age=contract['pending_max_age_s'],nearwall_threshold=contract['nearwall_gap_over_radius_threshold'])
    (S/'configs'/f'{name}.cfg').write_text('\n'.join(f'{k} {v}' for k,v in cfg.items())+'\n# NOT EXPERIMENTAL CONCENTRATION\n')
print(json.dumps(dict(size_capacity=proof,test_flux_cases=contract),indent=2))
