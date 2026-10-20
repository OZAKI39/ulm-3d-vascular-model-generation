"""Two native frictionless spherical contacts; independent work and reaction."""
from common import *
single=module('accepted_single_indenter',HA/'indenter.py')
APPROACH=single.APPROACH;PEAK=single.PEAK;RETRACT=single.RETRACT;CLEAR=single.CLEAR
GAP=single.GAP;DEPTH=single.DEPTH

class Actuator(single.Indenter):
    def __init__(self,ref,config,delay):
        self.n=int(ref['clot_count']);self.N=len(ref['X']);self.config=config
        self.c0=np.asarray(config['initial_center_m']);self.radius=config['radius_m'];self.k=config['penalty_K_N_m2'];self.direction=np.array(config['inward_normal']);self.delay=delay;self.id=config['fix_id']
    def center(self,t):return self.c0+self.direction*single.displacement(t-self.delay)
    def velocity(self,t):return self.direction*single.speed(t-self.delay)
    def command(self):
        def expr(a,b):
            a+=self.delay;b+=self.delay
            return f'((time>={b:.17g})+(time>{a:.17g})*(time<{b:.17g})*0.5*(1-cos(PI*(time-{a:.17g})/{b-a:.17g})))'
        displacement=f'({GAP:.17g}*{expr(0,APPROACH)}+{DEPTH:.17g}*{expr(APPROACH,PEAK)}-{GAP+DEPTH:.17g}*{expr(RETRACT,CLEAR)})'
        commands=[];position=[]
        for axis,c,direction in zip('xyz',self.c0,self.direction):
            if direction:
                name=self.id+'_'+axis;commands.append(f'variable {name} equal "{c:.17g}+{direction:.17g}*{displacement}"');position.append('v_'+name)
            else:position.append(f'{c:.17g}')
        commands.append(f'fix {self.id} clot indent {self.k:.17g} sphere '+ ' '.join(position)+f' {self.radius:.17g} units box')
        commands.append(f'fix_modify {self.id} energy no')
        return '\n'.join(commands)

class DistributedLoading:
    def __init__(self,ref,case):
        self.spec=json.loads((O/'LOAD_SITE_SELECTION.json').read_text());self.protocol=json.loads((O/'MULTISITE_LOADING_PROTOCOL.json').read_text())['cases'][case]
        self.sites=[Actuator(ref,config,delay) for config,delay in zip(self.spec['sites'],self.protocol['delays_s'])]
        self.clear_time=max(q.delay+CLEAR for q in self.sites)
        self.boundaries=sorted(set([0.,END,.0028]+[delay+t for delay in self.protocol['delays_s'] for t in [0,APPROACH,PEAK,RETRACT,CLEAR]]))
    def command(self):return '\n'.join(q.command() for q in self.sites)
    def evaluate(self,x,t):
        sites=[q.evaluate(x,t) for q in self.sites];forces=np.array([q['force'] for q in sites]);overlaps=np.array([q['overlap'] for q in sites]);energies=np.array([q['U'] for q in sites]);reactions=np.array([q['reaction'] for q in sites]);centers=np.array([q['center'] for q in sites]);velocities=np.array([q.velocity(t) for q in self.sites])
        return dict(force=forces.sum(0),U=float(energies.sum()),reaction=reactions.sum(0),overlap=overlaps.max(0),contact_atom_ids=np.flatnonzero(overlaps.max(0)>0)+1,site_forces=forces,site_overlaps=overlaps,site_energies=energies,site_reactions=reactions,site_centers=centers,site_velocities=velocities,site_contact_counts=np.sum(overlaps>0,axis=1))
    def works(self,x0,x1,t0,t1,order=4):return np.array([q.work(x0,x1,t0,t1,order) for q in self.sites])
    def work(self,x0,x1,t0,t1,order=4):return float(self.works(x0,x1,t0,t1,order).sum())
    def phase(self,t):return 'POST_WITHDRAWAL' if t>=self.clear_time else 'MULTISITE_ACTIVE_PROTOCOL'

def freeze():
    import datetime
    assert not (O/'LOAD_SITE_SELECTION.json').exists() and not (O/'MULTISITE_LOADING_PROTOCOL.json').exists()
    audit=json.loads((O/'CRACK_BRIDGE_AUDIT.json').read_text());ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);X=ref['X'][:n]
    old=json.loads((HISTORY/'INDENTER_GEOMETRY.json').read_text());radius=old['radius_m'];penalty=old['penalty_K_N_m2']
    # One horizon plus half a node spacing above/below the original notch.
    zvalues=np.unique(X[:,2]);z=float(zvalues[np.argmin(abs(zvalues-(HORIZON+DX/2)))])
    points=[np.array([X[:,0].min(),-1.5*DX,z]),np.array([X[:,0].max(),1.5*DX,-z])];normals=[np.array([1.,0.,0.]),np.array([-1.,0.,0.])]
    sites=[]
    for i,(point,direction) in enumerate(zip(points,normals)):
        assert np.linalg.norm(X-point,axis=1).min()<1e-14
        center=point-direction*(radius+GAP);initial_gap=np.linalg.norm(X-center,axis=1).min()-radius;assert initial_gap>=GAP-1e-14
        sites.append(dict(fix_id='indenter'+str(i),surface_point_m=point,initial_center_m=center,inward_normal=direction,radius_m=radius,penalty_K_N_m2=penalty,initial_gap_m=initial_gap,mechanism='Offset radial normal contact, frictionless; not prescribed shear traction',rationale='Upper left notch lip and lower right surviving ligament form an opening/bending couple across the original notch; Y offsets break exact mirror symmetry without a friction law'))
    separation=float(np.linalg.norm(points[1]-points[0]));cs=np.sqrt(G/RHO);wave_time=separation/cs;delay=round(.5*wave_time/(50e-6))*50e-6;assert abs(delay-.0004)<1e-14
    min_x_distance=points[1][0]-points[0][0]+2*radius-2*DEPTH;min_distance=float(np.sqrt(min_x_distance**2+(points[1][1]-points[0][1])**2+(points[1][2]-points[0][2])**2));assert min_distance>2*radius
    final=audit['cases']['uniform']['final']['regions'][0]
    dump(O/'LOAD_SITE_SELECTION.json',dict(status='FROZEN',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sites=sites,surface_site_separation_m=separation,separation_over_horizon=separation/HORIZON,minimum_sphere_center_separation_bound_m=min_distance,actor_overlap_impossible_for_frozen_paths=True,selection_evidence=dict(audit_sha256=sha(O/'CRACK_BRIDGE_AUDIT.json'),upper_flap_surviving_bridges=final['surviving_crossing_bonds'],long_bridges=final['long_bonds_ge_2dx'],right_of_notch_tip=final['right_of_notch_tip'],selection='External side faces bracketing the residual notch ligament; centers offset by one horizon plus half spacing in Z, 1.5 spacing in Y. No interior load source, no chosen bond failures.'),historical_radius_and_penalty_source=str(HISTORY/'INDENTER_GEOMETRY.json'),historical_geometry_sha256=sha(HISTORY/'INDENTER_GEOMETRY.json')))
    cases={name:dict(delays_s=offsets,clear_time_s=CLEAR+max(offsets)) for name,offsets in [('simultaneous',[0.,0.]),('staggered',[0.,delay])]}
    dump(O/'MULTISITE_LOADING_PROTOCOL.json',dict(status='FROZEN',label='MULTISITE_MECHANICAL_SURROGATE',cases=cases,approach_s=APPROACH,peak_s=PEAK,retract_start_s=RETRACT,return_s=CLEAR,initial_gap_m=GAP,normal_depth_m=DEPTH,delay_rationale=dict(shear_speed_m_s=cs,site_separation_shear_travel_s=wave_time,chosen_delay_s=delay,rule='Half estimated site-to-site shear travel time rounded to 50 us; frozen before any formal trajectory'),end_s=END,common_cohort_time_s=.0028,numerics=dict(base_dt_s=BASE,predictive_subdivision_floor_s=MIN_H,max_damage_increment=.1,relative_energy_gate=.01,numerical_damping=0),comparison='Same sites, material, trajectories and end time; only onset timing differs. Not equal work. Post-return observation and actual no-contact windows reported separately.'))
    print('LOAD_SITES_FROZEN',[(q['surface_point_m']*1e6).tolist() for q in sites],'delay_ms',delay*1e3)

if __name__=='__main__':freeze()
