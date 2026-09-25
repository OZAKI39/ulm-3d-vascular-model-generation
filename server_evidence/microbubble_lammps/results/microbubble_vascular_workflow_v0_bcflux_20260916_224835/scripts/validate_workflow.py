#!/usr/bin/env python3
"""Independent readback of actual CSVs; never promotes engine status to PASS."""
from pathlib import Path
import csv,json,sys,itertools,collections,hashlib,math
import numpy as np
from flow_geometry import Geometry,FrozenSampler
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S/'provenance/stable_rigid_source'));sys.path.insert(0,str(S/'inputs'))
from reference_rigid_sphere_lubrication import FlowGradientReference,dense_solve,coefficients
from sonovue_sampler import SonoVueDistribution
g=Geometry(S/'geometry/GEOMETRY_ARRAYS.npz');sampler=FrozenSampler(S/'fields/FROZEN_FLOW_FIELD_V0.h5');gradient=FlowGradientReference(S/'fields/FROZEN_FLOW_FIELD_V0.h5')
manifest=json.loads((S/'geometry/BOUNDARY_MANIFEST.json').read_text());section=dict(np.load(S/'geometry/INJECTION_SECTION.npz'));sizecap=json.loads((S/'validation/INLET_SIZE_CAPACITY.json').read_text())
Q=json.loads((S/'validation/AUTHORITATIVE_INLET_FLOW_AUDIT.json').read_text())['authoritative_flow_rate'];stream=np.loadtxt(S/'inputs/size_stream.txt');dist=SonoVueDistribution(S/'inputs/FROZEN_SONOVUE_HISTOGRAM.csv')
assert np.array_equal(stream[:,1],dist.inverse_cdf(stream[:,0])*.5e-6)
def rows(path):
    with path.open() as f:r=list(csv.DictReader(f))
    assert all(x['disclaimer']=='NOT EXPERIMENTAL CONCENTRATION' for x in r),str(path)
    return r
def num(r,k):return float(r[k])
def ident(r):return int(r['particle_id'])
def point(r,prefix=''):return np.array([float(r[prefix+k]) for k in ('x_m','y_m','z_m')])
def vel(r):return np.array([float(r[k]) for k in ('vx_m_s','vy_m_s','vz_m_s','omega_x','omega_y','omega_z')])
def safe_segment(a,b,radius,margin=1e-10,depth=0):
    da=g.distance(a)[0];db=g.distance(b)[0];length=np.linalg.norm(b-a)
    assert min(da,db)>=radius+margin-2e-14,('WALL_PENETRATION',min(da,db)-radius)
    if min(da,db)-length*.5>=radius+margin-2e-14:return
    assert depth<25,'INDEPENDENT_WALL_CERTIFICATION_DEPTH'
    m=(a+b)*.5;safe_segment(a,m,radius,margin,depth+1);safe_segment(m,b,radius,margin,depth+1)
def readcase(D):
    series=rows(D/'flux_timeseries.csv');ev=rows(D/'injection_events.csv');tr=rows(D/'TRAJECTORIES.csv');st=rows(D/'INTEGRATION_STAGES.csv');port=rows(D/'PORT_EVENTS.csv');near=rows(D/'NEAR_WALL_EVENTS.csv');pairs=rows(D/'PAIR_HYDRODYNAMIC_EVENTS.csv');hist=rows(D/'SOLVER_HISTORY.csv');retry=rows(D/'RETRY_HISTORY.csv')
    state=json.loads((D/'RUN_STATE.json').read_text());receipt=json.loads((D/'EXECUTION_RECEIPT.json').read_text());assert receipt['returncode']==0
    cfg={}
    if (D/'case.cfg').exists():
        assert hashlib.sha256((D/'case.cfg').read_bytes()).hexdigest()==receipt['config_sha256']
        for line in (D/'case.cfg').read_text().splitlines():
            if line and not line.startswith('#'):k,v=line.split();cfg[k]=v
    evstep=collections.defaultdict(list)
    for e in ev:evstep[int(e['step'])+1].append(e)
    exitsstep=collections.defaultdict(list)
    for e in port:exitsstep[int(e['step'])].append(e)
    source={};admitted={};rejected=set();active=set();pending=[];exits=[0,0,0,0];cumulative=0;max_active=0;max_pending=0;last_id=0
    for r in series:
        step=int(r['step']);dt=num(r,'dt_s');cumulative+=num(r,'target_number_flux')*dt
        assert abs(num(r,'authoritative_Q_in')/Q-1)<1e-14
        assert abs(num(r,'target_cumulative')-cumulative)<2e-10
        for e in evstep[step]:
            id=ident(e);radius=num(e,'radius_m');u=num(e,'source_u')
            if e['event']=='SOURCE_DRAWN':
                assert id==last_id+1;last_id=id;assert np.array_equal(stream[id-1],[u,radius]);source[id]=(u,radius);pending.append(id)
            else:
                assert source[id]==(u,radius),'CANDIDATE_IDENTITY_OR_SIZE_CHANGED'
                if e['event']=='SIZE_INADMISSIBLE_AT_INLET':
                    assert radius>sizecap['admissible_radius_upper_m'];rejected.add(id);pending.remove(id)
                elif e['event']=='ADMITTED':
                    assert pending[0]==id,'FIFO_REORDERED';pending.pop(0);admitted[id]=e;active.add(id)
                elif e['event']=='TRANSIENTLY_BLOCKED':assert pending[0]==id
                else:raise AssertionError(e)
        for e in exitsstep[step]:
            id=ident(e);assert id in active;active.remove(id);idx=[p['name'] for p in manifest['ports']].index(e['port']);exits[idx]+=1
        assert len(source)==len(rejected)+len(pending)+len(admitted)
        assert len(admitted)==len(active)+sum(exits)+int(r['N_terminal_failure'])
        expected={'N_source_drawn':len(source),'size_rejected_cumulative':len(rejected),'pending_count':len(pending),'admitted_cumulative':len(admitted),'active_count':len(active),'inlet_backflow_cumulative':exits[0],**{f'outlet_{i}_cumulative':exits[i+1] for i in range(3)}}
        assert all(int(r[k])==v for k,v in expected.items()),(step,expected,r)
        assert abs(num(r,'flux_deficit')-(cumulative-len(admitted)))<2e-10
        covered=len(source) if r['control_basis']=='SOURCE_POPULATION' else len(admitted)+len(pending)
        assert covered==math.floor(cumulative+1e-10),(step,covered,cumulative)
        max_active=max(max_active,len(active));max_pending=max(max_pending,len(pending))
    # Actual insertion/removal is checked against native LAMMPS counts every step.
    assert all(int(h['lammps_natoms'])==int(series[int(h['step'])]['active_count']) for h in hist)
    assert state['committed_LAMMPS_insertions']==len(admitted) and state['committed_LAMMPS_removals']==sum(exits)
    assert state['active']==len(active) and state['pending']==len(pending)
    pending_saved=rows(D/'PENDING_FINAL.csv');assert [ident(x) for x in pending_saved]==pending
    if state['reason']=='INJECTION_CAPACITY_EXCEEDED':
        assert len(pending)>int(cfg['max_pending']) or (len(pending)>0 and state['time_s']-num(pending_saved[0],'born_s')>float(cfg['pending_max_age']))
    for r in retry:assert int(r['inserted_count'])==int(r['rolled_back_count'])
    # All injection positions use native-equivalent VALID, inward flow, exact
    # finite-radius side-wall clearance and actual triangulated section.
    for id,e in admitted.items():
        x=point(e);a=num(e,'radius_m');status,u=sampler.query([x]);assert status[0]==0 and -u[0]@section['normal']>0
        assert g.distance(x)[0]>=a+sizecap['clearance_m']-2e-14
        assert abs((x-section['center'])@section['normal'])<1e-12
        found=False
        for t in section['triangles']:
            mat=np.column_stack([t[1]-t[0],t[2]-t[0]]);b=np.linalg.lstsq(mat,x-t[0],rcond=None)[0]
            if (b>=-1e-8).all() and b.sum()<=1+1e-8 and np.linalg.norm(mat@b-(x-t[0]))<1e-12:found=True;break
        assert found,'INJECTION_OUTSIDE_SECTION'
    # Geometric trajectory distances/region tags independent of C++ BVHs.
    trmap={};minimum_wall=1e100;region_ties=0;region_locators={}
    wall_classes=g.g['classes'][g.g['classes']<=2]
    for r in tr:
        x=point(r);a=num(r,'radius_m');d,_,idx=g.distance(x);minimum_wall=min(minimum_wall,d-a)
        assert abs(d-num(r,'wall_distance_m'))<2e-14
        if int(r['wall_region'])!=wall_classes[idx]:
            # At shared region edges the closest point can belong to two exact
            # triangle subsets. Verify the reported subset independently instead
            # of assuming VTK and C++ choose the same triangle in a distance tie.
            import vtk
            from flow_geometry import poly
            region=int(r['wall_region'])
            if region not in region_locators:
                loc=vtk.vtkStaticCellLocator();loc.SetDataSet(poly(g.g['points'],g.g['faces'][g.g['classes']==region]));loc.BuildLocator();region_locators[region]=loc
            closest=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.)
            region_locators[region].FindClosestPoint(x,closest,cell,sub,d2)
            assert abs(float(d2)**.5-d)<2e-14,'INCORRECT_WALL_REGION'
            region_ties+=1
        assert d-a>=1e-10-2e-14
        if r['state']=='ACTIVE':assert sampler.query([x])[0][0]==0
        key=(int(r['step']),ident(r));
        if key in trmap:assert np.array_equal(point(trmap[key]),x)
        trmap[key]=r
    bystage=collections.defaultdict(list)
    for r in st:bystage[(int(r['step']),int(r['stage']))].append(r)
    minimum_pair=1e100;max_mid_error=0.;max_end_error=0.;max_ref_error=0.;solves_checked=0;hydro_stages=set((int(r['step']),int(r['stage'])) for r in pairs)
    for (step,stage),records in sorted(bystage.items()):
        records.sort(key=ident);x=np.array([point(r) for r in records]);a=np.array([num(r,'radius_m') for r in records]);q=np.array([vel(r) for r in records]);dt=num(records[0],'dt_s');ids=[ident(r) for r in records]
        assert (sampler.query(x)[0]==0).all()
        starts=bystage[(step,0)];starts.sort(key=ident);base=np.array([point(r) for r in starts]);r0=np.array([vel(r) for r in starts])
        if stage==1:
            err=float(np.max(np.abs(x-(base+.5*dt*r0[:,:3]))));max_mid_error=max(max_mid_error,err);assert err<3e-20
            end=np.array([point(trmap[(step+1,id)]) for id in ids]);err=float(np.max(np.abs(end-(base+dt*q[:,:3]))));max_end_error=max(max_end_error,err);assert err<3e-20
            for i in range(len(ids)):safe_segment(base[i],x[i],a[i]);safe_segment(base[i],end[i],a[i])
            for i,j in itertools.combinations(range(len(ids)),2):
                r=base[i]-base[j];v=(end[i]-base[i])-(end[j]-base[j]);f=np.clip(-r@v/(v@v),0,1) if v@v else 0.;gap=np.linalg.norm(r+f*v)-a[i]-a[j];minimum_pair=min(minimum_pair,gap);assert gap>=-1e-12
        # Independent dense reference checks all hydrodynamic stages and a
        # deterministic sparse sample of uncoupled steps (same source model).
        if (step,stage) in hydro_stages or step%137==0 or step<5:
            U,O,E,_=gradient.query(x);allpairs=list(itertools.combinations(range(len(ids)),2));ref=dense_solve(x,a,U,O,E,allpairs,base=base,dt=dt*(.5 if stage==0 else 1))
            err=float(np.max(np.abs(q-ref)*np.c_[np.ones((len(a),3)),np.repeat(a[:,None],3,axis=1)]));max_ref_error=max(max_ref_error,err);assert err<2e-10,(D.name,step,stage,err);solves_checked+=1
    # Check finite-gap coefficients from the independently coded reference.
    for e in pairs:
        rr=bystage[(int(e['step']),int(e['stage']))];rr={ident(r):r for r in rr};a=num(rr[int(e['id_i'])],'radius_m');b=num(rr[int(e['id_j'])],'radius_m');c=coefficients(a,b,num(e,'gap_m'))
        for key,i in [('normal_resistance',0),('shear_resistance',1),('pump_resistance',2)]:assert abs(num(e,key)/c[i]-1)<2e-11
    # Ports: segment-plane direction and cap membership, independently using VTK
    # cell locator at the recorded intersection (caps never queried as walls).
    import vtk
    from flow_geometry import poly
    for e in port:
        portinfo=next(p for p in manifest['ports'] if p['name']==e['port']);n=np.array(portinfo['outward_unit_normal']);start=np.array([num(e,'start_'+k) for k in 'xyz']);end=np.array([num(e,'end_'+k) for k in 'xyz']);cross=np.array([num(e,'cross_'+k) for k in 'xyz']);f=num(e,'fraction')
        assert 0<f<=1 and (end-start)@n>0;assert np.max(np.abs(cross-(start+f*(end-start))))<3e-20
        reader=vtk.vtkSTLReader();reader.SetFileName(str(S/portinfo['stl']));reader.Update();loc=vtk.vtkStaticCellLocator();loc.SetDataSet(reader.GetOutput());loc.BuildLocator();q=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);d=vtk.reference(0.);loc.FindClosestPoint(cross,q,cell,sub,d);assert float(d)<1e-24
    # Biological near-wall transitions and right-censored residence use the same
    # explicitly discrete observation times, recomputed from distances above.
    active_near={};expected_near=[]
    for r in tr:
        id=ident(r);t=num(r,'time_s');d=g.distance(point(r))[0];a=num(r,'radius_m');isnear=int(r['wall_region'])==0 and (d-a)/a<=.2
        if isnear and id not in active_near:active_near[id]=t;expected_near.append((id,'ENTER',t,0.))
        if (not isnear or r['state']=='EXITED') and id in active_near:expected_near.append((id,'LEAVE',t,t-active_near.pop(id)))
    for id,t in sorted(active_near.items()):expected_near.append((id,'RIGHT_CENSORED',state['time_s'],state['time_s']-t))
    actual_near=[(ident(r),r['event'],num(r,'time_s'),num(r,'residence_s')) for r in near];assert len(expected_near)==len(actual_near)
    for a,b in zip(expected_near,actual_near):assert a[:2]==b[:2] and np.max(np.abs(np.array(a[2:])-b[2:]))<1e-12
    result=dict(status='PASS',case=D.name,engine_reason=state['reason'],accepted_steps=len(hist),source_drawn=len(source),size_rejected=len(rejected),admitted=len(admitted),active=len(active),pending=len(pending),max_active=max_active,max_pending=max_pending,outlet_exits=sum(exits[1:]),inlet_backflow_exits=exits[0],hydrodynamic_observations=len(pairs),independent_dense_solves=solves_checked,max_dense_reference_scaled_error_m_s=max_ref_error,max_RK2_midpoint_error_m=max_mid_error,max_RK2_endpoint_error_m=max_end_error,min_all_pair_swept_gap_m=minimum_pair,min_wall_gap_m=minimum_wall,biological_nearwall_events=len(near),rollback_attempts=len(retry),checks=['independent_flux_integral','SOURCE_ADMITTED_semantics','frozen_source_draws','FIFO','particle_accounting','native_LAMMPS_counts','VALID_support_injection','wall_clearance','full_pair_swept_safety','all_RK2_positions','sampled_independent_resistance_solve','port_geometry','nearwall_reconstruction'],limitations=['Geometric stall retains particles ACTIVE at last accepted state; terminal_failure=0 is not silent deletion','No outlet observed unless outlet_exits > 0'],disclaimer='NOT EXPERIMENTAL CONCENTRATION')
    (S/'validation'/f'{D.name}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True);return result
def main():
    names=sys.argv[1:];dirs=[S/'runs'/n for n in names] if names else sorted(p.parent for p in (S/'runs').glob('*/RUN_STATE.json'))
    results=[readcase(d) for d in dirs]
    (S/'validation/INDEPENDENT_WORKFLOW_VALIDATION.json').write_text(json.dumps(dict(status='PASS',cases=results),indent=2)+'\n')
if __name__=='__main__':main()
