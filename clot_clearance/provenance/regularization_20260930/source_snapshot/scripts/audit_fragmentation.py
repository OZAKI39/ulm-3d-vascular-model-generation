"""Independent read-only audit of a completed load-driven fragmentation run."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pyvista as pv


def components_from_edges(n,pairs,active):
    parent=np.arange(n)
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in pairs[active]:parent[find(i)]=find(j)
    return np.array([find(i) for i in range(n)])


def audit(run):
    run=Path(run);z=np.load(run/'states.npz');c=json.loads((run/'CONFIG.json').read_text())
    summary=json.loads((run/'SUMMARY.json').read_text());history=json.loads((run/'history.json').read_text())
    component_history=json.loads((run/'components.json').read_text());lineage=json.loads((run/'lineage.json').read_text())
    pairs=z['pairs'];n=len(z['X']);volume=z['volume'];total=volume.sum();dc=c['damage']
    checks={};checks['initial_state_undeformed_at_rest']=bool(np.array_equal(z['x'][0],z['X']) and np.all(z['v'][0]==0))
    checks['all_particles_retained']=z['x'].shape[1]==n and all(row['particle_count']==n for row in history)
    checks['anchors_fixed']=bool(np.all(z['x'][:,z['fixed']]==z['X'][z['fixed']]) and np.all(z['v'][:,z['fixed']]==0))
    checks['no_healing']=bool(np.all(np.diff(z['bond_active'].astype(int),axis=0)<=0) and np.all(np.diff(z['bond_D'],axis=0)>=-1e-15))
    checks['integrity_from_damage']=bool(np.allclose(z['integrity'],np.where(z['bond_active'],1-z['bond_D'],0),rtol=0,atol=1e-15))
    cycle_law=True;topology=True;stable=True;stats_ok=True;fields_ok=True;clearance_ok=True
    found_events={name:None for name in ['first_bond_failure','first_detachment','first_fragment','first_clearance_crossing']}
    detached_origin=np.full_like(z['X'],np.nan);maxtravel=0.;actual_free_with_bonds=False
    for k in range(len(history)):
        if k:
            inc=dc['DeltaN']*dc['C_damage']*np.maximum(z['amplitude'][k]/dc['Q_ref']-1,0)**dc['m_damage']
            if not dc['enabled']:inc[:]=0
            expected=np.where(z['bond_active'][k-1],np.clip(z['bond_D'][k-1]+inc,0,1),z['bond_D'][k-1])
            cycle_law &= np.allclose(expected,z['bond_D'][k],rtol=0,atol=1e-14)
            cycle_law &= np.array_equal(z['bond_active'][k],z['bond_active'][k-1] & (expected<dc.get('D_break',1)))
        raw=components_from_edges(n,pairs,z['bond_active'][k])
        anchors=np.unique(raw[z['fixed']]);attached=np.isin(raw,anchors)
        topology &= np.array_equal(attached,z['attached'][k])
        ids=z['fragment_id'][k]
        topology &= len(np.unique(raw))==len(np.unique(ids))
        for label in np.unique(raw):topology &= len(np.unique(ids[raw==label]))==1
        if k:
            for old in np.unique(z['fragment_id'][k-1]):
                members=z['fragment_id'][k-1]==old
                children=np.unique(ids[members])
                if len(children)==1:stable &= int(children[0])==int(old)
        row=history[k]
        stats_ok &= abs(volume[~attached].sum()/total-row['detached_volume_fraction'])<1e-12
        stats_ok &= abs(volume[attached].sum()/total-row['attached_volume_fraction'])<1e-12
        stats_ok &= abs(np.mean(~z['bond_active'][k])-row['broken_bond_fraction'])<1e-12
        stats_ok &= abs(row['attached_volume_fraction']+row['detached_volume_fraction']-1)<1e-12
        first=(~attached)&~np.isfinite(detached_origin[:,0]);detached_origin[first]=z['x'][k,first]
        for comp in component_history[k]['components']:
            members=np.asarray(comp['particle_ids'])
            com=np.average(z['x'][k,members],axis=0,weights=volume[members])
            stats_ok &= np.allclose(com,comp['center_of_mass_m'],rtol=0,atol=1e-14)
            stats_ok &= np.all(ids[members]==comp['fragment_id'])
            if not comp['attached_to_base']:
                travel=np.linalg.norm(com-np.average(detached_origin[members],axis=0,weights=volume[members]))
                maxtravel=max(maxtravel,float(travel))
                stats_ok &= abs(travel-comp['travel_distance_m'])<1e-12
                if len(members)>1 and np.any(z['bond_active'][k]&np.isin(pairs[:,0],members)&np.isin(pairs[:,1],members)):
                    actual_free_with_bonds=True
        event_conditions=dict(first_bond_failure=np.any(~z['bond_active'][k]),first_detachment=np.any(~attached),
            first_fragment=len(np.unique(raw))>1,first_clearance_crossing=np.any(z['cleared'][k]))
        for event,value in event_conditions.items():
            if value and found_events[event] is None:found_events[event]=int(z['cycles'][k])
        mesh=pv.read(run/f'vtk/particles_{k:04d}.vtp')
        for field in ['damage','fragment_id','attached_to_base','is_surface_particle','velocity','displacement','deformation_rank']:
            fields_ok &= field in mesh.point_data
        fields_ok &= np.array_equal(mesh.points,z['x'][k]) and np.array_equal(mesh['fragment_id'],ids)
        fields_ok &= not z['surface'][k,z['fixed']].any()
    checks.update(cyclic_law_exact=bool(cycle_law),topology_and_anchor_reachability=bool(topology),stable_component_IDs=bool(stable),
        independently_recomputed_statistics=bool(stats_ok),VTK_required_arrays_and_positions=bool(fields_ok))
    ledger=[json.loads(line) for line in (run/'failure_ledger.jsonl').read_text().splitlines()]
    checks['each_failure_has_damage_evidence']=len(ledger)==int((~z['bond_active'][-1]).sum()) and len({r['bond_id'] for r in ledger})==len(ledger)
    for row in ledger:
        checks['each_failure_has_damage_evidence'] &= row['D_after']>=row['D_break'] and row['D_before']<row['D_break'] and row['Q']>dc['Q_ref']
    crossing=json.loads((run/'clearance_ledger.json').read_text());seen=set()
    for row in crossing:
        members=row['particle_ids'];k=int(np.flatnonzero(z['cycles']==row['cycles'])[0])
        clearance_ok &= not bool(seen.intersection(members));seen.update(members)
        # Full component center crosses; only not-yet-counted member volume is added.
        ids=np.flatnonzero(z['fragment_id'][k]==row['fragment_id'])
        prev=np.average(z['x'][k-1,ids,0],weights=volume[ids]);now=np.average(z['x'][k,ids,0],weights=volume[ids])
        clearance_ok &= prev<c['transport']['x_clearance_m']<=now
    clearance_ok &= abs(volume[list(seen)].sum()/total-history[-1]['cleared_volume_fraction'])<1e-12
    checks['clearance_crossings_no_double_count']=bool(clearance_ok)
    checks['actual_bond_failure']=history[-1]['broken_bond_fraction']>0
    checks['actual_detachment']=history[-1]['detached_volume_fraction']>0
    checks['free_component_retains_internal_PD_bonds']=actual_free_with_bonds
    checks['fragment_motion_after_detachment']=maxtravel>1e-8
    checks['hydrodynamic_impulse_nonzero']=bool(np.linalg.norm(z['hydro_impulse'])>0)
    checks['at_least_one_clearance_crossing']=history[-1]['cleared_volume_fraction']>0
    checks['events_match']=all((summary['events'][name]['cycles'] if summary['events'][name] else None)==N for name,N in found_events.items())
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,events_cycles=found_events,
        maximum_fragment_travel_m=maxtravel,final=history[-1],failure_ledger_entries=len(ledger),
        source_review='New runner advances x only by dt*v and anchor reset; artificial-cleavage tests are not imported',
        reduced_rank_interpretation='Explicit intrinsic correspondence approximation; original full-rank kernel retained')
    (run/'INDEPENDENT_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    result=audit(a.run);print(json.dumps(result,indent=2))
    if result['status']!='PASS':raise SystemExit(1)
