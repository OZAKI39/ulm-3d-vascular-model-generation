from mechanics import *
import resource
assert not (O/'EXTERNAL_LEDGER_CHECK.json').exists()
g=json.loads((O/'PREFLIGHT_GATES.json').read_text());ns=Native('preflight/ledger');s=ns.dynamics();m=Mechanics(ns,s,g);rows=[];rows.append(m.check(s,s,0)[0]);p=s;times=[]
for k in range(8):
 ts=time.perf_counter();ns.l.command('run 1 pre no post no');s=ns.state();row,_,_,_=m.check(p,s,BASE);rows.append(row);ns.ext.accept(s);p=s;times.append(time.perf_counter()-ts)
assert all(row['token']==row['step'] for row in rows[1:])
assert not any(row['bulk_softened'] or row['interface_softened'] for row in rows)
save_state(O/'preflight/ledger_terminal.npz',s)
dump(O/'EXTERNAL_LEDGER_CHECK.json',dict(status='PASS',physical_steps=8,time_s=s['time'],rows=rows,max_independent_errors=m.maxerr,median_complete_checked_step_wall_s=float(np.median(times)),estimated_9320_steps_wall_s=float(np.median(times)*END/BASE),estimate_excludes_extra_damage_evaluation_and_event_restarts=True,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,load_uses_current_native_physical_time=True,anchors_not_integrated=True,wall_and_drag_work_separate=True))
ns.close();print('EXTERNAL_LEDGER_PASS',np.median(times),max(z['relative_R'] for z in rows),flush=True)
