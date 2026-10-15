from native import *
import resource
assert json.loads((O/'RESTART_DAMAGE_HISTORY_AUDIT.json').read_text())['status']=='PASS'
start=time.perf_counter();ref=np.load(PREV/'provenance/reference.npz');X=ref['X'];b=ref['bonds'];mask=X[:,0]<-30e-6;cross=mask[b[:,0]]!=mask[b[:,1]]
# Same real two-component cut; partial states are diagnostic, full damage is explicitly committed.
d=.4*cross;k=np.where(cross,.24*.15/(.24-.4*(.24-.15)),0.);ns=Native('postbreak_contact',contact=True);ns.import_state(k,d)
closed=X.copy();closed[:,0]+=np.where(mask,1.,-1.)*2e-6
partial=ns.eval(closed);partialextra=partial['force']-partial['force_cpp']
opened=X.copy();opened[:,0]+=np.where(mask,-1.,1.)*80e-6
ns.eval(opened);failed=ns.commit(0);repeated=ns.commit(0);closed_state=ns.eval(closed);fc=closed_state['force']-closed_state['force_cpp']
# Exact current pair_bpm_peri pairwise penalty, count each physical pair once.
i,j=np.triu_indices(341,1);dy=closed[i]-closed[j];rr=np.linalg.norm(dy,axis=1);active=rr<DX;bondset={tuple(sorted(map(int,bb))) for bb in b[~cross]};eligible=np.array([tuple((int(a),int(c))) not in bondset for a,c in zip(i,j)]);select=active&eligible
cs=18*K/(np.pi*DELTA**4);fm=-15*cs*(DX**3)**2*(rr[select]-DX)/DELTA
v=fm[:,None]*dy[select]/rr[select,None];expected=np.zeros_like(X);np.add.at(expected,i[select],v);np.add.at(expected,j[select],-v)
err=float(abs(fc-expected).max());scale=float(abs(expected).max());passed=float(abs(partialextra).max())<2e-19 and abs(failed['tdir'][:12047][cross]).max()==0 and abs(failed['tdir'][12047:][cross]).max()==0 and np.all(failed['d'][cross]==1) and np.all(failed['d'][~cross]==0) and failed['scalars'][3]==cross.sum() and np.array_equal(failed['scalars'][:6],repeated['scalars'][:6]) and scale>0 and err/scale<2e-10 and np.array_equal(closed_state['d'],failed['d'])
result=dict(status='PASS' if passed else 'FAIL',contact_distance_factor=1.0,source_pair_style='Unchanged locked pair_bpm_peri.cpp; symmetric volume-averaged penalty',partial_extra_pair_force_max_N=float(abs(partialextra).max()),committed_cross_failed_bonds=int(cross.sum()),actual_failure_events=int(failed['scalars'][3]),repeat_commit_no_extra_events=bool(np.array_equal(failed['scalars'][:6],repeated['scalars'][:6])),open_noncontact_cross_force_max_N=float(abs((failed['tdir'][:12047]+failed['tdir'][12047:])[cross]).max()),contact_pairs_once=int(select.sum()),contact_force_scale_N=scale,contact_force_error_max_N=err,contact_force_relative=err/scale,net_contact_force_N=float(np.linalg.norm(fc.sum(0))),closed_damage_unchanged=bool(np.array_equal(closed_state['d'],failed['d'])),native_active_bonds=closed_state['bonds'],scope='341-node diagnostic split/closure only; prescribed large opening used solely to exercise full failure, not strain validity or material calibration',wall_s=time.perf_counter()-start,RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
dump(O/'POSTBREAK_CONTACT_AUDIT.json',result);np.savez_compressed(O/'postbreak_contact.npz',reference_X=X,bonds=b,cross=cross,partial_x=partial['x'],partial_d=partial['d'],partial_pair_force_N=partialextra,opened_x=failed['x'],failed_d=failed['d'],closed_x=closed_state['x'],closed_d=closed_state['d'],contact_force_N=fc,expected_contact_force_N=expected,bond_ids=failed['bond_ids']);ns.close();print('CONTACT',result,flush=True);assert passed
