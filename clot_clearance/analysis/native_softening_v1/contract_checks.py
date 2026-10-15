from native import *
import resource
start=time.perf_counter();records=[]
# Wrong native restart version must produce a hard error, never zero history.
p=O/'partial_damage.restart';raw=p.read_bytes();assert raw.count(b'LPSVAR_RESTART_V1')==1
bad=O/'incompatible_version.restart';bad.write_bytes(raw.replace(b'LPSVAR_RESTART_V1',b'LPSVAR_RESTART_V0'))
l=lammps(cmdargs=['-log',str(O/'rejected_restart.log'),'-screen','none'],comm=COMM)
try:l.command(f'read_restart {bad}');raise AssertionError('incompatible restart accepted')
except Exception as e:
 text=str(e);assert 'Missing or incompatible variational restart history' in text,text;records.append(dict(case='incompatible_restart',status='PASS_EXPECTED_REJECTION',message=text))
finally:l.close()
# Diagnostic constructor validates collisions rather than aliasing two original IDs.
lines=(O/'provenance/reference.native').read_text().splitlines();a=lines[343].split();bb=lines[344].split();bb[1:3]=a[1:3];lines[344]=' '.join(bb);dup=O/'duplicate_identity.native';dup.write_text('\n'.join(lines)+'\n')
l=lammps(cmdargs=['-log',str(O/'rejected_duplicate.log'),'-screen','none'],comm=COMM)
try:
 l.commands_string(f' units si\natom_style bond\nnewton on off\nread_data {O/"provenance/body.data"}\nbond_style bpm/peri/variational {dup}\nbond_coeff 1 {K:.17g} {G} {DELTA} .15 .24');raise AssertionError('duplicate endpoints accepted')
except Exception as e:
 text=str(e);assert 'Duplicate endpoint pairs' in text,text;records.append(dict(case='duplicate_identity',status='PASS_EXPECTED_REJECTION',message=text))
finally:l.close()
ns=Native('commit_contract');X=np.load(PREV/'provenance/reference.npz')['X'];ns.eval(X*1.165);before=ns.commit(0);ns.eval(X*1.165+np.array([1e-6,0,0]));rc=ns.lib.lpsvar_action(ns.l.lmp,2,0,None,None);after=ns.snapshot();assert rc==-5 and np.array_equal(before['d'],after['d']) and np.array_equal(before['kappa'],after['kappa']) and np.array_equal(before['scalars'][1:6],after['scalars'][1:6]);records.append(dict(case='reused_token_different_coordinates',status='PASS_EXPECTED_REJECTION',return_code=rc,history_and_ledger_unchanged=True));ns.close()
ns=Native('legacy_C_API_rejection',legacy=True);buf=np.empty(10);rc=ns.lib.lpsvar_output(ns.l.lmp,8,ns.ptr(buf));assert rc==-99;ns.close();records.append(dict(case='old_style_not_softening_history',status='PASS_EXPECTED_REJECTION',return_code=rc))
dump(O/'CONTRACT_REJECTION_CHECKS.json',dict(status='PASS',records=records,wall_s=time.perf_counter()-start,RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss));print('CONTRACT_CHECKS_PASS',flush=True)
