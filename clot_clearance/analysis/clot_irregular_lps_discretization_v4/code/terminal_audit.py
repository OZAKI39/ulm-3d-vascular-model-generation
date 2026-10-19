"""Read the last legal native restart; no extra dynamics or parameter case."""
from common import *
from native import native
from current_surface import CurrentSurface
import ctypes
def main():
 case=O/'cases/H_IRREGULAR_LPS_V4';result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE'];last=json.loads((case/'LAST_LEGAL_CHECKPOINT.json').read_text());checkpoint=case/last['path'];saved=load_npz(checkpoint/'accepted.npz');a,_=native(O/'tests/terminal_restart',restart=checkpoint/'native.restart');s=a.dynamics();a.tool.restore_checkpoint(checkpoint/'surface.npz',s)
 errors={k:float(abs(s[k]-saved[k]).max()) for k in ['x','velocity','force','d','kappa']};assert errors['x']==0 and errors['velocity']==0 and errors['d']==0 and errors['kappa']==0 and errors['force']<1e-14
 count,raw=a.l.gather_bonds();actual=np.ctypeslib.as_array(raw,shape=(3*count,)).reshape(-1,3)[:,1:].copy();actual.sort(1);expected=a.ref['bonds'][s['d']<1]+1;N=a.N
 assert np.array_equal(np.sort(actual[:,0]*N+actual[:,1]),np.sort(expected[:,0]*N+expected[:,1]))
 material=load_npz(O/'BIPHASIC_FIELD.npz');assert np.array_equal(a.field(13),material['bond_s0']) and np.array_equal(a.field(14),material['bond_sc']);assert np.array_equal(a.field(9),a.ref['original_bond_ids']) and np.array_equal(a.field(15),a.ref['V']);assert np.array_equal(a.gather('rmass'),RHO*a.ref['V'])
 base=a.tool.base;oldvol=base.volume(a.ref['regular_X']);ratios=base.ref_volume/oldvol;assert ratios.min()>0 and np.isfinite(ratios).all()
 audit=dict(status='PASS',checkpoint=last,checkpoint_time_s=float(s['time']),matches_formal_accepted_end=abs(float(s['time'])-result['accepted_time_s'])<1e-14,restart_state_errors=errors,original_atom_bond_volume_threshold_identity_exact=True,native_live_bond_graph_exact=True,nonuniform_rmass_exact=True,L1_reference_signed_volume_ratio_min=float(ratios.min()),L1_reference_signed_volume_ratio_max=float(ratios.max()),no_additional_dynamic_steps=True)
 a.close();dump(O/'TERMINAL_RESTART_AUDIT.json',audit);print('TERMINAL_NATIVE_AUDIT_PASS',flush=True)
if __name__=='__main__':main()
