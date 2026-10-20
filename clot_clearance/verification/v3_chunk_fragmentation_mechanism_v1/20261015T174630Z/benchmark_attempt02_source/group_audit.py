"""Energy and force diagnostics on the exact historical 13-node identity set."""
from common import *
from graph_ops import topology,boundary
def star_energy(m,x,d):
 e,_=m.geometry(x);ee=e[m.edge];dd=d[m.edge];den=m.Ddir*dd+m.krefdir*(1-dd);w=m.krefdir*(1-dd)/den;h=dd/den
 num=np.bincount(m.center,weights=m.l*w*ee,minlength=m.n);a=1/m.A+np.bincount(m.center,weights=m.l**2*h,minlength=m.n)
 return .5*np.bincount(m.center,weights=m.Ddir*w*ee**2,minlength=m.n)+.5*num*num/a
def main():
 ref=reference();model=energy_model(ref);co=json.loads((O/'COHESIVE_ADHESIVE_COMPETITION.json').read_text());nodes=np.asarray(co['historical_maximum_fragment']['atom_ids'])-1;out,inside=boundary(ref,nodes);states=[s for _,s in keyframes()]+[load_npz(p) for p in (O/'derived/recovered_accepted').glob('*.npz')];states=sorted(states,key=lambda s:int(s['step']));rows=[]
 for s in states:
  v=s['velocity'][nodes];vcom=v.mean(0);vh=s['velocity_force_half'][nodes];vhrel=vh-vh.mean(0);U=star_energy(model,s['x'],s['d']);assert abs(U.sum()-model.evaluate(s['x'],s['d'])['U'])<1e-18
  lab,sz,att=topology(ref,s['d']);ul,cnt=np.unique(lab[nodes],return_counts=True)
  rows.append(dict(step=int(s['step']),time_s=float(s['time']),star_U_J=float(U[nodes].sum()),K_translation_J=float(.5*len(nodes)*RHO*V*np.dot(vcom,vcom)),K_internal_J=float(.5*RHO*V*np.sum((v-vcom)**2)),load_relative_power_W=float(np.sum(s['load'][nodes]*vhrel)),drag_relative_power_W=float(np.sum(s['drag'][nodes]*vhrel)),numeric_relative_power_W=float(np.sum(s['numeric'][nodes]*vhrel)),load_sum_N=s['load'][nodes].sum(0),drag_sum_N=s['drag'][nodes].sum(0),component_intersections=cnt,attached_nodes=int(att[lab[nodes]].sum()),internal_live=int(np.sum(s['d'][inside]<1)),external_live=int(np.sum(s['d'][out]<1))))
 dump(O/'LOCAL_FRAGMENT_COHESION.json',dict(label='DIAGNOSTIC ONLY on ACTUAL DETACHED identities',atom_ids=nodes+1,states=rows,interpretation='For a genuinely isolated single component the sum of its endpoint-star energies is its independent stored energy. Before separation or after splitting this is only the energy of the chosen identity set. K_internal is kinetic energy relative to this set COM. Power is instantaneous using the archived force-evaluation half velocity; these samples are not a reconstructed work integral. Uniform component-COM drag does not directly damp internal relative motion of an intact free component. This is a testable secondary-breakup mechanism, not a demonstrated causal replacement for the initial failure mechanism.'))
 print('GROUP_AUDIT_DONE')
if __name__=='__main__':main()
