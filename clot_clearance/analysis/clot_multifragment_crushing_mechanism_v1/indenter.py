"""Independent physical-SI oracle for LAMMPS fix indent (not a force callback)."""
from common import *
APPROACH=.0002; PEAK=.0012; RETRACT=.0014; CLEAR=.0024; END=.0036
GAP=25e-6; DEPTH=240e-6
def ramp(t,a,b):
 u=np.clip((t-a)/(b-a),0.,1.);return .5*(1-np.cos(np.pi*u))
def ramp_v(t,a,b):
 return np.pi/(2*(b-a))*np.sin(np.pi*(t-a)/(b-a)) if a<t<b else 0.
def displacement(t):return GAP*ramp(t,0,APPROACH)+DEPTH*ramp(t,APPROACH,PEAK)-(GAP+DEPTH)*ramp(t,RETRACT,CLEAR)
def speed(t):return GAP*ramp_v(t,0,APPROACH)+DEPTH*ramp_v(t,APPROACH,PEAK)-(GAP+DEPTH)*ramp_v(t,RETRACT,CLEAR)
def phase(t):return 'APPROACH' if t<APPROACH else 'INDENT' if t<PEAK else 'HOLD' if t<RETRACT else 'RETRACT' if t<CLEAR else 'FREE_EVOLUTION'
class Indenter:
 def __init__(self,ref):
  self.n=int(ref['clot_count']);self.N=len(ref['X']);self.config=json.loads((O/'INDENTER_GEOMETRY.json').read_text());self.c0=np.array(self.config['initial_center_m']);self.radius=self.config['radius_m'];self.k=self.config['penalty_K_N_m2']
 def center(self,t):return self.c0+np.array([0.,0.,-displacement(t)])
 def evaluate(self,x,t):return self.at(x,self.center(t))
 def at(self,x,c):
  dv=x[:self.n]-c;r=np.linalg.norm(dv,axis=1);overlap=np.maximum(self.radius-r,0.);assert np.all(r>0)
  force=np.zeros_like(x);force[:self.n]=self.k*overlap[:,None]**2*dv/r[:,None]
  return dict(force=force,U=float(self.k/3*np.sum(overlap**3)),reaction=-force.sum(0),overlap=overlap,contact_atom_ids=np.flatnonzero(overlap>0)+1,center=np.array(c))
 def command(self):
  def expr(a,b):return f'((time>={b:.17g})+(time>{a:.17g})*(time<{b:.17g})*0.5*(1-cos(PI*(time-{a:.17g})/{b-a:.17g})))'
  z=f'{self.c0[2]:.17g}-{GAP:.17g}*{expr(0,APPROACH)}-{DEPTH:.17g}*{expr(APPROACH,PEAK)}+{GAP+DEPTH:.17g}*{expr(RETRACT,CLEAR)}'
  return f'variable indent_z equal "{z}"\nfix indenter clot indent {self.k:.17g} sphere {self.c0[0]:.17g} {self.c0[1]:.17g} v_indent_z {self.radius:.17g} units box\nfix_modify indenter energy no'
 def work(self,x0,x1,t0,t1,order=4):
  # Integrate partial U/partial c dot dc on the actual discrete, piecewise-linear
  # accepted path. Split exactly at contact entry/exit, then Gauss integrate.
  c0=self.center(t0);dc=self.center(t1)-c0
  if not dc.any():return 0.
  a=x0[:self.n]-c0;b=(x1[:self.n]-x0[:self.n])-dc
  aa=np.einsum('ij,ij->i',b,b);bb=2*np.einsum('ij,ij->i',a,b);cc=np.einsum('ij,ij->i',a,a)-self.radius**2
  disc=bb*bb-4*aa*cc;lo=np.zeros(self.n);hi=np.ones(self.n);active=(aa>0)&(disc>0)
  lo[active]=np.maximum(0,(-bb[active]-np.sqrt(disc[active]))/(2*aa[active]));hi[active]=np.minimum(1,(-bb[active]+np.sqrt(disc[active]))/(2*aa[active]))
  active=(active&(hi>lo))|((aa==0)&(cc<0));idx=np.flatnonzero(active)
  if not len(idx):return 0.
  points,weights=np.polynomial.legendre.leggauss(order);half=(hi[idx]-lo[idx])/2;mid=(hi[idx]+lo[idx])/2;work=0.
  for p,w in zip(points,weights):
   s=mid+half*p;dv=a[idx]+s[:,None]*b[idx];rr=np.linalg.norm(dv,axis=1);delta=np.maximum(self.radius-rr,0)
   work+=float(np.sum(w*half*self.k*delta**2*(dv@dc)/rr))
  return work
