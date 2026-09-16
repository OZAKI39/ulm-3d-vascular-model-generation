"""Reference-only adapters to unmodified GPL RMBW. No upstream formulas copied.
The Lubrication API returns excess resistance; reconstruct total resistance by
adding the known isolated-sphere bulk diagonal, then solve six impulse systems.
Native units use fixed L*=1um and eta*=1mPa s to avoid upstream absolute COO
threshold 1e-12 discarding SI rotational entries. No cutoff or sign edits.
"""
from pathlib import Path
import os,sys,json,contextlib
import numpy as np
R=Path(__file__).resolve().parents[1]
W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['work'])
ROOT=W/'upstream/RigidMultiblobsWall'

def _sphere():
    sys.path[:0]=[str(ROOT/'sphere'),str(ROOT)]
    import sphere
    return sphere

@contextlib.contextmanager
def _runtime():
    p=W/'runtime/rmbw_sphere';p.mkdir(exist_ok=True,parents=True)
    table=p/'mobility.162-blob.dat'
    if not table.exists():table.symlink_to(ROOT/'sphere/mobility.162-blob.dat')
    old=os.getcwd();os.chdir(p)
    try:yield
    finally:os.chdir(old)

def _get_sphere(a,h,eta):
    with _runtime():
        sphere=_sphere();native=sphere.sphere_best_mobility_known(np.array([0.,0.,a+h]),eta,a)
    M=np.column_stack([native@np.eye(6)[:,k] for k in range(6)])
    return M,{'implementation':'sphere.sphere_best_mobility_known','native_units':'direct SI','native_basis_impulses':6,'table_center_height_over_a_range':[1.,10.],'table_evaluation':'extrapolation' if (a+h)/a>10 else 'interpolation','coupling_sign_correction_applied':False,'matrix_symmetrized':False,'all_blocks_available':True,'method':'Huang normal, Goldman/Faucheux parallel, archived 162-blob spline RR/TR'}

def _get_lub(a,h,eta):
    sys.path.insert(0,str(W/'build/rmbw'))
    import Lubrication_Class
    L=1e-6;eta0=1e-3;an=a/L;etan=eta/eta0;hn=(a+h)/L
    lc=Lubrication_Class.Lubrication(1e-5)
    data=[];rows=[];cols=[]
    lc.ResistCOO_wall([np.array([0.,0.,hn])],an,etan,0.,np.zeros(3),True,data,rows,cols)
    excess=np.zeros((6,6));excess[rows,cols]=data
    bulk=np.diag([6*np.pi*etan*an]*3+[8*np.pi*etan*an**3]*3)
    resistance=excess+bulk
    # Work-conjugate scaling q'=(V,a Omega), g'=(F,T/a).
    D=np.diag([1.]*3+[an]*3);invD=np.diag(1/np.diag(D))
    Rs=invD@resistance@invD
    cond=float(np.linalg.cond(Rs));assert cond<1e12,'Refuse unreliable resistance solve'
    Mnat=np.column_stack([invD@np.linalg.solve(Rs,invD@np.eye(6)[:,k]) for k in range(6)])
    # T,T: 1/(eta*L); mixed:1/(eta*L^2); R,R:1/(eta*L^3).
    row=np.array([1.]*3+[1/L]*3);M=Mnat*np.outer(row,row)/(eta0*L)
    return M,{'implementation':'Lubrication_Class.Lubrication.ResistCOO_wall(Sup_if_true=True) + isolated bulk diagonal','native_length_unit_m':L,'native_viscosity_unit_Pa_s':eta0,'native_radius':an,'native_viscosity':etan,'native_six_unit_loads':True,'SI_unit_impulse_equivalence':'Native column multiplied by analytical load/output unit conversion; linear solve','native_scaled_resistance_condition':cond,'debye_cut':1e-5,'debye_cut_active':h/a<1e-5,'wall_cutoff_argument':0.,'wall_cutoff_semantics':'upstream skips height<wall_cutoff; zero includes all positive heights','COO_absolute_drop_threshold_native':1e-12,'excess_resistance_native':excess.tolist(),'bulk_resistance_native':bulk.tolist(),'native_total_resistance':resistance.tolist(),'native_mobility':Mnat.tolist(),'upstream_rotational_excess_clip':'max(Rrotation-Rbulk,0), unmodified','sign_correction_applied':False,'matrix_symmetrized':False,'all_blocks_available':True}

def get_wall_mobility(radius_m,gap_m,viscosity_pa_s,implementation='sphere_semianalytical'):
    a,h,eta=map(float,(radius_m,gap_m,viscosity_pa_s))
    if min(a,h,eta)<=0:raise ValueError('Strictly positive physical inputs required')
    if implementation=='sphere_semianalytical':return _get_sphere(a,h,eta)
    if implementation=='lubrication':return _get_lub(a,h,eta)
    raise ValueError(implementation)
