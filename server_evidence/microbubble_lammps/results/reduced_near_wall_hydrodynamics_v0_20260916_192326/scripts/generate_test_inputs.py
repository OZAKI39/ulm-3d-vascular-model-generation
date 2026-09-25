"""Deterministic physical and invalid-input probes; no expected C++ answers."""
import csv
import json
import math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
grid=np.genfromtxt(ROOT/'raw/CF2003_REFERENCE_DENSE_GRID.csv',delimiter=',',names=True)['epsilon']
rows=[]
def add(a,mu,e,frame,g,group='extra',expected='OK',source=0,h=None):
    h=e*a if h is None else h
    # Decimal endpoints need a representable SI gap whose computed h/a is in
    # the domain. Record this <=1 ULP fixture adjustment, never clamp in kernel.
    if expected=='OK':
        if h/a<.001:h=np.nextafter(h,math.inf)
        if h/a>.2:h=np.nextafter(h,0.)
    rows.append([str(len(rows)),a,mu,h,*np.asarray(frame).T.ravel(),*g,source,expected,group])
for e in grid:add(1e-6,.001,e,np.eye(3),[100.,0.,0.],'dense')
rng=np.random.default_rng(20260916)
for i in range(96):
    q,_=np.linalg.qr(rng.normal(size=(3,3)))
    if np.linalg.det(q)<0:q[:,1]*=-1
    e=float(np.exp(rng.uniform(math.log(.001),math.log(.2))))
    a=float(10**rng.uniform(-7,-4));mu=float(10**rng.uniform(-4,-1))
    g=rng.normal(size=3)*10**rng.uniform(-3,4)
    if i==0:g=np.zeros(3)
    if i==1:g=30*q[:,2]
    for angle in [0.,.63,math.pi]:
        rotation=np.array([[math.cos(angle),-math.sin(angle),0],[math.sin(angle),math.cos(angle),0],[0,0,1]])
        add(a,mu,e,q@rotation,g,'rotated_'+str(i))
    for factor,name in [(0.,'zero'),(-1.,'reverse'),(2.,'double')]:
        add(a,mu,e,q,factor*g,name+'_'+str(i))
    add(2*a,3*mu,e,q,g,'scale_'+str(i))
def bad(**kw):
    args=dict(a=1e-6,mu=.001,e=.01,frame=np.eye(3),g=[1.,0.,0.],group='invalid')
    args.update(kw);add(**args)
for h in [0.,-1e-9]:bad(h=h,expected='INVALID_GAP')
for h in [float('nan'),float('inf'),-float('inf')]:bad(h=h,expected='INVALID_INPUT')
for a in [0.,-1e-6,float('nan'),float('inf')]:bad(a=a,h=1e-8,expected='INVALID_INPUT')
for mu in [0.,-.001,float('nan'),float('inf')]:bad(mu=mu,expected='INVALID_INPUT')
for e in [.000999,20.01]:bad(e=e,expected='OUTSIDE_RESISTANCE_DOMAIN')
for e in [.20001,1.,20.]:bad(e=e,expected='OUTSIDE_SHEAR_DOMAIN')
for source in [1,2,99]:bad(source=source,expected='UNSUPPORTED_SOURCE')
for g in [[float('nan'),0,0],[0,float('inf'),0]]:bad(g=g,expected='INVALID_INPUT')
for frame in [np.diag([1.,1.,-1.]),np.diag([1.,1.,2.]),np.zeros((3,3)),np.full((3,3),float('nan'))]:
    bad(frame=frame,expected='INVALID_FRAME')
with (ROOT/'raw/PHASE1_INPUTS.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['id','a','mu','h','t1x','t1y','t1z','t2x','t2y','t2z','nx','ny','nz','gx','gy','gz','source','expected_status','group']);w.writerows(rows)
(ROOT/'validation/INPUT_COVERAGE.json').write_text(json.dumps({'dense_reference_points':len(grid),'tilted_frame_and_transformed_cases':672,
    'invalid_cases':len(rows)-len(grid)-672,'seed':20260916,'fixture_rounding':'At epsilon endpoints only: nextafter gap by <=1 ULP when a multiplication/division rounds outside the closed domain. Kernel uses exact exterior rejection.'},indent=2)+'\n')
print('Input cases:',len(rows))
