"""Secondary diagnostic against visually checked O'Neill(1964) p73 F*,G*.
Eight paper alpha points; these are additional reference evaluations, never a
retry or replacement of the main 14-gap sweep. Original rotation F,G excluded
because Goldman et al. later corrected the Dean/O'Neill rotation calculations.
"""
from pathlib import Path
import json,csv,sys
import numpy as np
R=Path(__file__).resolve().parents[1];label=sys.argv[1]
if label=='RMBW_LUBRICATION':
 from rmbw_reference_adapter import get_wall_mobility as get
 def mobility(a,h,mu):return get(a,h,mu,'lubrication')
else:
 from pystokes_reference_adapter import get_wall_mobility as mobility
paper=[(3.,10.0677,1.0591,.00001),(2.,3.7622,1.1738,.00042),(1.,1.5431,1.5675,.01465),(.5,1.1276,2.1515,.07372),(.3,1.0453,2.6475,.14552),(.1,1.0050,3.7863,.34187),(.08,1.0032,4.0223,.38496),(.06,1.0018,4.3275,.44116)]
a=json.loads((R/'WALL_REFERENCE_CONVENTION.json').read_text())['radii_m']['d50'];mu=.001;out=[]
p=R/'raw'/(label+'_ONEILL_TABLE.json');assert not p.exists(),'No repeat of the first historical diagnostic'
for alpha,rounded_height,fs,gs in paper:
 e=float(np.cosh(alpha)-1);M,meta=mobility(a,a*e,mu);D=np.diag(np.sqrt([6*np.pi*mu*a]*3+[8*np.pi*mu*a**3]*3));B=D@M@D;C=np.linalg.solve(B,np.eye(6));Rsi=D@C@D
 # R maps imposed velocities to external holding loads; O'Neill uses fluid loads.
 fn=float(Rsi[0,0]/(6*np.pi*mu*a));gn=float(-Rsi[4,0]/(8*np.pi*mu*a*a))
 out.append({'implementation':label,'alpha':alpha,'epsilon':e,'printed_height_over_a':rounded_height,'height_method':'cosh(alpha), not rounded printed d/a','paper_Fstar':fs,'paper_Gstar':gs,'package_Fstar':fn,'package_Gstar':gn,'Fstar_relative_difference':abs(fn/fs-1),'Gstar_relative_difference':abs(gn/gs-1),'matrix_SI':M.tolist(),'metadata':meta,'scope':'historical independent translational resistance and induced torque; printed G* has low significant digits at far wall','source_pdf_page':7,'source_url':'https://doi.org/10.1112/S0025579300003508'})
p.write_text(json.dumps(out,indent=2)+'\n');print(label,[(x['epsilon'],x['Fstar_relative_difference'],x['Gstar_relative_difference']) for x in out])
