from pathlib import Path
import shutil
S=Path(__file__).resolve().parents[1];B=S.parent/'bcflux_20260916_224835'
for name in ['INLET_SIZE_CAPACITY.json','AUTHORITATIVE_INLET_FLOW_AUDIT.json']:shutil.copy2(B/'validation'/name,S/'validation'/name)
p=S/'scripts/validate_workflow.py';s=p.read_text().replace('import numpy as np','import numpy as np\nfrom independent_wall import segment_distance,check_wall_records')
a=s.index('def safe_segment(');b=s.index('def readcase(',a)
s=s[:a]+'''def safe_segment(a,b,radius,margin=1e-10,depth=0):
    d=segment_distance(g,a,b,radius+margin)
    assert d>=radius+margin-2e-14,('STOP_DEEP_PENETRATION',d-radius)
'''+s[b:]
s=s.replace('    state=json.loads(', '    used,offset,wallstats=check_wall_records(D,rows,g)\n    state=json.loads(',1)
s=s.replace("r0=np.array([vel(r) for r in starts])","r0=np.array([vel(r) for r in starts]); u0=np.array([used[(step,0,id)] for id in ids]); us=np.array([used[(step,stage,id)] for id in ids]); off0=np.array([offset.get((step,0,id),np.zeros(3)) for id in ids]); off1=np.array([offset.get((step,1,id),np.zeros(3)) for id in ids])")
s=s.replace('base+.5*dt*r0[:,:3]','base+.5*dt*u0+off0').replace('base+dt*q[:,:3]','base+dt*us+off1')
s=s.replace("    result=dict(status='PASS'", "    (S/'events'/f'{D.name}_WALL_PARTICLE_SUMMARY.json').write_text(json.dumps(wallstats,indent=2)+'\\n')\n    result=dict(wall_constraints=wallstats,status='PASS'")
p.write_text(s)
