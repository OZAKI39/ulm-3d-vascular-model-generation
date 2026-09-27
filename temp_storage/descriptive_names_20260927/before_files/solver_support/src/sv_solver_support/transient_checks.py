"""SV1.2 frozen production policy, native checkpoint audit and acceptance."""
import gzip,json,math,struct,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from sv_validation.provenance import sha256,write_json
from sv_validation.validation import require,steady_state
ROOT=Path(__file__).resolve().parents[3]
REPORT=ROOT/'reports/sv1_2';OUTPUT=ROOT/'outputs/sv1_2';CONFIG=ROOT/'configs/sv1_2';LOG=ROOT/'logs/sv1_2'
def load(name,stage='sv1_2'):
    return json.loads((ROOT/'reports'/stage/(name+'.json')).read_text())
def checkpoint_audit(path,step,dt,ranks=4):
    path=Path(path);require(path.suffix=='.bin','VTU is not native restart state')
    data=path.read_bytes();prefix=struct.pack('<3i',ranks,1,1);offsets=[];start=0
    while True:
        pos=data.find(prefix,start)
        if pos<0:break
        if pos+56<=len(data):
            head=struct.unpack_from('<8i3d',data,pos)
            if head[4:8]==(0,4,0,step) and math.isclose(head[8],step*dt,rel_tol=1e-12):offsets.append(pos)
        start=pos+1
    require(len(offsets)==ranks and offsets[0]==0,'Incomplete native MPI checkpoint headers')
    stride=offsets[1];require(offsets==[i*stride for i in range(ranks)],'Invalid MPI record spacing')
    records=[]
    for rank,offset in enumerate(offsets):
        header=struct.unpack_from('<8i3d',data,offset);nodes=header[3];end=offset+56+64*nodes
        require(end<=len(data),'Truncated native integration history')
        values=np.frombuffer(data,dtype='<f8',count=8*nodes,offset=offset+56)
        require(np.isfinite(values).all() and np.isfinite(header[8:]).all(),'Nonfinite checkpoint state')
        records.append({'rank':rank,'offset':offset,'nodes':nodes,'step':header[7],'time_s':header[8],
                        'stored_timer':header[9],'equation_initial_norm':header[10],
                        'Y_values':4*nodes,'A_values':4*nodes,'finite':True})
    require(stride==56+64*max(r['nodes'] for r in records),'Native record length mismatch')
    require(end==len(data),'Unexpected trailing checkpoint bytes')
    return {'status':'PASS','path':str(path),'sha256':sha256(path),'size_bytes':len(data),'records':records,
            'complete_integration_history':True,'fields':['velocity and pressure Y_n','time derivative A_n','step','physical time','eq.iNorm'],
            'format':'pinned native write_restart/read_restart: 7-int stamp, int step, 3 doubles, Y, A; no displacement in this rigid fluid case'}

def canonical_xml(path,allow_restart=False,allow_extension=False):
    root=ET.parse(path).getroot()
    if allow_restart:root.find('.//Continue_previous_simulation').text='false'
    if allow_extension:root.find('.//Number_of_time_steps').text='400'
    def canon(e):return (e.tag,tuple(sorted(e.attrib.items())),(e.text or '').strip(),tuple(canon(c) for c in e))
    return canon(root)
def validate_xml(path,extension=False):
    original=ROOT/'configs/sv1_1/sv_flow_petsc.xml'
    require(canonical_xml(original)==canonical_xml(path,True,extension),'FROZEN_INPUT_CHANGED: XML')
    root=ET.parse(path)
    require(root.findtext('.//Continue_previous_simulation').strip().lower() in ('true','false'),'Invalid restart flag')
    require(int(root.findtext('.//Number_of_time_steps'))==(800 if extension else 400),'Unauthorized step target')
    return True
def validate_frozen():
    manifest=load('frozen_input_manifest')
    for item in manifest['files']:
        path=Path(item['path']);require(path.is_file() and sha256(path)==item['sha256'],'FROZEN_INPUT_CHANGED: '+str(path))
    config=json.loads((CONFIG/'petsc_options.json').read_text())
    require(config['PETSC_OPTIONS']==load('petsc_short_execution','sv1_1')['PETSC_OPTIONS'],'FROZEN_INPUT_CHANGED: PETSc options')
    require(config['mpi_ranks']==4 and config['OMP_NUM_THREADS']==1,'FROZEN_INPUT_CHANGED: MPI/OMP')
    return True
def interval_metrics(measure,u,previous_u,current,previous,dt_steps=10):
    require(current['step']-previous['step']==dt_steps,'Saved states are not consecutive')
    return {'step':current['step'],'previous_step':previous['step'],'time_s':current['time_s'],
            'E_u':measure.velocity_l2(u-previous_u)/(measure.velocity_l2(u)+np.finfo(float).tiny),
            'E_Q':max(abs(current['signed_outward_boundary_flows_m3_s'][r]-previous['signed_outward_boundary_flows_m3_s'][r])/measure.Q for r in current['signed_outward_boundary_flows_m3_s'])}
def steady_gate(intervals):
    if len(intervals)<5:return False
    return steady_state([r['E_u'] for r in intervals],[r['E_Q'] for r in intervals])
def consecutive_passes(intervals):
    count=0
    for r in reversed(intervals):
        if not (math.isfinite(r['E_u']) and math.isfinite(r['E_Q']) and r['E_u']<=1e-5 and r['E_Q']<=1e-6):break
        count+=1
    return count
def extension_eligible(run,qc):
    if run['linear_failures'] or run['nonlinear_failures'] or run['stop'] or qc['last_step']!=400:return False
    intervals=qc['intervals']
    if len(intervals)<10 or not all(math.isfinite(s['epsilon_mass']) and s['velocity_finite'] and s['pressure_finite'] for s in qc['states']):return False
    # Stable/declining final five intervals relative to preceding five, allowing
    # threshold-level roundoff floors without disguising growth above a gate.
    return all(max(r[k] for r in intervals[-5:])<=max(limit,max(r[k] for r in intervals[-10:-5])) for k,limit in [('E_u',1e-5),('E_Q',1e-6)])
