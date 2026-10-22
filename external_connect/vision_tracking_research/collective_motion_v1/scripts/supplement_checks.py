"""Additional explicit stress/transport checks, no threshold changes or model selection."""
import sys,json,time,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
from distribution import CandidateDetector
from methods import piv,cpd
from affine import parameters,empty
from state import *
import cv2

z=np.load(ROOT/'results/synthetic/test/7202_mixed.npz');base=z['frames'];H=z['H'];det=CandidateDetector();rows=[]
for case in ['local_defocus','extra_dark_distractors']:
    ims=base.copy()
    for i,im in enumerate(ims):
        if case=='local_defocus':ims[i,:,:128]=cv2.GaussianBlur(im,(0,0),5)[:,:128]
        else:
            rng=np.random.default_rng(1943+i)
            if i>=8:
                for point in rng.integers(15,240,(15,2)):cv2.circle(ims[i],tuple(point),3,85,-1)
    np.savez_compressed(ROOT/f'results/supplement_{case}.npz',frames=ims)
    for n in P['windows']:
        end=15;begin=end-n+1;dt=(n-1)*.0333;hh=H[end]@np.linalg.inv(H[begin]);truth=parameters(hh[:2,:2],hh[:2,2],dt,np.array([128.,128.]));oa,ga=det(ims[begin]);ob,gb=det(ims[end])
        for method,r in [('PIV',piv(ga,gb,dt,'local',oa,ob,multipass=True)),('CPD',cpd([o['xy'] for o in oa],[o['xy'] for o in ob],dt,ga.shape))]:
            r.pop('vectors',None);rows.append({'case':case,'window':n,'truth':truth,'method':method,'detected_counts':[len(oa),len(ob)],**r})
dump(ROOT/'results/supplement_stress.json',{'scope':'post-benchmark diagnostic only; no retuning; primary test results preserved','rows':rows})

mock=MockConsumer();received=time.monotonic();base=make_state('FAULT_INJECTION','mock_s',1,1.,received,[0,0,32,32],{'group_center_xy':[16,16],'field_validity':{'group_center_xy':True}},empty('NO_SUPPORT'),{},3)
base['motion']['values']['translation_vx']=2.;base['motion']['valid']['translation_vx']=True
mock.send(copy.deepcopy(base));time.sleep(.02)
mock.send(copy.deepcopy(base));time.sleep(.02)
bad=copy.deepcopy(base);bad['frame_id']=2;bad['source_PTS']=.5;mock.send(bad);time.sleep(.02)
bad=copy.deepcopy(base);bad['frame_id']=3;bad['source_PTS']=1.2;bad['host_receive_timestamp']=received-1;mock.send(bad);time.sleep(.15)
mock.close();reasons=[r['reason'] for r in mock.rows];assert len(mock.rows)==4 and mock.rows[0]['fresh'] and not any(r['fresh'] for r in mock.rows[1:]);assert mock.watchdog_events
dump(ROOT/'results/mock_fault_tests.json',{'PASS':True,'consumer':mock.rows,'watchdog_events':mock.watchdog_events,'null_input_not_zero':bool(np.isnan(tensor_input(mock.states[-1])[0]).all()),'hardware_executed':False})

# Formal minimal JSON Schema: runtime adds statistical diagnostics without lossy coercion.
nullable_number={'type':['number','null']};nullable_vec={'type':['array','null'],'items':{'type':'number'},'minItems':2,'maxItems':2}
motion={'type':'object','required':['values','valid','reasons','observability_status'],'properties':{
    'values':{'type':'object','required':NAMES,'properties':{k:nullable_number for k in NAMES},'additionalProperties':False},
    'valid':{'type':'object','required':NAMES,'properties':{k:{'type':'boolean'} for k in NAMES},'additionalProperties':False},
    'reasons':{'type':'object','required':NAMES,'properties':{k:{'type':['string','null']} for k in NAMES}},
    'observability_status':{'enum':['NOT_OBSERVABLE','PARTIAL','INTERVAL_AFFINE_ESTIMATE']}}}
for k in NAMES:
    motion.setdefault('allOf',[]).append({'if':{'properties':{'valid':{'properties':{k:{'const':False}}}}},'then':{'properties':{'values':{'properties':{k:{'type':'null'}}}}}})
action={'type':['object','null'],'properties':{'pzt_channel':{'enum':['X-','X+','Y-','Y+',None]},'actual_action_state':{'type':['string','null']},'amplitude':nullable_number,'amplitude_unit':{'type':['string','null']},'phase_deg':nullable_number,'frequency_Hz':nullable_number,'command_time':nullable_number,'execution_time':nullable_number,'time_sync_quality':{'type':['string','null']},'acoustic_pressure_calibrated':{'type':'boolean'}}}
schema={'$schema':'https://json-schema.org/draft/2020-12/schema','title':'M4 collective observation; no device actions','type':'object','required':['schema_version','source','session_id','frame_id','available_capture_timestamp','host_receive_timestamp','processing_finished_timestamp','source_PTS','ROI','calibration_status','image_quality','distribution','motion','time_window','age_ms','stale_or_discontinuous','action_event','hardware_control_executed','training_ready'],
 'properties':{'schema_version':{'const':'M4-collective-1.0'},'source':{'type':'string'},'session_id':{'type':'string'},'frame_id':{'type':'integer'},'available_capture_timestamp':nullable_number,'source_PTS':nullable_number,'host_receive_timestamp':{'type':'number'},'processing_finished_timestamp':{'type':'number'},'ROI':{'type':'array','minItems':4,'maxItems':4,'items':{'type':'integer'}},'distribution':{'type':'object','properties':{'group_center_xy':nullable_vec,'distribution_orientation':nullable_number}},'motion':motion,'action_event':action,'hardware_control_executed':{'const':False},'training_ready':{'const':False}}}
dump(ROOT/'configs/collective_state.schema.json',schema)
dump(ROOT/'evidence/state_example.json',json.loads((ROOT/'results/causal_replay_fast.json').read_text())['records'][10]['state'])
print('supplement and mock checks complete')
