from sv13n_support import *
def test_native_safe_stop_after_first_complete_rule():
 d=accepted('gpu_steady_candidate');r=read('gpu_steady_stop_request')
 assert r['kind']=='STEADY' and r['detail']['first_full_steady_step']==d['first_full_steady_step']
 assert 'STOP_SIM=0' in r['mechanism'] and d['safe_stop'] and d['normal_exit']
 assert d['stop_step']>=d['first_full_steady_step'] and d['checkpoint']['step']==d['stop_step']
 assert not d['production_changed']
def test_incomplete_and_nonfinite_checkpoint_rejected(tmp_path):
 import struct
 data=struct.pack('<8i3d',1,1,1,2,0,4,0,10,10e-7,0.,1.)+struct.pack('<16d',*range(16))
 p=tmp_path/'state.bin';p.write_bytes(data)
 assert checkpoint_one_rank(p,10,1e-7)['complete_integration_history']
 p.write_bytes(data[:-1])
 with pytest.raises(GateError):checkpoint_one_rank(p,10,1e-7)
 p.write_bytes(data[:-8]+struct.pack('<d',float('nan')))
 with pytest.raises(GateError):checkpoint_one_rank(p,10,1e-7)
