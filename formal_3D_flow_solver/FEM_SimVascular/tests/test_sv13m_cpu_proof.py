from sv13m_support import *
from sv_validation.sv13m import *
def test_matching_cpu_proof():
 d=accepted('CPU_PROOF_20_1R_acceptance');flow_gate(d,gpu=False,steps=20)
 gpu=accepted('GPU_PROOF_20_acceptance');assert d['xml_sha256']==gpu['xml_sha256']
 assert d['command'][-2]==gpu['command'][-2]
