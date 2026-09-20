from sv13n_support import *
def test_real_vascular_stability():
 d=accepted('REAL_VASCULAR_GPU_acceptance');flow_gate(d,gpu=True)
 assert d['mpi_ranks']==d['OMP_NUM_THREADS']==d['GPU_count']==1 and d['flows_finite'] and d['wall_noslip_pass']
 assert {r['step'] for r in d['history']['linear_solves']}==set(range(1,d['steps_completed']+1))
def test_scientific_configuration_preserved():
 import xml.etree.ElementTree as ET
 paths=[ROOT/'outputs/sv1_3n'/n/'solver.xml' for n in ('OLD_PETSC_CPU_PROOF_20','REAL_VASCULAR_GPU')]
 roots=[ET.parse(p).getroot() for p in paths]
 for r in roots:
  r.find('.//Number_of_time_steps').text='CAP'
  r.find('.//Increment_in_saving_VTK_files').text='IO_ONLY'
 def canonical(e):return e.tag,tuple(sorted(e.attrib.items())),(e.text or '').strip(),tuple(map(canonical,e))
 assert canonical(roots[0])==canonical(roots[1])
