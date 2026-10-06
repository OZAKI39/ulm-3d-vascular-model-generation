from common import *
import importlib.util
import numpy as np

def main():
 tests=[]
 def check(n,v,**extra):tests.append(dict(name=n,status='PASS' if v else 'FAIL',**extra))
 for directory in [PREVIOUS,NATIVE]:
  manifest=read(directory/'OUTPUT_SHA256.json');check(directory.name+'_unchanged_artifacts',all(sha(directory/p)==h for p,h in manifest.items()),files=len(manifest))
 inherited={}
 for filename in ['COMMON_GEOMETRY_HARNESS_RESULTS.json','P1_REGRESSION_RESULTS.json']:
  v=read(PREVIOUS/filename);check(filename+'_prior_status',v['status']=='PASS');inherited[filename]={'sha256':sha(PREVIOUS/filename),'status':v['status'],'tests':len(v['tests']),'evidence_mode':'Unchanged validated input/source/object/result hashes; not counted as newly rerun tests'}
 boundaries=read(PREVIOUS/'BOUNDARY_MIXED_GEOMETRY_AUDIT.json');check('inherited_strong_and_Cauchy_68_checks',all(t['status']=='PASS' for t in boundaries['tests']));inherited['boundary']={'engineering_passed':len(boundaries['tests']),'unused_weak_gate_failures_preserved':boundaries['scientific_gate_failed'],'sha256':sha(PREVIOUS/'BOUNDARY_MIXED_GEOMETRY_AUDIT.json')}
 patched=read(PREVIOUS/'COMMON_GEOMETRY_PATCH_IDENTITY.json');check('validated_patched_object_SHA',sha(PREVIOUS/'runtime/patched.o')==patched['objects']['runtime/patched.o']['object_sha256'])
 check('validated_patched_fluid_source_SHA',sha(PREVIOUS/'runtime/fluid_patched.cpp')==patched['objects']['runtime/patched.o']['source_sha256'])
 # Minimal current CPU regression, not a new representation campaign.
 spec=importlib.util.spec_from_file_location('sampler',NATIVE/'p2_sampler.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 nodes=np.vstack([np.eye(3),np.zeros(3)]);X,C,_=mod.elevate_affine(nodes,np.array([[0,1,2,3]]));u=X**2;p=X[:,2]
 sample=mod.P2Sampler(X,C,u,p);pts=np.array([[.1,.2,.3],[.25,.25,.25],[.5,.1,.1]])
 got=sample.sample(pts);exactG=np.array([np.diag(2*q) for q in pts]);err=max(float(abs(got['u']-pts**2).max()),float(abs(got['G']-exactG).max()))
 check('CPU_quadratic_sampler_minimal_regression',got['inside_lumen'].all() and err<1e-11,error=err,threshold=1e-11)
 old=np.load(NATIVE/'P2_CURVED_GEOMETRY_SAMPLES.npz');new=np.load(HERE/'runtime/PIPE_CURVED_ELEVATION.npz')
 differences={}
 for label,cell in [('central_wall',143654),('inlet_rim',35265)]:
  a=old[label+'_curved_nodes'];b=new['points_m'][new['tetra10'][cell]]
  if a.shape!=b.shape:a=a.T
  delta=float(abs(a-b).max());differences[label]=delta;check(label+'_retained_geometry_reproduced',delta<1e-12,max_abs_m=delta)
 # All known original input / P1 baseline results rechecked, including final VTU/BIN.
 inputs=read(HERE/'INPUT_SHA256.json');pw=BASE/'production_workflow'
 for key,h in read(pw/'FIELD_VALIDATION.json')['source_sha256'].items():
  path=pw/key;check('P1_field_SHA:'+key,sha(path)==h);inputs[str(path.relative_to(BASE))]=h
 identity=read(pw/'GPU_MANUFACTURED_RUN_IDENTITY.json');saved=identity['steady_saved_state']
 for key in ['final_VTU','final_native_checkpoint']:
  path=pw/saved[key];h=saved[key+'_sha256'];check('P1_'+key+'_unchanged',sha(path)==h);inputs[str(path.relative_to(BASE))]=h
 check('all_inherited_inputs_unchanged',all(sha(BASE/p)==h for p,h in inputs.items()),files=len(inputs));put('INPUT_SHA256.json',inputs)
 put('REGRESSION_RESULTS.json',dict(status='PASS' if all(t['status']=='PASS' for t in tests) else 'FAIL',tests=tests,inherited=inherited,prior_CPU_sampler_tests=read(NATIVE/'P2_FIELD_SAMPLER_TESTS.json')['status'],prior_Poiseuille_representation='PASS_IN_AFFINE_GEOMETRY; reused only, not assumed exact on curved P2',curved_reference_differences_m=differences))
 print('Regression',len(tests),'checks',sum(t['status']=='FAIL' for t in tests),'failed')
if __name__=='__main__':main()
