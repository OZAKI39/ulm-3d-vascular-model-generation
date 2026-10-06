from common import *
import subprocess,shutil,csv,collections,xml.etree.ElementTree as ET
CASES={'pipe':'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu','BraVa':'cases/balanced_pressure_final'}

def main():
 rt=HERE/'runtime';audit={};inputs=read(HERE/'INPUT_SHA256.json');rows=[]
 for name,rel in CASES.items():
  original=BASE/rel;xml=original/'run/solver.xml';root=ET.parse(xml).getroot();eq=root.find('Add_equation');configs=[]
  case=rt/name;run=case/'run';run.mkdir(parents=True,exist_ok=True)
  if not (case/'SV_MESH').exists():(case/'SV_MESH').symlink_to(original/'SV_MESH',target_is_directory=True)
  shutil.copyfile(xml,run/'solver.xml');inputs[str(xml.relative_to(BASE))]=sha(xml)
  for bc in eq.findall('Add_BC'):
   fields={e.tag:e.text.strip() for e in bc};record=dict(face=bc.attrib['name'],xml_type=fields['Type'],fields=fields,xml=ET.tostring(bc,encoding='unicode'),external_files={})
   for key,value in fields.items():
    if key.endswith('_file_path'):
     src=xml.parent/value;dst=run/value;shutil.copyfile(src,dst);record['external_files'][key]={'path':str(src.relative_to(BASE)),'sha256':sha(src)};inputs[str(src.relative_to(BASE))]=sha(src)
   configs.append(record)
  audit[name]=dict(solver_xml=str(xml.relative_to(BASE)),sha256=sha(xml),equation_attributes=eq.attrib,equation_scalar_fields={e.tag:e.text.strip() for e in eq if not list(e) and e.text},BCs=configs)
  put('TARGET_BC_CONFIG_AUDIT.json',audit);put('INPUT_SHA256.json',inputs)
  for variant in ['original','trace']:
   prefix=rt/(name+'_'+variant)
   with (rt/(name+'_'+variant+'.log')).open('w') as log:
    subprocess.run([str(rt/('dispatch_'+variant)),'solver.xml',str(prefix)],cwd=run,env=native_env(),stdout=log,stderr=subprocess.STDOUT,check=True)
  assert (rt/(name+'_original_mixed_checks.csv')).read_bytes()==(rt/(name+'_trace_mixed_checks.csv')).read_bytes(),'observation changed mixed output'
  shutil.copyfile(rt/(name+'_trace_mixed_checks.csv'),HERE/(name+'_MIXED_BOUNDARY_CHECKS.csv'))
  assert (rt/(name+'_original_properties.csv')).read_bytes()==(rt/(name+'_trace_properties.csv')).read_bytes(),'observation changed output'
  props=list(csv.DictReader((rt/(name+'_trace_properties.csv')).open()));trace=list(csv.DictReader((rt/(name+'_trace_trace.csv')).open()))
  shutil.copyfile(rt/(name+'_trace_properties.csv'),HERE/(name+'_NATIVE_PROPERTIES.csv'));shutil.copyfile(rt/(name+'_trace_trace.csv'),HERE/(name+'_NATIVE_DISPATCH_TRACE.csv'))
  for cfg in configs:
   prop=next(p for p in props if p['face']==cfg['face']);stages={}
   for stage in ['actual_P1','mixed_P2_replay']:
    count=collections.Counter(t['function'] for t in trace if t['stage']==stage and t['face']==cfg['face'])
    # nn::get_nnx rows encode eNoN instead of face; these target sweeps must have none.
    inverses=collections.Counter(t['face'] for t in trace if t['stage']==stage and t['function']=='get_nnx')
    stages[stage]={'function_call_count':dict(count),'pressure_inverse_count_in_sweep':inverses.get('4',0),'velocity_inverse_count_in_sweep':inverses.get('10',0)}
   expected={'Dir':'set_bc_dir_l','Neu':'set_bc_neu_l','Traction':'set_bc_trac_l'}[cfg['xml_type']]
   match=all(v['function_call_count'].get(expected)==1 and not v['function_call_count'].get('set_bc_dir_wl',0) and not v['pressure_inverse_count_in_sweep'] for v in stages.values()) and prop['weakDir']=='0'
   rows.append(dict(case=name,face=cfg['face'],xml_type=cfg['xml_type'],native_properties=prop,stages=stages,static_expected=expected,static_runtime_agree=match,weak_dir_entered=any(v['function_call_count'].get('set_bc_dir_wl',0)>0 for v in stages.values()),common_volume_geometry_used_in_boundary=False,common_volume_geometry_scope='Not applicable: no volume assembly called; previously validated patched fluid object linked',result='PASS' if match else 'FAIL'))
 result=dict(status='BC_DISPATCH_ESTABLISHED' if all(x['static_runtime_agree'] for x in rows) else 'BC_DISPATCH_NOT_ESTABLISHED',rows=rows,pipe_status='PIPE_TARGET_BC_SET_CLEAR_OF_WEAK_DIRICHLET' if all(not x['weak_dir_entered'] for x in rows if x['case']=='pipe') else 'PIPE_WEAK_PATH_PRESENT',BraVa_status='BRAVA_TARGET_BC_SET_CLEAR_OF_WEAK_DIRICHLET' if all(not x['weak_dir_entered'] for x in rows if x['case']=='BraVa') else 'BRAVA_FUTURE_WEAK_DIRICHLET_BLOCKER_PRESENT',CFD_steps=0,linear_solves=0,nonlinear_solves=0,parser='Unchanged native Parameters/read_files/read_bc + distribute using byte-identical XML and unchanged full original mesh/data; face_ini and bc_ini executed',mixed_replay='Native parsed objects retain all flags; synthetic positive curved TET10/TRI6 geometry, first six parsed profile values, no new BC types; isolates order-dependent dispatch; not a full P2 BraVa mesh or boundary accuracy test',observation_output_bitwise_unchanged=True)
 put('TARGET_BC_RUNTIME_DISPATCH.json',result);print(result['status'])
if __name__=='__main__':main()
