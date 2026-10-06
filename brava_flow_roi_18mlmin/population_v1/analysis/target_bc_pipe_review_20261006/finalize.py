from common import *
import csv,re,datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import static_audit

def main():
 static_audit.main()
 runtime=read(HERE/'TARGET_BC_RUNTIME_DISPATCH.json');cfg=read(HERE/'TARGET_BC_CONFIG_AUDIT.json');reg=read(HERE/'REGRESSION_RESULTS.json');mesh=read(HERE/'P2_CURVED_MESH_QA.json');resource=read(HERE/'P2_RESOURCE_PREFLIGHT.json')
 for name,c in cfg.items():
  for bc in c['BCs']:
   r=next(r for r in runtime['rows'] if r['case']==name and r['face']==bc['face'])
   bc['native_resolved']=dict(bType=int(r['native_properties']['bType']),weakDir=False,Weakly_applied_XML_present='Weakly_applied' in bc['fields'],Weakly_applied_default=False,General_profile_setup_early_return=(bc['fields'].get('Time_dependence')=='General' or bc['xml_type']=='Traction'))
 put('TARGET_BC_CONFIG_AUDIT.json',cfg)
 # Expose exact enum names as well as integer bitmasks.
 enum=(SOURCE/'consts.h').read_text().split('enum class BoundaryConditionType',1)[1].split('};',1)[0]
 flags={int(i):name for name,i in re.findall(r'(bType_\w+)\s*=\s*(\d+)',enum)}
 tests=list(reg['tests'])
 def check(name,condition,**extra):tests.append(dict(name=name,status='PASS' if condition else 'FAIL',**extra))
 for row in runtime['rows']:
  mask=int(row['native_properties']['bType']);row['internal_type_flags']=[name for bit,name in flags.items() if mask&(1<<bit)]
  check(row['case']+'_'+row['face']+'_native_dispatch_matches_static',row['static_runtime_agree'])
  check(row['case']+'_'+row['face']+'_no_pressure_boundary_inverse',all(x['pressure_inverse_count_in_sweep']==0 for x in row['stages'].values()))
  check(row['case']+'_'+row['face']+'_pressure_slots_unchanged',float(row['native_properties']['pressure_changed'])==0.)
 for name in cfg:
  with (HERE/(name+'_MIXED_BOUNDARY_CHECKS.csv')).open() as f:
   for row in csv.DictReader(f):
    for k,v in row.items():
     if k!='face':check(name+'_'+row['face']+'_'+k,float(v)<=1e-11,value=float(v),limit=1e-11)
  inlet=next(r for r in runtime['rows'] if r['case']==name and r['face']=='INLET')['native_properties']
  if name=='BraVa':
   check('BraVa_native_impose_flux_profile_integral_one',abs(float(inlet['flux_integral'])-1)<1e-11)
   check('BraVa_native_zero_out_perimeter_22_nodes',int(inlet['gx_zero'])==22)
  else:check('pipe_General_keeps_nodal_provider_no_gx_normalization',int(inlet['general_dof'])==3 and int(inlet['gx_size'])==0)
 check('trace_does_not_modify_original_numerical_outputs',runtime['observation_output_bitwise_unchanged'])
 for k,v in mesh['checks'].items():check('mesh_'+k,v)
 check('resource_graph_P1_reproduces_native_measured_NNZ',resource['P1_scalar_structural_nnz']==7719184)
 check('resource_decision_follows_fixed_route_no_OOM_rule',resource['status']=='PREDICTED_MEMORY_EXCEEDS_CGROUP')
 runtime['weak_dispatch_definition']='set_bc_dir_w wrapper is called and skips all eight target boundaries; zero set_bc_dir_wl calls. Wrapper entry alone does not apply weak Dirichlet.'
 put('TARGET_BC_RUNTIME_DISPATCH.json',runtime)
 readiness={
 'volume_common_geometry_validated_unchanged':reg['inherited']['COMMON_GEOMETRY_HARNESS_RESULTS.json']['status']=='PASS',
 'P1_P1_regression_unchanged':reg['inherited']['P1_REGRESSION_RESULTS.json']['status']=='PASS',
 'pipe_runtime_dispatch_established':runtime['status']=='BC_DISPATCH_ESTABLISHED',
 'pipe_clear_of_unfixed_weak_path':runtime['pipe_status']=='PIPE_TARGET_BC_SET_CLEAR_OF_WEAK_DIRICHLET',
 'strong_Dirichlet_regression_pass':all(t['status']=='PASS' for t in tests if 'strong' in t['name']),
 'prescribed_Cauchy_regression_pass':reg['inherited']['boundary']['engineering_passed']==68,
 'complete_curved_mesh_QA_pass':mesh['status']=='PASS',
 'CPU_sampler_regression_pass':reg['status']=='PASS' and reg['prior_CPU_sampler_tests']=='PASS',
 'production_P1_baseline_hash_unchanged':all(t['status']=='PASS' for t in reg['tests'] if 'P1_' in t['name']),
 'fixed_GPU_route_resource_preflight_pass':resource['status']!='PREDICTED_MEMORY_EXCEEDS_CGROUP'}
 gate=dict(status='READY_FOR_ONE_P2P1_PIPE_GPU_CFD' if all(readiness.values()) else 'NOT_READY_FOR_P2P1_PIPE_GPU_CFD',checks=readiness,nine_implementation_prerequisites_pass=all(v for k,v in readiness.items() if k!='fixed_GPU_route_resource_preflight_pass'),stop_reason='PREDICTED_MEMORY_EXCEEDS_CGROUP',authority='Current user section9.3: if predicted resources exceed cgroup stop before launch; do not OOM-test',weak_boundary_patch_needed=False,CFD_runs=0,science_classification='P2_PIPE_NOT_EVALUABLE')
 put('P2_PIPE_READINESS_GATE.json',gate)
 # Dispatch figure and exact table source.
 columns=['Case / face','XML / internal flags','Actual native leaf path','Mode','Weak calls\nP1 / P2']
 display=[];csvrows=[]
 order={'INLET':0,'WALL':1,'OUTLET_01':2,'OUTLET_02':3,'OUTLET_03':4}
 sortedrows=sorted(runtime['rows'],key=lambda r:(r['case']!='pipe',order[r['face']]))
 for r in sortedrows:
  mode={'Dir':'Strong Dirichlet','Neu':'Natural pressure','Traction':'Prescribed Cauchy'}[r['xml_type']]
  flags_text='/'.join(f.removeprefix('bType_') for f in r['internal_type_flags'])
  native=r['static_expected']
  if r['xml_type']=='Neu':native+='\n→ b_assem_neu_bc'
  if r['xml_type']=='Dir':native='set_bc_dir\n→ '+native
  else:native='set_bc_neu\n→ '+native
  display.append([r['case']+' '+r['face'],r['xml_type']+' | '+str(r['native_properties']['bType'])+'\n'+flags_text,native,mode,'0 / 0'])
  csvrows.append(dict(case=r['case'],face=r['face'],xml_type=r['xml_type'],internal_mask=r['native_properties']['bType'],internal_flags=flags_text,native_leaf=r['static_expected'],weak_P1=0,weak_P2=0,mode=mode))
 with (HERE/'TARGET_BC_DISPATCH_FIGURE.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(csvrows[0]));w.writeheader();w.writerows(csvrows)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11})
 fig,ax=plt.subplots(figsize=(16,8.4));ax.axis('off')
 table=ax.table(cellText=display,colLabels=columns,colWidths=[.17,.26,.29,.18,.10],cellLoc='left',bbox=[0,.08,1,.82]);table.auto_set_font_size(False);table.set_fontsize(10.5)
 for (r,c),cell in table.get_celld().items():
  cell.set_edgecolor('#cbd5e1')
  if r==0:cell.set_facecolor('#18334a');cell.set_text_props(color='white',weight='bold')
  elif c==4:cell.set_facecolor('#d6f0e6');cell.set_text_props(weight='bold',color='#185743')
  else:cell.set_facecolor('#f1f5f9' if r%2==0 else '#ffffff')
 ax.set_title('Target boundary dispatch: neither case enters weak Dirichlet assembly',loc='left',fontsize=18,pad=18)
 fig.text(.04,.055,'Full original XML / mesh / data → native parser + setup; mixed replay retains parsed flags.\nset_bc_dir_w is visited as a dispatcher; set_bc_dir_wl and pressure inverse-map call counts are zero.',fontsize=11,color='#334155')
 fig.subplots_adjust(left=.04,right=.97,top=.9,bottom=.12);fig.savefig(HERE/'01_TARGET_BC_DISPATCH.png',dpi=300);plt.close(fig)
 # Resource figure deliberately is not a fabricated actual-field comparison.
 models=resource['models'];labels=['RSS × structural NNZ','RSS × independent DOF','LU arrays; no overlap*','LU arrays; observed ASM overlap'];vals=[m['GiB'] for m in models]
 fig,(ax,bx)=plt.subplots(1,2,figsize=(15,5.5),gridspec_kw={'width_ratios':[1.5,1]})
 ax.barh(labels,vals,color=['#487ca8','#487ca8','#adb5bd','#d68748']);ax.axvline(resource['remaining_to_cgroup_limit_GiB'],color='#b73535',ls='--',label='Available to cgroup limit');ax.axvline(resource['cgroup_limit_GiB'],color='#111827',ls=':',label='Total cgroup limit')
 for i,v in enumerate(vals):ax.text(v+1,i,f'{v:.1f} GiB',va='center',fontsize=10)
 ax.set_xlim(0,max(vals)*1.17);ax.set_xlabel('Planning estimate, not measured P2 memory');ax.invert_yaxis();ax.legend(loc='lower right',fontsize=9)
 arr=np.load(HERE/'P2_MESH_QA_SAMPLES.npz');bx.hist(arr['perturbation_bound'],bins=50,color='#368475');bx.set_yscale('log');bx.set_xlabel('Whole-element Jacobian perturbation bound');bx.set_ylabel('TET10 cells');bx.set_title('174,785 / 174,785 certified positive',fontsize=12)
 bx.text(.40,.92,f'Max bound {mesh["max_perturbation_bound"]:.5f} < 1\nNo corners moved; no refinement',transform=bx.transAxes,va='top',fontsize=10,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
 fig.suptitle('Implementation gates pass; unchanged GPU route fails resource preflight',fontsize=16)
 fig.text(.03,.015,'* No-overlap estimate is sensitivity only: the frozen route uses ASM overlap 2. All values are conditional estimates. No CFD launched.',fontsize=10)
 fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(HERE/'02_READINESS_RESOURCE_AUDIT.png',dpi=300);plt.close(fig)
 with (HERE/'RESOURCE_ESTIMATE_FIGURE.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['model','GiB','role']);w.writerows((m['name'],m['GiB'],'sensitivity' if 'no overlap' in m['name'] else 'planning') for m in models)
 for f in HERE.glob('*.png'):
  with Image.open(f) as im:check(f.name+'_300dpi',all(abs(x-300)<.1 for x in im.info['dpi']),dpi=list(im.info['dpi']))
 put('TEST_RESULTS.json',dict(status='PASS' if all(t['status']=='PASS' for t in tests) else 'FAIL',passed=sum(t['status']=='PASS' for t in tests),failed=sum(t['status']=='FAIL' for t in tests),tests=tests,inherited_checks=reg['inherited'],scientific_gates='NOT_RUN: no actual P2 CFD field',readiness=gate['status'],unused_weak_library_defect_preserved=True,tests_are_not_physical_certification=True))
 put('RESULT_IDENTITY.json',dict(status='P2_PIPE_NOT_EVALUABLE',NEXT_RECOMMENDATION='P2P1_PIPE_RESULT_NOT_EVALUABLE',readiness=gate['status'],stop_reason=gate['stop_reason'],pipe_BC=runtime['pipe_status'],BraVa_BC=runtime['BraVa_status'],scope='Target BC audit + single deterministic curved elevation and prelaunch assessment',local_start_HEAD='c802f0f7a569b9f4b7cc63a1d96d4d91c181be0e',current_source_commit=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json').get('source_commit'),common_patch_identity_sha256=sha(PREVIOUS/'COMMON_GEOMETRY_PATCH_IDENTITY.json'),dispatch_build_identity_sha256=sha(HERE/'DISPATCH_BUILD_IDENTITY.json'),curved_mesh_sha256=mesh['curved_mesh_sha256'],production_modified=False,CFD_runs=0,BraVa_CFD=0,particle_runs=0,population_runs=0,weak_patch=False,full_solver_binary_built=False,high_order_CFD_output='NOT_RUN',P2_actual_field_metrics='NOT_RUN',automatic_followup=False))
 (HERE/'OPEN_RESULTS.html').write_text('''<!doctype html><meta charset="utf-8"><title>Target BC audit and P2 pipe readiness</title><style>body{font:18px system-ui;max-width:1400px;margin:32px auto;padding:0 24px;color:#18334a}img{width:100%}a{color:#176694}</style><h1>Target BC audit and P2 pipe readiness</h1><p>Pipe and BraVa: no weak Dirichlet assembly. Full curved mesh QA passed. CFD not launched: resource preflight stopped the frozen GPU route.</p><p><a href="P2_TARGET_BC_AND_PIPE_CFD_REVIEW_ZH.md">中文主报告</a> · <a href="P2_PIPE_READINESS_GATE.json">Readiness gate</a> · <a href="P2_RESOURCE_PREFLIGHT.json">Resource estimate and assumptions</a></p><img src="01_TARGET_BC_DISPATCH.png"><img src="02_READINESS_RESOURCE_AUDIT.png"><p>These are dispatch and readiness figures. No actual P2 CFD field or trajectory result is claimed.</p>''')
 print('Final engineering tests',len(tests),'failed',sum(t['status']=='FAIL' for t in tests))
if __name__=='__main__':main()
