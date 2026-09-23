"""Read existing evidence and prepare thesis figures; no trajectory integration."""
from pathlib import Path
import json, hashlib, shutil, subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'THESIS_SUMMARY'
def load(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
stages={s:f'particle_3d/reports/particle{s}/PARTICLE{s}_VALIDATION.json' for s in ['0','1','2','3','4','5','6','6_5','7','8']}
stages.update({'8_full3d':'particle_3d/outputs/particle8_full3d/PARTICLE8_FULL3D_VALIDATION.json','8_1':'particle_3d/outputs/particle8_1/PARTICLE8_1_VALIDATION.json'})
records={k:load(v) for k,v in stages.items()}
inventory=[]
for folder in [ROOT/'particle_3d/reports',ROOT/'particle_3d/outputs',ROOT/'particle_3d/contracts']:
 for p in folder.rglob('*'):
  rel=str(p.relative_to(ROOT))
  if p.is_file() and 'particle8_2' not in rel and (p.name.endswith('REVIEW.md') or p.name.endswith('VALIDATION.json') or '/literature/' in rel or '/contracts/' in rel):
   inventory.append(dict(path=rel,bytes=p.stat().st_size,sha256=sha(p)))
(OUT/'audit/evidence_inventory.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2)+'\n')
(OUT/'audit/preexisting_git_status.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
# Exact source values used in prose/tables. This is an extraction, not a new experiment.
keys={'0':['node_count','tetra_count','max_affine_velocity_error','max_affine_pressure_error','max_affine_gradient_error','max_node_velocity_error','max_node_pressure_error'],
'1':['single_mb_diameter_um','validation_timesteps_s','exit_time_s'],
'2':['accepted_D_statistics','accepted_V_statistics','validation_population_N'],
'3':['wall_triangle_count','max_synthetic_sphere_gap_error_m','max_synthetic_ellipsoid_gap_error_m','max_synthetic_capsule_gap_error_m','rbc_volume_conservation_error'],
'5':['dynamic_viscosity_pa_s','real_near_wall_mb_gap_m','real_near_wall_attenuation_ratio'],
'6_5':['real_old_min_gap_m','real_matched_p5_min_gap_m','real_v1_min_gap_m','real_v1_h_lower_m','real_handoff_event_count','time_accuracy_convergence'],
'7':['inlet_Q_m3_s','outlet_Q_m3_s','mb_number_rate_s_inv','mean_mb_event_interval_s','rbc_volume_rate_m3_s','mass_balance'],
'8':['real_geometry','restart'],
'8_1':['scheduled','admitted','completed','end_reasons','outlet_counts','outlet_frozen_volume_flux_m3_s','outlet_frozen_volume_flux_fraction','completed_path_length_m','completed_residence_time_s','clock','acquisition_birth_window_s','physical_sample_count','statuses','user_acceptance']}
values={s:{'source':stages[s], 'values':{k:records[s][k] for k in ks}} for s,ks in keys.items()}
values['1']['values']['case_lengths_m']=[x['trajectory_length_m'] for x in records['1']['real_trajectory_cases']]
(OUT/'audit/extracted_values.json').write_text(json.dumps(values,ensure_ascii=False,indent=2)+'\n')
selection=[
(2,'particle_3d/reports/particle0/figures/00_frozen_input_overview.png','真实几何与边界'),
(3,'particle_3d/reports/particle0/figures/01_affine_field_validation.png','解析线性场验证算例'),
(4,'particle_3d/reports/particle1/figures/05_real_vessel_single_mb_trajectory.png','真实几何单微泡验证轨迹'),
(5,'particle_3d/reports/particle2/figures/01_rbc_diameter_volume_distribution.png','统计几何模型验证样本'),
(6,'particle_3d/reports/particle2/figures/06_simple_shear_jeffery_orbits.png','解析剪切验证算例'),
(7,'particle_3d/reports/particle3/figures/04_hard_contact_sliding_validation.png','人工壁面验证算例'),
(8,'particle_3d/reports/particle4/figures/02_mb_mb_head_on_contact.png','双球接触验证算例'),
(9,'particle_3d/reports/particle6_5/figures/03_sphere_wall_regularized_resistance.png','近场解析响应验证算例'),
(10,'particle_3d/reports/particle6_5/figures/01_nearfield_activation_weight.png','项目激活函数'),
(11,'particle_3d/reports/particle6_5/figures/07_real_fem_subnanometer_replay.png','真实场同初态局部重放'),
(12,'particle_3d/reports/particle7/figures/02_flux_weighted_inlet_sampling.png','人工与真实入口采样验证'),
(13,'particle_3d/reports/particle8/figures/particle8_07_restart_parity.png','人工生命周期重启验证'),
(14,'particle_3d/outputs/particle8_1/figures/particle8_1_fig_01_full_geometry_trajectories.png','真实几何完整轨迹库'),
(15,'particle_3d/outputs/particle8_1/figures/particle8_1_fig_04_cumulative_localization_views.png','真实轨迹驻留时间累积')]
figs=[]
for n,rel,kind in selection:
 dest=OUT/f'figures/fig_{n:02d}.png';shutil.copyfile(ROOT/rel,dest)
 figs.append(dict(number=n,source=rel,destination=str(dest.relative_to(OUT)),evidence_type=kind,sha256=sha(dest),copy_exact=True))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig,ax=plt.subplots(figsize=(12,6));ax.set(xlim=(0,12),ylim=(0,6));ax.axis('off')
boxes=[(0.2,4.4,'Frozen 3D flow\nGeometry and boundaries'),(4.2,4.4,'Local field sampling\nVelocity and gradients'),(8.2,4.4,'Particle geometry\nMicrobubbles and RBCs'),(8.2,2.5,'Finite size and contact\nWall / pair constraints'),(4.2,2.5,'Near-field resistance\nContinuum handoff'),(.2,2.5,'Continuous inlet flux\nAdmission and waiting'),(.2,.6,'Complete MB journeys\nIndependent trajectories'),(4.2,.6,'Residence accumulation\nULM-style visualization'),(8.2,.6,'Open research questions\nRBC passage / branches')]
for i,(x,y,label) in enumerate(boxes):
 ax.add_patch(FancyBboxPatch((x,y),3.6,1.05,boxstyle='round,pad=.08',facecolor='#eef5f8' if i<8 else '#fff0e6',edgecolor='#346b82' if i<8 else '#ad6a3e',linewidth=1.3));ax.text(x+1.8,y+.525,label,ha='center',va='center',fontsize=11)
for i in range(8):
 x,y,_=boxes[i];xx,yy,_=boxes[i+1]
 if yy==y: start=(x+3.68,y+.525) if xx>x else (x-.08,y+.525);end=(xx-.1,yy+.525) if xx>x else (xx+3.7,yy+.525)
 else:start=(x+1.8,y-.08);end=(xx+1.8,yy+1.16)
 ax.annotate('',xy=end,xytext=start,arrowprops=dict(arrowstyle='->',lw=1.5,color='#406172'))
fig.tight_layout();fig.savefig(OUT/'figures/fig_01.png',dpi=220);fig.savefig(OUT/'figures/fig_01.pdf');plt.close(fig)
figs.append(dict(number=1,source='本次按既有研究逻辑绘制；非数值结果',destination='figures/fig_01.png',evidence_type='方法路线示意',sha256=sha(OUT/'figures/fig_01.png'),copy_exact=False))
d=records['8_1'];labels=list(d['outlet_counts']);q=[d['outlet_frozen_volume_flux_fraction'][x]*100 for x in labels];c=[d['outlet_counts'][x]/d['completed']*100 for x in labels]
fig,ax=plt.subplots(figsize=(10,5));xs=list(range(3));width=.32
ax.bar([x-width/2 for x in xs],q,width,label='Frozen volume flux fraction',color='#3e8195')
ax.bar([x+width/2 for x in xs],c,width,label='Fraction of completed trajectories',color='#b8724f')
for i in xs:
 ax.text(i-width/2,q[i]+2,f'{q[i]:.2f}%',ha='center');ax.text(i+width/2,c[i]+2,f'{c[i]:.0f}%\n(n={d["outlet_counts"][labels[i]]})',ha='center',fontsize=10)
ax.set(xticks=xs,xticklabels=labels,ylabel='Fraction (%)',ylim=(0,122));ax.legend(loc='upper left',frameon=False)
fig.text(.5,.02,'Different denominators: volume flow vs. admitted and numerically completed cohort',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.055,1,1));fig.savefig(OUT/'figures/fig_16.png',dpi=220);fig.savefig(OUT/'figures/fig_16.pdf');plt.close(fig)
figs.append(dict(number=16,source=stages['8_1'],destination='figures/fig_16.png',evidence_type='既有真实结果重绘；无新轨迹计算',sha256=sha(OUT/'figures/fig_16.png'),copy_exact=False))
(OUT/'audit/figure_manifest.json').write_text(json.dumps(sorted(figs,key=lambda x:x['number']),ensure_ascii=False,indent=2)+'\n')
print('Evidence files:',len(inventory),'Main figures:',len(figs))
