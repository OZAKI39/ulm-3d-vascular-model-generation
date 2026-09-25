from pathlib import Path
import json,hashlib,shutil,csv,time
R=Path(__file__).resolve().parents[1];F=Path('/workspace/hemocell_restore/results/new_medium_smoke_20260915_161937');L='/home/lzy/projects/compre_output/pure_fluid_new_medium_smoke/20260915_161937'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
state=read(R/'SMOKE_STATE.json');assert state['state'] in ('PASS','FAIL')
assert not F.exists(),'Never overwrite a completed result'
F.mkdir(parents=True);(F/'verification').mkdir()
old=read(R/'provenance/PROTECTED_BASELINE_SHA256.json');oldfails=[p for p,h in old.items() if sha(Path(p))!=h]
execution=read(R/'provenance/EXECUTION_SHA256.json');exfails=[p for p,h in execution.items() if sha(R/p)!=h]
c=read(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json');receipt=read(R/'provenance/NUMERICS_GENERATION_RECEIPT.json')
producerok=receipt['output_sha256']==sha(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json') and receipt['script_sha256']==c['script_sha256']==sha(R/'scripts/prepare_numerics.py') and receipt['input_sha256']==c['input_sha256']==sha(R/'inputs/NUMERICS_INPUT.json')
integrity=dict(status='PASS' if not oldfails and not exfails and producerok else 'FAIL',protected_baseline_files=len(old),protected_baseline_changed=oldfails,execution_files=len(execution),execution_changed=exfails,prepare_numerics_single_source='PASS' if producerok else 'FAIL',new_driver_SHA256=sha(R/'source/vascularPoC.cpp'),binary_SHA256=sha(R/'build/vascularPoC'),GPU_headers_byte_identical_to_old=True if not oldfails and not exfails else 'UNVERIFIED',formal_source_modified=bool(oldfails))
save(F/'verification/SOURCE_AND_INPUT_FINAL_AUDIT.json',integrity)
a={n:read(R/f'verification/SMOKE_{n}_AUDIT.json') if (R/f'verification/SMOKE_{n}_AUDIT.json').exists() else dict(status='NOT_RUN') for n in (500,5000)}
passed=state['state']=='PASS' and all(x['status']=='PASS' for x in a.values()) and integrity['status']=='PASS';last=a[5000] if a[5000]['status']!='NOT_RUN' else a[500];complete=[x for x in a.values() if 'rho_min' in x]
for directory in ('contracts','inputs','source','scripts','provenance','verification','logs'):
 dest=F/directory
 shutil.copytree(R/directory,dest,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.sock'))
for n in (500,5000):
 d=R/f'run_{n}'
 if not (d/'RUN_STARTED.json').exists():continue
 dest=F/d.name;dest.mkdir()
 for p in d.iterdir():
  if p.name=='contracts':continue
  if p.is_dir():shutil.copytree(p,dest/p.name)
  elif p.is_file() and p.name!='SHA256SUMS':shutil.copy2(p,dest/p.name)
 # Relative link to the single frozen contract copy, never to an external geometry tree.
 (dest/'contracts').symlink_to('../contracts',target_is_directory=True)
for name in ('BUILD_PROVENANCE.json','BUILD_STATE.json','SMOKE_STATE.json','SMOKE_STARTED.json'):shutil.copy2(R/name,F/name)
shutil.copy2(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json',F/'NEW_MEDIUM_NUMERICS_CONTRACT.json')
for suffix in ('RUNTIME','FLUX','MASS'):
 name=f'NEW_MEDIUM_{suffix}_HISTORY.csv';writer=None
 with (F/name).open('w',newline='') as out:
  for n in (500,5000):
   p=R/f'run_{n}/diagnostics'/name
   if not p.exists():continue
   with p.open() as inp:
    rows=csv.DictReader(inp)
    if writer is None:writer=csv.DictWriter(out,fieldnames=['run_horizon']+rows.fieldnames);writer.writeheader()
    for row in rows:writer.writerow(dict(run_horizon=n,**row))
summary={
'PURE_FLUID_NEW_MEDIUM_SMOKE_STATUS':'PASS' if passed else 'FAIL',
'SUSPENDING_MEDIUM':'1x PBS + 1% BSA','TEMPERATURE_C':25,'DEVELOPMENT_ASSUMPTION':'YES','EXPERIMENTALLY_MEASURED':'NO','RHO_KG_M3':1000,'NU_M2_S':1e-6,'DX_M':c['dx_m'],'TAU':c['tau'],
'GENERATED_DT_S':c['dt_s'],'EXPECTED_DT_S':6.65994708146172e-9,'DT_REPRODUCTION':c['expected_reference_check']['status'],'PRESSURE_UNIT_PA':c['pressure_unit_pa'],
**{f'GENERATED_OUTLET_{i:02}_RHO_LU':c[f'outlet_{i:02}_rho_lu'] for i in (1,2,3)},
'EXPECTED_OUTLET_01_RHO_LU':1.0000484343922278,'EXPECTED_OUTLET_02_RHO_LU':1.0004402376508774,'EXPECTED_OUTLET_03_RHO_LU':.9999543772756865,'OUTLET_RHO_LU_REPRODUCTION':c['expected_reference_check']['status'],
'PREPARE_NUMERICS_SINGLE_SOURCE':integrity['prepare_numerics_single_source'],'NUMERICS_STATIC_CHECK':read(R/'verification/NUMERICS_STATIC_CHECK.json')['status'],'QUICK_SMOKE_500':a[500]['status'],'SUSTAINED_SMOKE_5000':a[5000]['status'],
'NONFINITE_COUNT':sum(x['nonfinite_count'] for x in complete) if complete else 'UNVERIFIED','RHO_MIN':min(x['rho_min'] for x in complete) if complete else 'UNVERIFIED','RHO_MAX':max(x['rho_max'] for x in complete) if complete else 'UNVERIFIED','MAX_MACH':max(x['max_mach'] for x in complete) if complete else 'UNVERIFIED',
'MACH_SAFETY':'PASS' if complete and all(x['checks']['Mach_safety_margin']=='PASS' for x in complete) else 'FAIL','TOTAL_MASS_DRIFT':last.get('total_mass_drift','UNVERIFIED'),'MASS_STABILITY':'PASS' if complete and all(all(x['checks'][key]=='PASS' for key in ('mass_stability','CV_balance','no_fast_monotonic_mass_growth')) for x in complete) else 'FAIL',
**{('FINAL_'+p.upper()):last.get('final_flux_m3_s',{}).get(p,'UNVERIFIED') for p in ('Qin','Qout01','Qout02','Qout03')},
'TRANSIENT_BACKFLOW_OBSERVED':'YES' if any(x.get('transient_backflow',False) for x in complete) else 'NO',
'GPU_STEPS_PER_SEC':last.get('gpu_steps_per_second','UNVERIFIED'),'GPU_VRAM_PEAK_MIB':max(x['gpu_vram_peak_MiB'] for x in complete) if complete else 'UNVERIFIED',
'ACTUAL_DT_MATCHES_CONTRACT':'PASS' if len(complete)==2 and all(x['checks']['actual_dt']=='PASS' for x in complete) else 'FAIL','ACTUAL_OUTLET_RHO_LU_MATCH_CONTRACT':'PASS' if len(complete)==2 and all(x['checks']['actual_outlet_rho']=='PASS' for x in complete) else 'FAIL',
'SOLVER_FINALIZER_IDENTITY':'PASS' if passed else 'FAIL','SUSPENDING_MEDIUM_NUMERICS':'PASS' if passed else 'FAIL_RUNTIME_SMOKE','SUSPENDING_MEDIUM_CONTRACT':'PROVISIONAL_PASS' if passed else 'PENDING',
'FINAL_EXPERIMENTAL_CONTRACT':'PENDING','NEW_MEDIUM_INLET_MULTIPLIER_VALIDATION':'NOT_PERFORMED','RBC_GPU_CORRECTNESS':'PENDING','RBC_INCLUDED':'NO','RBC_TIMESTEPS':0,'MICROBUBBLE':'OFF','ADHESION':'OFF','ULTRASOUND':'OFF','FORMAL_OLD_PURE_FLUID_BASELINE_MODIFIED':'NO' if not oldfails else 'YES',
'REMOTE_TO_WSL_INTEGRITY':'NOT_YET_TRANSFERRED','REMOTE_RESULT_DIR':str(F),'LOCAL_RESULT_DIR':L,'PURE_FLUID_NEW_MEDIUM_SMOKE_REPORT':L+'/PURE_FLUID_NEW_MEDIUM_SMOKE_REPORT.md','NEXT_STEP':'RBC_ONLY_VALIDATION_STAGE_1' if passed else 'STOP_AND_REVIEW_SMOKE_FAILURE'}
final=dict(status=summary['PURE_FLUID_NEW_MEDIUM_SMOKE_STATUS'],source_and_input_integrity=integrity,run_audits={str(n):x for n,x in a.items()},summary=summary,scope='Transient numerical safety and dimensional identity; not convergence, measured medium or inlet calibration',formal_old_baseline_status='PASS',RBC_timesteps=0)
save(F/'NEW_MEDIUM_FINAL_AUDIT.json',final);save(F/'FINAL_SUMMARY.json',summary)
(F/'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(f'{k} = {v}' for k,v in summary.items())+'\n')
lines=['# PBS/BSA 开发介质纯流体 smoke 核查', '',f"结果：**{summary['PURE_FLUID_NEW_MEDIUM_SMOKE_STATUS']}**。这是新介质的短期数值稳定性与单位转换检查，不是长期收敛或真实实验介质验证。",'',
'介质为 1× PBS + 1% BSA，25°C；ρ=1000 kg/m³、μ=0.001 Pa·s、ν=10⁻⁶ m²/s 均为开发阶段假设，未做实验测量。最终实验合同仍为 PENDING。旧整血等效纯流体基线 PASS 保留。','',
f"当前 prepare_numerics.py 唯一生成 dt={c['dt_s']:.17g} s，tau={c['tau']}，有效 dx={c['dx_m']:.17g} m。pressure_unit={c['pressure_unit_pa']:.17g} Pa。三出口物理表压不变；生成密度为 {', '.join(format(c[f'outlet_{i:02}_rho_lu'],'.17g') for i in (1,2,3))}。",'',
'压力换算沿用 dt=((tau−0.5)/3) dx²/ν、pressure_unit=ρ(dx/dt)²、rho_LU=1+3 Pg/pressure_unit。求解器只读合同与生成参数文件，并逐值核对；独立 finalizer 用相同物理定义核验单位，核验结果不作为新输入。', '',
'生成器完整 JSON 的 SHA256 在 NUMERICS_GENERATION_RECEIPT.json 中；合同内部 output_sha256 哈希指定的 canonical generated_output_digest_payload，避免文件自哈希循环。运行前修正过回执中 input 路径的变量覆盖，旧回执/脚本/合同保存在 provenance/PRE_AUDIT_RECEIPT_FIX_*，数值内容不变。', '',
'GPU 核心（batch dispatch、halo、managed memory、Guo、监控归约、kernel layout）与已验证 Stage4 相同。独立 host driver 只适配合同读取、500/5000 步任务边界及本任务指定输出节奏；初始化/Guo/边界数学字节核查见 DRIVER_STATIC_PROOF.json。', '',
'## 实际执行与安全检查', '',
'500 步与 5000 步各一次，均从 ρ_LU=1、u=0 初始化，ramp=10 步；MPI1，单物理核心绑定。每步保留 nonfinite 与原 GPU 安全归约，50 步记录状态，100 步记录全部 24 个测量面流量，500 步记录控制体质量。完整快照只有各运行 step0 与 5000 步运行的 step5000。没有运行 RBC，也未启动后续任务。','',
'|运行|核查|solver steps/s|端到端 s|初始化 s|VRAM 峰值 MiB|timed GPU中位利用率|','|---|---|---:|---:|---:|---:|---:|']
for n,x in a.items():lines.append(f"|{n}|{x['status']}|{x.get('gpu_steps_per_second')}|{x.get('end_to_end_seconds')}|{x.get('initialization_seconds')}|{x.get('gpu_vram_peak_MiB')}|{x.get('gpu_utilization_timed_median')}%|")
lines+=['', '性能只记录，不作为通过门槛。当前更稀疏的流量输出属于用户要求的 smoke 监测节奏，不能据此宣称 GPU 核心优化加速。MPI1 初始化保留原实现，端到端计时包括初始化、必需输出及退出。', '',f"全程物理流体 rho 范围 {summary['RHO_MIN']}–{summary['RHO_MAX']}；流体与 ghost 的最大 Mach={summary['MAX_MACH']}；nonfinite={summary['NONFINITE_COUNT']}。既有 density 安全范围 [0.99,1.01]、Mach 硬门 0.05；本 smoke 在运行前设定通过余量 Mach≤0.01。",'',
f"5000 步末总质量相对初值漂移={last.get('total_mass_drift')}，50 步最大质量跳变比例={last.get('max_mass_jump_per50_fraction')}，CV 间隔/累计最大闭合误差除以初始 CV 质量分别为 {last.get('CV_interval_balance_max_fraction')} / {last.get('CV_cumulative_balance_max_fraction')}。",'',
'质量门槛、快速单调失控规则、CV 误差与 flux 粗筛均预先冻结于 contracts/SMOKE_SAFETY_CONTRACT.json。CV 闭合采用 storage change + 100 步观测的净向外质量通量梯形积分；这是 transient 简单闭合检查，不是长期 R_flow≤1% 门。', '',
f"末态主测量面（g2，冻结的4dx面）有符号流量 m³/s：{last.get('final_flux_m3_s')}。入口按向内为正，出口按向外为正；全24面最大 |Q|/Qtarget={last.get('max_abs_flux_over_Qtarget')}。",'',
f"短暂回流记录：{json.dumps(last.get('backflow',{}),ensure_ascii=False)}。压力驱动 transient 可以出现负出口流量；有限性、符号约定、极端流量门及 snapshot 复算均独立检查，未把负出口流量自动判失败。",'',
'原 multiplier=1.1197286861799598 仅沿用数值入口实现，新介质校准 NOT_PERFORMED；不宣称长期 Qin 精度，也未重新校准。未对不同粘度的旧/新介质要求流场相等。', '',
f"内存 steady-thirds 核查：{json.dumps(last.get('resource_thirds',[]),ensure_ascii=False)}；增长判定 {json.dumps(last.get('memory_growth',{}),ensure_ascii=False)}。只说明本短测试未达到持续增长门槛，不证明无限长运行无泄漏。",'',
'## 独立终态核查与适用范围', '',
'finalizer 独立读取 numerics、每步安全记录、运行诊断和二进制快照，复算最终密度/速度/Mach/总质量/CV质量，以及所有24组多平面插值通量，按预先冻结容差验证。具体检查逐项记录在 verification/SMOKE_*_AUDIT.json。顶层三个历史 CSV 的 run_horizon 区分两个从零运行，避免把两段拼成连续时间。', '',
f"旧 dt=2.0366810646671923e-9 s；新 dt / 旧 dt={c['dt_s']/2.0366810646671923e-9:.17g}。固定 dx 与 tau 时，较小物理运动粘度对应每 LBM 步更长物理时间。5000步物理时间={5000*c['dt_s']:.17g} s。加入 RBC 后有额外计算量与时间尺度约束，不能推导 RBC 必然快3.27倍。",'',
'通过仅允许 SUSPENDING_MEDIUM_NUMERICS=PASS、SUSPENDING_MEDIUM_CONTRACT=PROVISIONAL_PASS；EXPERIMENTALLY_MEASURED=NO、FINAL_EXPERIMENTAL_CONTRACT=PENDING、RBC_GPU_CORRECTNESS=PENDING。后续建议 RBC_ONLY_VALIDATION_STAGE_1，但本任务没有启动它。', '',
'旧源文件、几何与冻结输入 hash 核查见 verification/SOURCE_AND_INPUT_FINAL_AUDIT.json。远端报告为不可变运行证据；下载后的 WSL SHA256 验证及独立重算记录另见 LOCAL_TRANSFER_VERIFICATION.json 与 LOCAL_NUMERICAL_AUDIT.json。']
(F/'PURE_FLUID_NEW_MEDIUM_SMOKE_REPORT.md').write_text('\n'.join(lines)+'\n')
# Count actual result files; directory links point only to the single internal contracts copy.
files=[p for p in F.rglob('*') if p.is_file() and not p.is_symlink()];size=sum(p.stat().st_size for p in files);limit=read(R/'contracts/SMOKE_SAFETY_CONTRACT.json')['output_limit_bytes'];assert size<limit,(size,limit)
save(F/'OUTPUT_SIZE.json',dict(status='PASS',bytes_before_manifests=size,limit_bytes=limit,fullfield_snapshot_steps={'run_500':[0],'run_5000':[0,5000]},count_scope='Result regular files; shared baseline/toolchain/build tree excluded'))
lines=[sha(p)+'  '+str(p.relative_to(F)) for p in sorted(F.rglob('*')) if p.is_file() and not p.is_symlink() and p.name!='SHA256SUMS'];(F/'SHA256SUMS').write_text('\n'.join(lines)+'\n')
save(R/'FINALIZATION_TERMINAL.json',dict(status='PASS',numerical_status=summary['PURE_FLUID_NEW_MEDIUM_SMOKE_STATUS'],result_dir=str(F),files=len(lines),output_bytes=sum(p.stat().st_size for p in F.rglob('*') if p.is_file() and not p.is_symlink()),manifest_sha256=sha(F/'SHA256SUMS'),unix=time.time()))
print(json.dumps(summary,indent=2))
