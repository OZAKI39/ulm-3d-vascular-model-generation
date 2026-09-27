from pathlib import Path
import json,hashlib,shutil,time,subprocess,os
R=Path(__file__).resolve().parents[1];F=Path('/workspace/hemocell_restore/results/rbc_stage1_20260915_165632');L='/home/lzy/projects/compre_output/rbc_stage1/20260915_165632'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
fit=read(R/'geometry/RBC_GEOMETRY_FIT_AUDIT.json');ind=read(R/'geometry/INDEPENDENT_GEOMETRY_FINAL_AUDIT.json');contract=read(R/'RBC_STAGE1_CONTRACT.json');build=read(R/'RBC_STAGE1_BUILD_PROVENANCE.json');medium=read(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json');wall=read(R/'provenance/WALL_INTERACTION_AUDIT.json');res=read(R/'provenance/RESOURCE_PREFLIGHT.json')
assert fit['RBC_GEOMETRIC_FIT']=='FAIL' and ind['status']=='PASS' and build['status']=='PASS'
protected=read(R/'provenance/PROTECTED_SOURCE_SHA256.json');changes=[p for p,h in protected.items() if sha(Path(p))!=h];assert not changes,changes
orig=read(R/'provenance/ORIGINAL_BUNDLE_MANIFEST.json');native=[x for x in orig['files'] if x['path'].startswith('source/native/')];nativechanges=[x['path'] for x in native if sha(R/x['path'])!=x['sha256']];assert not nativechanges
receipt=read(R/'provenance/NUMERICS_GENERATION_RECEIPT.json');single=sha(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json')==receipt['output_sha256'] and sha(R/'provenance/prepare_numerics.py')==medium['script_sha256']==receipt['script_sha256'] and sha(R/'provenance/NUMERICS_INPUT.json')==medium['input_sha256']==receipt['input_sha256'];assert single
assert not list(R.glob('CASE_*/RUN_STARTED.json')) and not (R/'MPI1_SANITY/RUN_STARTED.json').exists()
active=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  if os.readlink(p/'exe').startswith(str(R/'build')):active.append(p.name)
 except OSError:pass
assert not active
assert not F.exists();F.mkdir(parents=True)
for name in ('source','scripts','contracts','provenance','logs','geometry','verification','MPI1_SANITY','CASE_A','CASE_B','CASE_C'):
 shutil.copytree(R/name,F/name,ignore=shutil.ignore_patterns('__pycache__','*.sock'))
for name in ('RBC_STAGE1_CONTRACT.json','RBC_STAGE1_METRICS.csv','BUILD_STATE.json','RBC_STAGE1_BUILD_PROVENANCE.json'):shutil.copy2(R/name,F/name)
(F/'bin').mkdir();shutil.copy2(R/'build/rbc_geometry_probe',F/'bin/rbc_geometry_probe')
# Preserve original build receipt and append an explicit final scope, without claiming a coupled timestep executable exists.
finalbuild=read(F/'RBC_STAGE1_BUILD_PROVENANCE.json');finalbuild.update(production_coupled_timestep_binary='NOT_BUILT_GEOMETRY_HARD_GATE_FAILED',native_geometry_binary=str(F/'bin/rbc_geometry_probe'),geometry_generation='PASS',MPI1_SANITY='NOT_RUN',HDF5_runtime_particle_output_test='NOT_RUN');save(F/'RBC_STAGE1_BUILD_PROVENANCE.json',finalbuild)
sourceaudit=dict(status='PASS',protected_source_and_medium_files_checked=len(protected),protected_changes=changes,native_restored_files_checked=len(native),native_changes=nativechanges,source_lineage=read(R/'provenance/SOURCE_LINEAGE_AUDIT.json'),prepare_numerics_single_source='PASS',actual_RBC_solver_dt_verification='NOT_RUN',actual_RBC_solver_outlet_density_verification='NOT_RUN',coupled_run_start_records=0,active_task_executables=active,formal_baseline_modified=False)
save(F/'verification/FINAL_SOURCE_INPUT_AUDIT.json',sourceaudit)
summary={
'RBC_ONLY_VALIDATION_STAGE_1':'FAIL','PLATFORM':'VAST RTX4090','SUSPENDING_MEDIUM_NUMERICS':'PASS','SUSPENDING_MEDIUM':'1x PBS + 1% BSA','RHO_KG_M3':1000,'NU_M2_S':1e-6,'TAU':1.,'DT_S':medium['dt_s'],'PREPARE_NUMERICS_SINGLE_SOURCE':'PASS','MOUSE_RBC_MODEL':'MOUSE_RBC_DEV_V0','RBC_TARGET_VOLUME_UM3':50,'RBC_EFFECTIVE_DIAMETER_UM':fit['mesh']['nominal_diameter_um'],'RBC_ACTUAL_DISCRETE_VOLUME_UM3':fit['mesh']['actual_volume_um3'],'RBC_GEOMETRIC_FIT':'FAIL','WALL_INTERACTION_IMPLEMENTATION':'ABSENT','WALL_PENETRATION_METRIC':'VERIFIED','WALL_PENETRATION_METRIC_SCOPE':'STATIC_GEOMETRY_ONLY',
'MPI1_SANITY':'NOT_RUN_GEOMETRY_GATE_FAILED','CASE_A_SINGLE_CENTER':'NOT_RUN_GEOMETRY_GATE_FAILED','CASE_A_STEPS':0,'CASE_B_SINGLE_NEAR_WALL':'BLOCKED_WALL_INTERACTION_IMPLEMENTATION','CASE_B_STEPS':0,'CASE_C_SMALL_PACK':'NOT_RUN_GEOMETRY_GATE_FAILED','CASE_C_RBC_COUNT':'NOT_SELECTED_NO_SAFE_SPAWN','CASE_C_STEPS':0,'MPI4_CORRECTNESS':'UNVERIFIED_NOT_RUN',
'MAX_RELATIVE_VOLUME_ERROR':'UNVERIFIED_NOT_RUN','MAX_RELATIVE_AREA_ERROR':'UNVERIFIED_NOT_RUN','TRIANGLE_INVERSION_COUNT':'UNVERIFIED','WALL_PENETRATION':'UNVERIFIED','NONFINITE_COUNT':'UNVERIFIED_NOT_RUN','MAX_RBC_NODE_VELOCITY':'UNVERIFIED_NOT_RUN','MAX_RBC_NODE_FORCE':'UNVERIFIED_NOT_RUN','RUNTIME_SAFETY':'NOT_RUN','FLUID_SAFETY':'NOT_RUN_RBC_COUPLED_CASES','SOLVER_FINALIZER_IDENTITY':'NOT_RUN','GEOMETRY_FINALIZER_IDENTITY':'PASS','TOTAL_STAGE1_OUTPUT_GIB':'PENDING_SEAL','REMOTE_TO_WSL_INTEGRITY':'NOT_YET_TRANSFERRED',
'HUMAN_RBC_REVIEW':'PENDING','SMALL_PACK_IS_NOT_HCT_VALIDATION':'YES','NEW_MEDIUM_INLET_MULTIPLIER_VALIDATION':'NOT_PERFORMED','MOUSE_SPECIFIC_MEMBRANE_MECHANICS':'UNVERIFIED','MOUSE_SPECIFIC_INTERIOR_VISCOSITY':'UNVERIFIED','EXPERIMENTAL_MOUSE_RBC_CALIBRATION':'PENDING','RBC_GPU_CORRECTNESS':'PENDING','RBC_ONLY_VALIDATION_STAGE_2':'PENDING','MICROBUBBLE_WITH_RBC':'PENDING','READY_FOR_MICROBUBBLE':'NO','PHYSICS_PARAMETERS_TUNED':'NO','FORMAL_PURE_FLUID_BASELINE_MODIFIED':'NO','MICROBUBBLE':'OFF','ADHESION':'OFF','ULTRASOUND':'OFF','RBC_TIMESTEPS':0,'INCOMPLETE_FULL_UPLOAD_EXISTS':'YES' if res['incomplete_full_upload_exists'] else 'NO','INCOMPLETE_FULL_UPLOAD_GIB':res['incomplete_full_upload_GiB'],'INCOMPLETE_FULL_UPLOAD_DELETED':'NO',
'REMOTE_RESULT_DIR':str(F),'LOCAL_RESULT_DIR':L,'RBC_STAGE1_REPORT':L+'/RBC_STAGE1_REPORT.md','NEXT_STEP':'USER_REVIEW_BEFORE_RBC_STAGE2'}
final=dict(status='FAIL',reason='RBC_GEOMETRIC_FIT_FAILED',additional_blockers=['FIXED_NOMINAL_SCALE_NATIVE_MESH_VOLUME_DIFFERS_FROM_50_UM3','WALL_INTERACTION_ABSENT_FOR_PRODUCTION_GUO_SETUP'],geometry_audit=fit,independent_geometry_audit=ind,source_and_input_audit=sourceaudit,case_statuses=contract['cases'],solver_finalizer_identity='NOT_RUN',runtime_claims='No RBC creation/IBM/MPI1-sanity/MPI4 correctness claim; only native library compilation and static reference mesh generation performed',RBC_timesteps=0)
save(F/'RBC_STAGE1_FINAL_AUDIT.json',final)
m=fit['mesh'];b=fit['best_pose'];dxum=medium['dx_m']*1e6
geom=f'''# RBC 几何适配硬门报告

RBC_GEOMETRIC_FIT = FAIL

本轮使用固定的原生 HemoCell RBC_FROM_SPHERE 网格。名义缩放 cbrt(50/90)={m['nominal_scale']:.17g}，名义直径 {m['nominal_diameter_um']:.9f} µm；实际顶点最大跨度 {m['maximum_vertex_extent_um']:.9f} µm。网格为 {m['vertices']} 顶点、{m['triangles']} 三角形，三个轴向跨度为 {m['axis_extents_um']} µm，最薄轴跨度约1.886 µm。未修改形状、网格分辨率或机械常数。

## 固定缩放的实际体积

原生离散网格实际体积={m['actual_volume_um3']:.12f} µm³，面积={m['actual_area_um2']:.12f} µm²。与名义目标50 µm³的相对差异={ind['relative_volume_difference_from_nominal_target']:.9%}。50/90为用户给出的名义缩放关系，原生离散网格不是精确90 µm³的参考体积，因此不能把固定缩放后的几何声称为精确50 µm³。未为了消除该差异重新设定scale、半径或网格。本结果不是运行中的体积漂移。

## 真实网格放置搜索

闭合 STL 与 Step2 closed voxel mask 保持冻结。候选中心取既有lumen中每2 LU网格点，共{fit['coarse_centers']}个；131个确定性法向、每个0/30/60度旋转，共{fit['coarse_poses']}个粗筛姿态。粗筛用真实网格的极值顶点与固定子集，再对候选检查全部642顶点；对48个候选按预先固定的8轮平移/旋转规则细化。它只调整预运行放置姿态，不调整细胞形状、尺寸或物理参数。

CASE A/C要求所有顶点有至少2dx的初始壁面间隙，并且RBC曲面与血管曲面没有相交。2dx是IBM核支撑对应的放置安全余量，不是允许穿墙的容差。穿墙容差为几何零：signed distance>0即越出闭合lumen。

没有找到可接受放置。最佳诊断候选仍有{b['outside_vertices']}/642个顶点在lumen外，最小间隙={b['min_distance']:.12f} LU={b['min_distance']*dxum:.12f} µm，检测到{b['wall_triangle_contacts']}个三角形接触对。这个候选仅供失败诊断，绝不是已冻结可运行的spawn。

此结论来自有界确定性搜索，不是对所有连续姿态的数学不可嵌入证明。当前没有获得可审核的安全初态，因此按本任务硬门停止；不能通过放宽间隙、缩小RBC或事后移动初始位置来宣称PASS。

## 独立壁面判据复核

以原始闭合STL的负signed distance为正向内壁间隙；闭合面包含端口cap，使用它对离开封闭lumen作保守筛查。保留native 0.001 LU inflate，不在诊断中扩大原始lumen。用既有inside格点及已知outside点验证sign；确认STL闭合且无非流形边。

另一个独立脚本以全曲面solid-angle winding number判定inside/outside，以独立nearest-triangle locator重算距离，以VTK mass properties重算体积/面积。与搜索代码一致，GEOMETRY_FINALIZER_IDENTITY=PASS。这只验证静态几何证据，不是SOLVER_FINALIZER_IDENTITY=PASS。

## ParaView 文件

geometry/FROZEN_LUMEN_LU.vtp 与 geometry/BEST_RIGID_PLACEMENT.vtp 使用同一LU坐标系，可共同打开查看。后者明确为失败诊断姿态；转换到物理坐标需使用contracts中的origin与dx。没有成功case，未伪造initial/mid/final运行帧。HUMAN_RBC_REVIEW=PENDING。
'''
(F/'RBC_GEOMETRY_FIT_REPORT.md').write_text(geom)
report=f'''# RBC_ONLY_VALIDATION_STAGE_1：当前生产版本

**RBC_ONLY_VALIDATION_STAGE_1 = FAIL。** 停在真实RBC网格几何适配硬门。MPI1 sanity及CASE A/B/C全部未运行，RBC_TIMESTEPS=0。另有当前Guo血管配置缺少RBC-wall作用实现，以及固定名义缩放的实际离散体积偏离50 µm³两项需审阅的问题。

## 本轮实际完成的工作

在Vast RTX4090实例建立独立工作目录{R}。从只读minimal restore恢复并SHA256验证1512个原生文件；完整构建HemoCell CPU/MPI静态库和新的原生RBC网格生成程序。没有使用旧服务器未知binary，没有修改HemoCell core或纯流体基线。

HemoCell commit=5a410848bd5c57d5ae1c171112e78eab4a82e650；Palabos commit=05712164d940a42e06afdd705249912fa0c49f14。额外从已校验Palabos原始源码包应用既有HemoCell patch重建，与1421个Palabos文件逐一一致。编译器、MPI、CMake、并行HDF5、flags与SHA256见RBC_STAGE1_BUILD_PROVENANCE.json。并行HDF5配置与链接通过；RBC HDF5运行输出管线未进行测试。

构建范围必须区分：原生CPU库和geometry construction executable为PASS；没有发布或声称完成可执行IBM timestep的生产case binary。几何硬门失败后，不继续建立依赖有效spawn的耦合case入口，也没有MPI1 sanity或MPI4运行证据。原生网格生成在服务器完成（0 timestep），几何搜索和独立复核在WSL已有NumPy/VTK环境完成。

## 介质合同保持

先前PBS/BSA pure-fluid smoke与SUSPENDING_MEDIUM_NUMERICS=PASS保留；ρ=1000 kg/m³、ν=1e-6 m²/s、tau=1、dt={medium['dt_s']:.17g} s。三个出口rho_LU为{medium['outlet_01_rho_lu']:.17g}、{medium['outlet_02_rho_lu']:.17g}、{medium['outlet_03_rho_lu']:.17g}。直接复用已经验证的生成合同，prepare_numerics仍是唯一生成源，本轮不再生成另一套dt或出口密度。副本与生成脚本、输入、回执的SHA256相符。

本轮没有RBC耦合solver，因此新的RBC solver实际读取dt/出口密度一致性是NOT_RUN，不能借用纯流体PASS来冒充RBC运行PASS。介质仍为1×PBS+1%BSA、25°C的开发假设，EXPERIMENTALLY_MEASURED=NO，最终实验合同PENDING。入口倍率1.1197286861799598未重新校准，新介质校准仍NOT_PERFORMED。

## 几何失败证据

细节见RBC_GEOMETRY_FIT_REPORT.md。固定scale={m['nominal_scale']:.17g}，名义直径{m['nominal_diameter_um']:.9f}µm。原生网格642顶点/1280三角形，实际体积{m['actual_volume_um3']:.12f}µm³（比名义50低约9.907%），面积{m['actual_area_um2']:.12f}µm²。该差异是初始化几何差异，不是运行volume error。

{fit['coarse_poses']}个确定性粗筛姿态及预定连续细化未找到安全初态。最佳失败姿态有154个顶点在lumen外、最深越界约0.597581µm，并有186个曲面接触对。未自动缩小RBC，也未调整网格、kLink/kArea/kBend/kVolume、粘度或tau。模型常数15/5/80/20仅被冻结为本轮开发定义；膜力、IBM、体积与面积时间演化均未验证。

## CASE B 的独立阻塞

HemoCell库确实存在enableBoundaryParticles、populateBoundaryParticles、applyBoundaryRepulsionForce。现有血管实现没有启用该路径。其壁面粒子生成依赖Dynamics::isBoundary()，而当前血管是BGK/NoDynamics配合Guo off-lattice processors，不提供该格点壁面粒子集合。流体Guo no-slip不能视作RBC排斥力。

所以WALL_INTERACTION_IMPLEMENTATION=ABSENT（限定当前production vascular setup），CASE_B=BLOCKED_WALL_INTERACTION_IMPLEMENTATION。未臆造force/cutoff，未把Guo换成bounce-back，也未修改HemoCell核心来绕过此条件。库中接口存在不等于当前血管应用具备已验证的壁面相互作用。

## 阶段与未验证项

|阶段|结果|实际步数|
|---|---|---:|
|Native CPU/MPI library + reference mesh build|PASS|0|
|真实RBC几何fit|FAIL|0|
|MPI1 sanity|NOT_RUN_GEOMETRY_GATE_FAILED|0|
|CASE A / MPI4|NOT_RUN_GEOMETRY_GATE_FAILED|0|
|CASE B / MPI4|BLOCKED_WALL_INTERACTION_IMPLEMENTATION，同时几何门未通过|0|
|CASE C / MPI4|NOT_RUN_GEOMETRY_GATE_FAILED；未选定2–4个安全位置|0|

没有以MPI1回退替代MPI4。没有RBC count/ownership、IBM耦合、nonfinite时间序列、force/velocity runaway、面积/体积漂移、动态triangle inversion或运行中穿墙的PASS证据。对应数值为UNVERIFIED/NOT_RUN，不能填0冒充验证。wall metric仅静态验证；runtime wall penetration为UNVERIFIED。

RBC_STAGE1_CONTRACT.json记录用户指定阈值、计划步数/核数、监测与输出政策，以及因无安全spawn而故意保留的null位置。任何未来运行都必须先解决几何定义/适配与壁面实现问题，重新冻结可审阅的有效初态和运行合同；本轮失败合同禁止启动。

## 资源与归档

资源预检剩余磁盘约{res['free_disk_bytes']/2**30:.3f}GiB，超过20GiB门槛。旧incomplete full upload仍存在，实际占用{res['incomplete_full_upload_GiB']:.9f}GiB（有重复硬链接，逻辑文件大小合计约{res['incomplete_full_upload_logical_GiB']:.6f}GiB），未删除。结果只包含代码/构建证据、冻结合同、紧凑CSV、网格及失败姿态、日志和核查报告；不含整个临时build树，无周期checkpoint或大量ParaView帧。输出大小见OUTPUT_SIZE.json，小于1GiB。

报告和SHA256结果位于{F}，随后下载到{L}。远端原始报告/manifest保持不可变；WSL下载SHA256与独立静态几何复核记录见LOCAL_TRANSFER_VERIFICATION.json、LOCAL_GEOMETRY_AUDIT.json与LOCAL_SHA256SUMS。

## 保留的限制与后续

HUMAN_RBC_REVIEW=PENDING。SMALL_PACK_IS_NOT_HCT_VALIDATION=YES。mouse membrane mechanics与interior viscosity为UNVERIFIED；实验mouse RBC标定、RBC_GPU_CORRECTNESS、RBC Stage2、MICROBUBBLE_WITH_RBC均PENDING；READY_FOR_MICROBUBBLE=NO。microbubble、adhesion、ultrasound全部OFF。

NEXT_STEP=USER_REVIEW_BEFORE_RBC_STAGE2；这里首先需要审阅并解决Stage1失败，绝不表示已允许进入Stage2。本任务没有自动运行任何Hct群体、GPU IBM、GPU particle mechanics或微泡任务。
'''
(F/'RBC_STAGE1_REPORT.md').write_text(report)
(F/'geometry/README.md').write_text('此目录的BEST_RIGID_PLACEMENT.vtp为明确不合格的失败诊断姿态，不是可运行spawn。与FROZEN_LUMEN_LU.vtp共同打开；单位LU。没有RBC运行帧。\n')
# Generate final summaries then measure their actual final file footprint.
summary['TOTAL_STAGE1_OUTPUT_GIB']=sum(p.stat().st_size for p in F.rglob('*') if p.is_file())/2**30
save(F/'FINAL_SUMMARY.json',summary);(F/'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(f'{k} = {v}' for k,v in summary.items())+'\n')
size=sum(p.stat().st_size for p in F.rglob('*') if p.is_file());assert size<2**30
save(F/'OUTPUT_SIZE.json',dict(status='PASS',bytes_before_manifest=size,GiB_before_manifest=size/2**30,limit_bytes=2**30,scope='All compact result regular files including native sources and new geometry executable; no external links or full build tree'))
lines=[sha(p)+'  '+str(p.relative_to(F)) for p in sorted(F.rglob('*')) if p.is_file()];(F/'SHA256SUMS').write_text('\n'.join(lines)+'\n')
size=sum(p.stat().st_size for p in F.rglob('*') if p.is_file());assert size<2**30
save(R/'FINALIZATION_TERMINAL.json',dict(status='PASS',stage_status='FAIL',reason='GEOMETRY_HARD_GATE_FAILED',result_dir=str(F),manifest_sha256=sha(F/'SHA256SUMS'),file_count=len(lines),output_bytes=size,output_GiB=size/2**30,solver_steps=0,unix=time.time()))
print(json.dumps(summary,indent=2))
