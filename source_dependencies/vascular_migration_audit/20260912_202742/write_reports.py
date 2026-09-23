"""Build review-only reports from the observed local code, configuration and data."""
from pathlib import Path
import json, ast, collections, re, hashlib, csv

S=Path('/home/lzy/projects/ulm_3D_vascular')
T=Path('/home/lzy/projects/hemocell_starter')
A=Path(__file__).resolve().parent
def load(n):return json.loads((A/n).read_text())
def dump(n,x):(A/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def ev(path,needle=None,kind='CONFIRMED_FROM_CODE'):
    p=S/path;lines=p.read_text().splitlines();i=next((i for i,l in enumerate(lines,1) if needle and needle in l),1)
    return {'file':str(p),'line':i,'context':lines[i-1][:1000],'classification':kind}
def link(e):return f"[{Path(e['file']).relative_to(S)}:{e['line']}](<{e['file']}:{e['line']}>)"
def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n|'+'|'.join('---' for _ in headers)+'|\n'+''.join('| '+' | '.join(str(c).replace('|','/').replace('\n','; ') for c in row)+' |\n' for row in rows)

idx=load('code_structure_index.json');git=load('source_git_worktree_assessment.json')
stats=load('source_size_and_format_summary.json');headers=load('representative_data_headers.json')
surface='outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611'
model='outputs/model_generate/ultraliser_anchor003274_20260907_192331'
pre='outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628'
base='outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901'
prom='outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637'

# Each source module has a static import/function basis; deeper API compatibility remains untested.
inventory=[]
special={
 'utils/mesh/quality.py':('independent','ADAPT','只读拓扑/几何度量；完整结果结构硬编码 _um，需显式单位接口'),
 'utils/cfd_flow/geometry.py':('mixed','ADAPT','CellEntityIds 分区、边界法向与米制 STL；需解除固定 1 入 3 出和 FlowInputs 耦合'),
 'utils/cfd_flow/dimensionless_geometry_kernel.py':('mixed','ADAPT','无量纲射线三角形核可提取；同模块仍含旧 q 值/阶段账本逻辑'),
 'utils/cfd_flow/physical_port_flux.py':('mixed','ADAPT','连续截面与积分有价值；场输入、平面契约和旧求解器解析需要适配'),
 'utils/cfd_flow/steady_export.py':('mixed','ADAPT','规则中心点场转六面体 VTK；保留物理坐标，解除旧场命名约定'),
 'utils/cfd_flow/steady_state.py':('mixed','ADAPT','物理时间窗统计可复用；当前直接导入固定 Tau1 端口/阈值/流量契约'),
 'utils/cfd_preprocess/one_d_flow.py':('independent','ADAPT','SciPy 稀疏一维阻力网络；物理假设和边界参数需要单独批准'),
 'utils/cfd_lumen/ultraliser_backend.py':('independent','ADAPT','血管 H5 适配及外部几何工具执行；不是血流 solver，但有路径/程序依赖'),
 'utils/cfd_surface_prepare/vmtk_runner.py':('independent','ADAPT','VMTK 外部几何工具请求；当前 Windows 环境路径在 WSL 下不能照搬'),
 'utils/sampling/feature_scaling.py':('independent','KEEP_AS_IS','纯特征缩放；连同独立 ScalerState 数据结构复用'),
 'utils/rodent_vasculature/geometry.py':('independent','KEEP_AS_IS','NumPy 中心线弧长/重采样数学逻辑；调用者必须遵守已声明单位'),
 'utils/cfd_flow/validated_contract.py':('specific','REFERENCE_ONLY','旧 Tau1 数值尺度与硬编码参数；禁止直接映射 HemoCell'),
 'utils/cfd_flow/restart_decode.py':('specific','REFERENCE_ONLY','旧 APES/Lua/LSB 重启、PDF 分布及宏观场解析'),
 'doc_visualize_v1.py':('independent','ADAPT','读取保存的 SWC/ROI 结果并生成文档图；路径和输出选择需适配'),
 'doc_visualize_v2.py':('independent','ADAPT','基于几何生成过程与已有数据的展示；专用 UI 和路径需适配'),
 'doc_visualize_v3.py':('mixed','ADAPT','旧稳态 VTU/物理截面/壁面力字段可视化；需新的 HemoCell 字段契约'),
}
for row in idx:
    rel=row['path'];p=S/rel;text=p.read_text();tree=ast.parse(text)
    if rel in special:dependency,cls,reason=special[rel]
    elif rel.startswith('utils/cfd_flow/') or rel=='s5_cfd_flow_solve.py':
        dependency='specific';cls='REWRITE' if rel.endswith(('pipeline.py','production.py','apes.py','s5_cfd_flow_solve.py','config.py')) else 'REFERENCE_ONLY';reason='当前 imports/API 指向 APES/Musubi、旧数值契约或其结果；仅功能需求可作为未来接口参考'
    else:dependency='independent';cls='ADAPT';reason='与旧血流 solver 无直接依赖；保留算法前需适配项目数据结构、路径、单位和副作用接口'
    funcs=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and not n.name.startswith('_')]
    evidence=[ev(rel)]
    for f in row['imports'][:4]:evidence.append({'file':str(p),'line':f['line'],'context':f['text'],'classification':'CONFIRMED_FROM_CODE'})
    for f in row['functions'][:2]:evidence.append({'file':str(p),'line':f['line'],'context':f['name']+': calls '+', '.join(f['calls'][:12]),'classification':'CONFIRMED_FROM_CODE'})
    imports=[x['text'] for x in row['imports']]
    record={'source_path':str(p),'component':rel,'purpose':ast.get_docstring(tree) or 'Module exports/imports; see evidence',
        'entrypoint':bool(row['main_guard_lines']),'solver_dependency':dependency,
        'input':[f.name+'('+ast.unparse(f.args)+')' for f in funcs[:6]],
        'output':[f.name+' -> '+(ast.unparse(f.returns) if f.returns else 'UNVERIFIED') for f in funcs[:6]],
        'units':'UNIT_UNVERIFIED at module boundary; see units_inventory.json for evidenced contracts',
        'coordinate_system':'COORDINATE_SYSTEM_UNVERIFIED unless a specific transform is cited',
        'dependencies':imports,'experimental_relevance':'No measurement provenance established by module name; see experimental_assets.json',
        'migration_class':cls,'migration_reason':reason,'migration_decision_classification':'INFERRED',
        'modify_required':cls in ('ADAPT','REWRITE'),'target_candidate':'候选：独立几何工具层 / tools/geometry_qc；真实求解案例才考虑 examples 或 cases；本轮未定目录',
        'migration_status':'DISCOVERED','evidence':evidence,
        'unknowns':['未执行模块、依赖安装、测试或 HemoCell 接口试接；静态检查不能证明运行兼容性']}
    if rel.startswith(('utils/cfd_lumen/','utils/sampling/','utils/rodent_vasculature/')):record['units']='程序约定位置/半径为 um；fMOST 输入位置先按 spacing 缩放，原始半径物理校准 UNIT_UNVERIFIED'
    if dependency=='specific':record['target_candidate']='REFERENCE_ONLY 或未来独立重写的 HemoCell 应用层；不得放入 core/Palabos core'
    inventory.append(record)

def artifact(path,purpose,cls,units='UNIT_UNVERIFIED',coords='COORDINATE_SYSTEM_UNVERIFIED',evidence=None,experimental='UNVERIFIED',dependency='independent'):
    assert (S/path).exists(),path
    inventory.append({'source_path':str(S/path),'component':path,'purpose':purpose,'entrypoint':False,'solver_dependency':dependency,
        'input':[],'output':[],'units':units,'coordinate_system':coords,'dependencies':[],
        'experimental_relevance':experimental,'migration_class':cls,'modify_required':cls in ('ADAPT','REWRITE'),
        'target_candidate':'待评审的数据/配置候选；本轮不复制','migration_status':'DISCOVERED',
        'evidence':evidence or [{'file':str(S/path),'context':'Current local file/directory exists; see data_file_inventory.json / source_dataset_archives.json','classification':'CONFIRMED_FROM_DATA'}],
        'migration_decision_classification':'INFERRED','unknowns':['最终迁移范围、目录和兼容性尚未确定']})
for c in load('configuration_inventory.json'):
    artifact(c['path'],'已实际解析的配置；默认值不等于实验条件','REFERENCE_ONLY' if 'cfd_flow' in c['path'] or c['path']=='cfd_lumen_config.yaml' else 'ADAPT',
        units='参见单位/物理参数表，禁止自动映射 HemoCell',evidence=[ev(c['path'],kind='CONFIRMED_FROM_CONFIG')],dependency='specific' if 'cfd_flow' in c['path'] else 'independent')
artifact('test_data','ROI003274 小型输入与局部/全局映射 fixture','DATA_ONLY','内部数组显式 _um；原始物理校准仍待核实')
for path,purpose in [('outputs/sampling','ROI 库、节点边映射、裁切端口、代表选择记录'),('outputs/rodent_vasculature','原始与分析 SWC、标注 mask、图结构、来源清单'),(model,'当前本地成功重建的原始管腔 um/m 几何、H5、QC'),(surface+'/geometry','带实体及端口标签的派生封帽表面'),(surface+'/boundaries','四个端口 STL 与实体映射清单'),(surface+'/bc','原始及延长段修正边界参数'),('outputs/reference_validation','半径保真及历史验证证据')]:
    artifact(path,purpose,'DATA_ONLY' if '/bc' not in path and 'reference_validation' not in path else 'REFERENCE_ONLY','代码/记录区分 um 与 m；必须按具体文件取单位')
for path,purpose in [(base,'接受的旧 Base 重启和数值验证凭证'),(prom,'已存在的 promotion replay 与可视化验收记录'),('outputs/documentation','教学或说明性派生几何，不能当真实 CFD 输入'),('references','论文、讲稿、图示及历史说明，部分内容过期')]:artifact(path,purpose,'REFERENCE_ONLY',dependency='specific' if 'cfd_flow' in path else 'independent')
for path in load('source_dataset_archives.json')['families']:
    artifact(path,'原始/公开数据集候选，逐格式保留；归档成员未解压，不能声称已可直接使用','DATA_ONLY',experimental='EXPERIMENTAL_ASSET candidate; dataset provenance/calibration not globally verified')
for path in ['.git','.codex_tmp','.pytest_cache','.ruff_cache','__pycache__','tmp','Ultraliser/build-wsl','vessel_model/T - A high-resolution dataset of mouse brain vasculature/BVLab-Annotation','vessel_model/T - NNE2/MCRInstaller.exe','vessel_model/T - NNE2/NNE2.exe']:
    if (S/path).exists():artifact(path,'版本元数据、缓存、编译产物或打包运行时；不进入 HemoCell 源码','DO_NOT_MIGRATE')
artifact('Ultraliser','外部几何依赖源码身份参考；当前子仓库 dirty，不能声称完全官方未修改','REFERENCE_ONLY')
artifact('external_reference/LBPM','独立旧多相 LBM 源码参考，未发现正式 s1–s5 链调用；不进入 HemoCell core','REFERENCE_ONLY',dependency='specific')
artifact('patches','旧 Musubi/Seeder 补丁及算法历史','REFERENCE_ONLY',dependency='specific')
dump('migration_inventory.json',inventory)

units=[]
def unit(quantity,value,u,file,needle,meaning,confidence='CONFIRMED_FROM_CODE'):
    e=ev(file,needle,confidence);units.append({'quantity':quantity,'value':value,'unit':u,'source_file':e['file'],'source_line_or_context':str(e['line'])+': '+e['context'],'meaning':meaning,'confidence':confidence,'evidence':e})
unit('原始 fMOST SWC 坐标','原始三列 × diag(1,1,2)','raw voxel index -> code-contract um','utils/rodent_vasculature/swc_io.py','points_um=points_voxel','程序把 xyz 当像素坐标乘配置；不是文件自证的物理标定')
unit('体素 spacing',[1,1,2],'um（配置）','configs/swc_roi_generate.yaml','spacing_xyz_um','配置声明；代表 TIFF ResolutionUnit=1、无长度标定，实验依据 UNIT_UNVERIFIED','CONFIRMED_FROM_CONFIG')
unit('原始 SWC 半径','第六列原样保留','程序称 um；物理单位 UNIT_UNVERIFIED','utils/rodent_vasculature/swc_io.py','radius_raw_um=array','没有乘体素 spacing；不能因变量后缀认定原始单位')
unit('ROI 范围',[80,80,120],'um（配置）','configs/swc_roi_generate.yaml','size_um:','采样几何参数，不是 PDMS 尺寸','CONFIRMED_FROM_CONFIG')
unit('重采样间距',1.0,'um（配置）','configs/swc_roi_generate.yaml','resample_step_um','当前 smooth_centerlines=false；参数存在不等于每步启用','CONFIRMED_FROM_CONFIG')
unit('重建半径系数',0.91,'dimensionless','configs/swc_stl_model_generate.yaml','radius_scale:','数值重建进料补偿，仅 H5 feed；原始半径保留','CONFIRMED_FROM_CONFIG')
unit('H5 第四列','2 * source_radius * 0.91','diameter_um','utils/cfd_lumen/ultraliser_backend.py','2.0 * float(feed_radii','实际 H5 属性也写明 diameter_um')
unit('Ultraliser 体素分辨率',6.0,'voxels/um','configs/swc_stl_model_generate.yaml','voxels_per_micron','重建工具参数，不是采集分辨率','CONFIRMED_FROM_CONFIG')
unit('STL 长度变换',1e-6,'um -> m','utils/cfd_lumen/ultraliser_qc.py','vertices=np.asarray','只证明这条生成链的比例；不能推广到任意 STL')
unit('三角面目标边长',0.25913916380971913,'um','configs/cfd_surface_prepare.yaml','target_edge_length_um','派生表面重网格目标','CONFIRMED_FROM_CONFIG')
unit('代表 ROI 半径范围',[1.0,2.7279],'um（保存元数据）',model+'/input/metadata.json','source_radius_min_um','派生直径为 2–5.4558 um；不是已验证的最小流体通道或 RBC 可通行性','CONFIRMED_FROM_DATA')
unit('一维入口平均速度',0.7,'mm/s','configs/cfd_preprocess.yaml','mean_velocity_mm_s','旧结构根节点边界假设；EXPERIMENTAL 未证明','CONFIRMED_FROM_CONFIG')
unit('一维入口体积流量',7.693508475538942e-16,'m3/s',pre+'/qc/run_summary.json','cut_port_inlet_total_m3_s','历史 ROI 流量，不同于后期固定 CFD target','CONFIRMED_FROM_DATA')
unit('后期 CFD 体积流量',2.7369132390905703e-15,'m3/s','configs/cfd_flow.yaml','target_volume_flow_m3_s','旧 Tau1 baseline target；不是本轮测量','CONFIRMED_FROM_CONFIG')
unit('后期 CFD 质量流量',2.890180380479642e-12,'kg/s','configs/cfd_flow.yaml','target_mass_flow_kg_s','旧 CFD 参数','CONFIRMED_FROM_CONFIG')
unit('出口表压',[14.544978101274268,132.20454922317552,-13.700626673311461],'Pa','configs/cfd_flow.yaml','outlet_gauge_pressures_pa','延长段修正/数值基线；负表压不是绝对负压','CONFIRMED_FROM_CONFIG')
unit('结构叶压力参考',0.0,'Pa gauge','configs/cfd_preprocess.yaml','leaf_pressure_pa','结构叶零表压假设','CONFIRMED_FROM_CONFIG')
unit('密度',1056.0,'kg/m3','configs/cfd_flow.yaml','density_kg_m3','旧均匀牛顿流体参数，不是当前样本测量','CONFIRMED_FROM_CONFIG')
unit('运动黏度',3.27e-6,'m2/s','configs/cfd_flow.yaml','kinematic_viscosity_m2_s','旧求解器/模型物性','CONFIRMED_FROM_CONFIG')
unit('动力黏度',0.00345312,'Pa s',pre+'/qc/run_summary.json','dynamic_viscosity_pa_s','历史记录：rho * nu，不是新增 HemoCell 赋值','CONFIRMED_FROM_DATA')
unit('体积黏度',2.18e-6,'m2/s','configs/cfd_flow.yaml','bulk_viscosity_m2_s','旧 LBM 参数','CONFIRMED_FROM_CONFIG')
unit('CFD 网格间距',2e-7,'m','configs/cfd_flow.yaml','dx_m:','旧 Base 0.2 um 笛卡尔格，不是 HemoCell 网格决策','CONFIRMED_FROM_CONFIG')
unit('CFD 时间步',2.038735983690112e-9,'s','utils/cfd_flow/validated_contract.py','return float(dx_m) ** 2','由 dx²/(6ν) 推导的旧 Tau1 数值步长','INFERRED')
unit('数值压力偏置',3387510.7199999993,'Pa','utils/cfd_flow/validated_contract.py','Return the LBM numerical pressure offset','rho * cs² * (dx/dt)²；不能解释为生理压力','INFERRED')
unit('历史时间步',2.44140625e-8,'s','utils/cfd_flow/validated_contract.py','HISTORICAL_DT_S','明确 regression-only，不能与当前 Tau1 混用')
nne=load('nne2_acquisition_metadata.json')[0]
unit('NNE2 代表采集 XY 像素间距',1.16279069767442,'um/pixel',nne['path'],'micronsPerPixel_XAxis','显式采集 XML 元数据；只适用于该栈，不能用于 fMOST','CONFIRMED_FROM_DATA')
unit('通用 STL/NIfTI/Schmid/归档内部单位','UNVERIFIED','UNIT_UNVERIFIED','utils/schmid_pkl/loader.py','coordinates =','部分变量有 _um 后缀，但独立源单位/压力/时间证明不足','UNVERIFIED')
dump('units_inventory.json',units)
physics=[]
for u in units:
    if any(k in u['quantity'] for k in ['速度','流量','压力','表压','密度','黏度','时间步','网格间距']):
        physics.append({'parameter':u['quantity'],'old_value':u['value'],'old_unit':u['unit'],'old_meaning':u['meaning'],'source_file':u['source_file'],
            'source_context':u['source_line_or_context'],'likely_origin':'OLD_SOLVER_PARAMETER','status':u['confidence'],'migration_to_hemocell':'NOT_PERFORMED; requires separate physical and lattice contract'})
dump('physics_parameters.json',physics)

coords=[
 ('origin','活动 fMOST 链没有显式平移；局部 ROI 保留全局坐标。物理世界原点未知。','CONFIRMED_FROM_CODE','utils/sampling/roi_extraction.py','positions = np.asarray'),
 ('X/Y/Z','SWC 列 3/4/5 -> xyz；TIFF 数组为 zyx；解剖轴含义未证明。','CONFIRMED_FROM_CODE','utils/rodent_vasculature/tiff_io.py','array_axis_order'),
 ('axis direction / handedness','缩放矩阵 diag(1,1,2) 不主动翻轴；输入是右手/左手、是否 LPS/RAS 均 COORDINATE_SYSTEM_UNVERIFIED。','UNVERIFIED','utils/rodent_vasculature/swc_io.py','points_um=points_voxel'),
 ('rotate / translate / center','正式分区代码保留几何位置并写 translation_applied=False；可视化相机旋转不能当几何变换。','CONFIRMED_FROM_CODE','utils/cfd_flow/geometry.py','translation_applied'),
 ('scale','fMOST 像素乘 spacing；后续 um -> m 乘 1e-6；这两次缩放不可漏做或重复做。','CONFIRMED_FROM_CODE','utils/cfd_lumen/ultraliser_qc.py','vertices=np.asarray'),
 ('legacy image/world','main.py 的旧 STL/voxel 分支假设输入 LPS，write_nifti_mask 使用 X/Y 取反 affine 写 RAS 和 micron。输入 LPS 本身未被任意 STL 头证明。','CONFIRMED_FROM_CODE','utils/io.py','def lps_to_ras_affine'),
 ('NNE2','TIFF stack 使用 vdb/XML spacing；XY 缩放经取整后保存 actual_x/actual_y；未推广为统一解剖坐标。','CONFIRMED_FROM_CODE','utils/nne2/stack_io.py','actual_x ='),
 ('HemoCell lattice mapping','轴选择、原点/边距、dx、STL 缩放策略、端口面到格点边界尚未建立。','UNVERIFIED','configs/cfd_flow.yaml','dx_m:')]
dump('coordinate_inventory.json',[{'item':k,'meaning':v,'classification':c,'evidence':ev(f,n,c)} for k,v,c,f,n in coords])

risks={
 'UNIT_RISK':('HIGH','fMOST spacing 和原始半径单位主要来自代码/配置；泛用 STL 与代表 NIfTI 文件缺少单位证明；um/m 副本共存。'),
 'COORDINATE_RISK':('HIGH','活动 xyz/zyx 与旧 LPS→RAS 分支不同；解剖轴、原点及 HemoCell 格点变换未确定。'),
 'BOUNDARY_LABEL_RISK':('HIGH','标签存在，但角色是假设；STL 不携带 VTP 的 CellEntityIds/port_id，HemoCell 映射尚未实现。'),
 'GEOMETRY_TOPOLOGY_RISK':('MEDIUM','历史 QC 通过不等于本轮重做拓扑验证，也不等于 HemoCell/RBC 分辨率、可通行性和边界可用性通过。'),
 'SOLVER_COUPLING_RISK':('HIGH','Musubi Lua/LSB、D3Q19 PDF、控制器、Tau1 固定契约不能直接用于 HemoCell；默认 FRESH_STEADY 实际被代码拒绝。'),
 'EXPERIMENTAL_PROVENANCE_RISK':('HIGH','存在采集图像和 XML，但没有确认 PDMS、样本实测入口波形/流量、动静脉身份或原始 SWC 半径校准。'),
 'DEPENDENCY_RISK':('HIGH','VMTK Windows 路径在当前 WSL 不成立；requirements 与 /usr/bin/python3 环境不同；Git LFS 缺失；嵌套仓库 dirty。'),
 'DATA_DUPLICATION_RISK':('HIGH','19 个 SWC 和 19 个 ROI 批次、不同模型/边界版本并存；同名原始 image 相同而 mask/SWC 不同，不能按名字去重。')}
dump('risk_register.json',[{'risk':k,'level':v[0],'reason':v[1],'classification':'INFERRED'} for k,v in risks.items()])

assets=[
 {'asset':'fMOST image / mask / SWC 原始样本','path':str(S/'vessel_model/T - A high-resolution dataset of mouse brain vasculature/raw_data/analysis_data/analysis_data'),
  'asset_type':'EXPERIMENTAL_ASSET','classification':'CONFIRMED_FROM_DATA','unit_evidence':'TIFF 192³，ResolutionUnit=1，无物理标定；SWC 无单位头；配置 spacing 不能代替采集证据','provenance_status':'UNVERIFIED acquisition/sample/calibration linkage','migration_class':'DATA_ONLY'},
 {'asset':'NNE2 microscopy stacks and acquisition XML','path':str(S/'vessel_model/T - NNE2/hana_stk'),
  'asset_type':'EXPERIMENTAL_ASSET','classification':'CONFIRMED_FROM_DATA','unit_evidence':'67 个 XML 有 micronsPerPixel；见 nne2_acquisition_metadata.json','provenance_status':'Specific metadata confirmed; registration/biological identity not fully audited','migration_class':'DATA_ONLY'},
 {'asset':'Schmid 网络和 RBC 轨迹字典','path':str(S/'vessel_model/Schmid/NW1_results'),
  'asset_type':'SCIENTIFIC_REFERENCE_ASSET','classification':'CONFIRMED_FROM_DATA','unit_evidence':'UNIT_UNVERIFIED; no unpickling performed','provenance_status':'不能把 pressure/flow/RBC_trajectories 当本地实测数据','migration_class':'REFERENCE_ONLY'},
 {'asset':'PDMS geometry / channel dimensions / measured inlet waveform','path':None,'asset_type':'EXPERIMENTAL_ASSET','classification':'UNVERIFIED',
  'status':'NOT_FOUND in authored text/path and representative metadata audit','unit_evidence':'UNIT_UNVERIFIED','provenance_status':'PDF contents and compressed archive internals were not exhaustively searched'}]
dump('experimental_assets.json',assets)

# Data-family profiles supply real code/config relationships without falsely assigning units to every file.
profiles=[]
for family in load('data_family_summary.json'):
    f=family['family'];code=[];unittext='UNIT_UNVERIFIED';coordinate='COORDINATE_SYSTEM_UNVERIFIED';value='REFERENCE_ONLY'
    if f=='vessel_model/T - A high-resolution dataset of mouse brain vasculature':code=[ev('utils/rodent_vasculature/catalog.py'),ev('utils/rodent_vasculature/swc_io.py','points_um=')];unittext='fMOST active sample uses configured voxel→um mapping; raw calibration UNIT_UNVERIFIED';value='DATA_ONLY'
    elif f=='vessel_model/T - NNE2':code=[ev('utils/nne2/catalog.py','loadmat'),ev('utils/nne2/stack_io.py','micronsPerPixel_XAxis')];unittext='Per-stack XML/vdb evidence; no global uniform scale';value='DATA_ONLY'
    elif f=='vessel_model/Schmid':code=[ev('utils/schmid_pkl/loader.py','def load_schmid_input')];value='DATA_ONLY, numerical fields REFERENCE_ONLY'
    elif f in ('outputs/sampling','test_data'):code=[ev('utils/cfd_lumen/roi_io.py'),ev('utils/sampling/sampling_io.py')];unittext='Saved _um arrays and local/global IDs; raw calibration still bounded';value='DATA_ONLY'
    elif f=='outputs/rodent_vasculature':code=[ev('utils/rodent_vasculature/pipeline.py','def run_rodent')];unittext='voxel xyz + program-contract um; auxiliary TIFF zyx';value='DATA_ONLY'
    elif f=='outputs/model_generate' or f=='outputs/reference_validation':code=[ev('utils/cfd_lumen/ultraliser_qc.py','def write_geometry_outputs')];unittext='specific generated STL um versus m; H5 diameter_um';value='DATA_ONLY / QC REFERENCE_ONLY'
    elif f=='outputs/cfd_surface_prepare':code=[ev('utils/cfd_surface_prepare/vmtk_pipeline.py','def run_vmtk'),ev('utils/cfd_flow/geometry.py','def partition_surface')];unittext='specific um VTP and STL plus m copy; preserve tags and sidecars';value='DATA_ONLY / boundary physics REFERENCE_ONLY'
    elif f=='outputs/cfd_preprocess':code=[ev('utils/cfd_preprocess/pipeline.py','def run_cfd_preprocess')];unittext='Pa, m3/s, um named and calculated in code';value='REFERENCE_ONLY, port geometry ADAPT'
    elif f=='outputs/cfd_flow':code=[ev('utils/cfd_flow/pipeline.py','def run_cfd_flow'),ev('utils/cfd_flow/restart_decode.py')];unittext='mixed physical SI, lattice units and old restart-specific metadata';value='REFERENCE_ONLY; not importable HemoCell restart'
    elif f in ('scratch','bundled_annotation_runtime','Ultraliser/build-wsl'):value='DO_NOT_MIGRATE'
    profiles.append({'family':f,'format_counts':family['format_counts'],'representatives':family['examples'],'referenced_by_code':code,
        'unit_evidence':unittext,'coordinate_evidence':coordinate,'migration_value':value,
        'confidence':'CONFIRMED_FROM_DATA','migration_decision_classification':'INFERRED','unknowns':['No claim that all family members share the active sample provenance or units']})
dump('data_family_profiles.json',profiles)

# Imports are read statically; installing or importing source-project modules is outside this audit.
dependencies=[];import_to_dist={'PIL':'Pillow','yaml':'PyYAML','skimage':'scikit-image','vtk':'vtk','numpy':'numpy','scipy':'scipy'}
allimports=collections.defaultdict(list)
for r in idx:
    tree=ast.parse((S/r['path']).read_text())
    for n in ast.walk(tree):
        names=[x.name.split('.')[0] for x in n.names] if isinstance(n,ast.Import) else [n.module.split('.')[0]] if isinstance(n,ast.ImportFrom) and n.module and n.level==0 else []
        for name in names:allimports[name].append({'file':r['path'],'line':n.lineno})
for line in (S/'requirements.txt').read_text().splitlines():
    if not line.strip():continue
    name=re.split(r'[<>=!~]',line)[0];roots=[k for k in allimports if import_to_dist.get(k,k)==name]
    dependencies.append({'dependency':name,'declared_constraint':line,'category':'Python package','evidence':sum((allimports[k] for k in roots),[])[:12],
        'installed_version':'UNVERIFIED except explicitly noted /usr/bin/python3 capability check','installed_this_audit':False})
for name,category,context in [('VMTK','mesh tool','tools/run_vmtk_flowextension.py; configs/cfd_surface_prepare.yaml Windows Python/prefix'),('Ultraliser ultraVessMorpho2Mesh','external binary','utils/cfd_lumen/ultraliser_backend.py; nested Git commit 3e4b0eee... and dirty worktree'),('CMake / GNU C++ compiler / C++17','system library','Ultraliser/CMakeLists.txt:27/33/54; CMake >=3.5, GNU compiler check >=9.4, C++17; comment gcc 8.4 is stale'),('OpenMP / TIFF / HDF5 / Eigen3 / GLM / FMT / ZLIB / BZip2','system library','Ultraliser/CMakeLists.txt:41-48 includes corresponding dependency modules; no build invoked'),('APES Seeder / Musubi / Treelm / harvesting / MPI','old solver','configs/cfd_flow.yaml; paths exist, not executed or rebuilt'),('VTK/PyVista + rendering backend','visualization tool','doc_visualize_v*.py; plotting requires display/offscreen policy'),('nibabel / TIFF-Pillow','medical-imaging tool','utils/io.py; rodent_vasculature/tiff_io.py; units need distinct contracts'),('Git LFS','external binary','Native git status failed because git-lfs not found')]:
    dependencies.append({'dependency':name,'category':category,'context':context,'installed_this_audit':False,'classification':'CONFIRMED_FROM_CODE' if name!='Git LFS' else 'CONFIRMED_FROM_DATA'})
dump('dependency_inventory.json',dependencies)

first={'name':'READ_ONLY_TRIANGLE_TOPOLOGY_QC','source_module':str(S/'utils/mesh/quality.py'),
    'functions':['polydata_arrays','edge_topology','duplicate_triangle_count'],
    'input':'Already-selected triangular VTP or in-memory faces (N,3) integer vertex IDs; use a separate VTK reader without modifying source data',
    'output':'boundary-edge count, non-manifold-edge count, duplicate-triangle count, file identity and explicit unit/coordinate provenance status',
    'units':'Dimensionless topology counts; no physical unit inference or conversion','coordinate_relation':'Connectivity-only counts are independent of xyz orientation/translation/length scale',
    'migration_class':'ADAPT','classification':'INFERRED','candidate_target':'tools/geometry_qc/ or separate py_scripts/geometry_qc.py; candidates only',
    'core_changes_required':False,'limitations':['Not a full watertight/self-intersection/normal/volume proof','STL reader may weld vertices; freeze reader semantics','Full measure_mesh_quality _um outputs must wait for explicit unit metadata'],
    'future_tests':['closed tetrahedron','one removed face','three faces sharing one edge','duplicate triangle with reversed index order','coordinate transform leaves counts unchanged'],
    'executed_or_migrated':False,'next_action':'review audit results'}
dump('first_minimal_migration_unit.json',first)

# Existence is confirmed from current files; physical meaning remains separately bounded.
data_rows=load('data_file_inventory.json')
def data_example(suffix):
    return next(r['relative_path'] for r in data_rows if r['relative_path'].endswith(suffix))
features=[]
for name,path,meaning in [
 ('lumen surface',model+'/geometry/lumen_surface_um.vtp','生成链声明 um；历史 QC 证据，本轮未重建'),
 ('STL',model+'/geometry/lumen_surface_m.stl','同链 m 副本；任意其他 STL 不继承此单位'),
 ('surface mesh',surface+'/geometry/cfd_surface_vmtk_tps_boundarynormal_crossseam_um.vtp','实读带 CellEntityIds/port_id 的三角表面'),
 ('volume mesh',prom+'/flow/production_steady_flow_field.vtu','旧规则中心点场导出为六面体 VTU；代码见 steady_export.py，未当通用 HemoCell 输入'),
 ('centerline',data_example('analysis_swc_single_component.npz'),'SWC 中心线/父子关系的保存表示，非实测速度'),
 ('mask',next(r['path'] for r in headers if '/mask/' in r['path']),'代表 TIFF mask 为 192³；物理标定未证'),
 ('voxel',data_example('teaching_voxel_domain.vti'),'教学派生 VTI；不是当前正式血管求解格点'),
 ('inlet surface',surface+'/boundaries/inlet.stl','显式 patch，对应 ID4；角色为结构方向假设'),
 ('outlet surface',surface+'/boundaries/outlet_01.stl','共 outlet_01/02/03 三份，ID3/5/2'),
 ('wall surface',data_example('geometry_solver_m/wall.stl'),'旧 CFD 分区输出；由 VTP 非端口实体提取'),
 ('ROI',data_example('test_data/sampling/roi_library/raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274.npz'),'代表 ROI 有 109 节点/108 边；local/global 映射保留'),
 ('branch labels',data_example('graphs/directed_branches.csv'),'图分支标签保存为 CSV/GraphML/VTP；不是动静脉生理标签')]:
    assert (S/path).is_file(),path
    features.append({'feature':name,'existence':'FOUND','classification':'CONFIRMED_FROM_DATA','path':str(S/path),'meaning':meaning})
features.append({'feature':'target region','existence':'UNVERIFIED','classification':'UNVERIFIED','path':None,'meaning':'已确认计算采样 ROI；未确认独立实验靶区、PDMS 区域或生物学靶区配准'})
dump('geometry_feature_inventory.json',features)
boundary_rows=list(csv.DictReader((S/surface/'boundaries/boundary_manifest.csv').open()))
dump('boundary_labels.json',{'source':str(S/surface/'boundaries/boundary_manifest.csv'),'classification':'CONFIRMED_FROM_DATA','records':boundary_rows,'physiological_role_status':'UNVERIFIED; code/config assign structural ASSUMED_INLET/ASSUMED_OUTLET','vtp_header_evidence':'representative_data_headers.json'})
dump('historical_validation_scope.json',{'classification':'CONFIRMED_FROM_DATA','evidence':[ev('docs/CFD_FLOW.md','## Current scientific scope','CONFIRMED_FROM_DATA'),ev(prom+'/qc/run_summary.json','VALIDATED_BASE_PROMOTION_REPLAY','CONFIRMED_FROM_DATA')],'reported_scope':'Base/Coarse accepted; two-grid sensitivity only; Fine steady and formal three-grid GCI incomplete; WSS deferred','this_audit_reran_validation':False})

entryrows=[
 ['s1_swc_roi_generate.py:148','load_swc_roi_yaml_config → run_rodent_vasculature_pipeline → run_sampling_from_rodent_run','原始 SWC/TIFF/Mask → analysis SWC/图 → ROI NPZ、CSV','CONFIRMED_FROM_CODE'],
 ['s2_swc_stl_model_generate.py:70','load_swc_stl_yaml_config → resolve_sampling_run/select ROI → run_ultraliser_reconstruction','保存 ROI → SWC/H5 → Ultraliser → STL/VTP/QC','CONFIRMED_FROM_CODE'],
 ['s3_cfd_1D_data_preprocess.py:42','load_cfd_preprocess_config → run_cfd_preprocess → solve_global_flow / transfer_all_boundaries','全局图+ROI → 一维压力/流量及假定端口','CONFIRMED_FROM_CODE'],
 ['s4_cfd_surface_prepare.py:39','load_surface_prepare_config → run_vmtk_surface_prepare → local_cut / VMTK / tagged caps','局部切面+延长+重网格+标签及边界侧车','CONFIRMED_FROM_CODE'],
 ['s5_cfd_flow_solve.py:17','load_cfd_flow_config → run_cfd_flow → promotion replay','FRESH_STEADY 分支实际抛授权错误；replay 可解析旧输出并可选 smoke','CONFIRMED_FROM_CODE'],
 ['main.py:170','run_pipeline / run_hierarchical_graph_pipeline','旧 STL 清理→体素→中心线图；默认 unprocess_stl 路径不存在','CONFIRMED_FROM_CODE'],
 ['doc_visualize_v1/v2/v3.py','main / generate_visualizations / render_*','保存结果可视化；v3 明确依赖旧稳态 VTU/物理截面字段','CONFIRMED_FROM_CODE']]
geometry='原始 analysis SWC + TIFF/Mask（辅助） → xyz×spacing、保留原始半径 → 单连通 analysis SWC/层级图 → ROI 裁切及 CUT_PORT/TRUE_TERMINAL → ROIRecord/NPZ → canonical SWC + diameter H5（feed radius×0.91） → Ultraliser → lumen_surface_um.STL/VTP + lumen_surface_m.STL → 一维端口传递 → VMTK 局部切面/延长/重网格/封帽 → CellEntityIds/port_id VTP + 四端口 STL/CSV/JSON → 旧 APES 分区与笛卡尔 mesh。'
unknowns=['fMOST 原始 SWC 半径是否已是物理 um，以及 [1,1,2] 间距的原始采集/论文对应证据',
 '各数据集的物理原点、轴方向、手性及其相互配准；不能把 legacy LPS/RAS 声明套到全部输入',
 '端口生理身份、实测入口流向/流量/波形、PDMS 通道及标定数据',
 '一维 ROI Q=7.6935e-16 与后期 target=2.7369e-15 的完整物理重标定依据；不能拼接不同阶段参数',
 'HemoCell 的几何单位契约、格点映射、边界实现与 RBC 所需有效流道分辨率',
 '多份兼容采样/模型输出中应冻结哪一个批次；null 自动解析不能代替选择决定',
 'VMTK 在本地 Linux 的可用环境及依赖版本组合；两个嵌套 Git 仓库 dirty 的具体源码/换行差异未深审',
 '通用源 STL/OBJ/NIfTI/Schmid PKL 和未解压 ZIP 的单位、授权、适用性及实验来源',
 '历史表面人工评审是否完成；本轮只读取现存 QC，不重跑自交/网格修复或 solver']
dump('unknown_items.json',[{'item':x,'classification':'UNVERIFIED'} for x in unknowns])

md=f'''# 1. Executive summary

本轮为 **DISCOVERY ONLY**，事实源仅为当前 WSL 的 `{S}`。HemoCell 仅观察结构和已有 Git 状态。未迁移、复制几何、运行源工程脚本、测试、solver、构建、安装、提交或 push。

**CONFIRMED_FROM_CODE / CONFIRMED_FROM_CONFIG / CONFIRMED_FROM_DATA**：当前主线是 s1–s5 分段流程，保存几何、局部/全局拓扑、显式端口标签和旧 APES 数值证据。**INFERRED**：适合先提取只读拓扑检查，再评审单位/坐标/边界契约；不能直接迁旧 CFD 配置。

事实标签：`CONFIRMED_FROM_CODE`=实际代码；`CONFIRMED_FROM_CONFIG`=配置声明而非实验真值；`CONFIRMED_FROM_DATA`=当前文件内容/头/元数据；`INFERRED`=推断/建议；`UNVERIFIED`=证据不足。单位未知使用 `UNIT_UNVERIFIED`，坐标未知使用 `COORDINATE_SYSTEM_UNVERIFIED`。本报告的迁移类别和风险等级均为 INFERRED。

范围：根目录深度 2 起步；发现重要目录后再定向读取。文件系统统计覆盖全树；数据扫描仅格式/大小/关联信息，共 {len(load('data_file_inventory.json')):,} 条数据路径，未加载 28 GiB 数据体。静态解析 {len(idx)} 个主项目 Python 模块；完整读取 7 份配置；只检查代表性 TIFF/NPZ/H5/NIfTI 头和一个 1.64 MB 带标签 VTP；67 个采集 XML 与 ZIP 中央目录只读。PDF 全文、巨型压缩包内容、PKL 反序列化及第三方源码全面审阅未执行。

# 2. Source project identity

- CONFIRMED_FROM_DATA：path/realpath 均为 `{S}`，与 `{T}` 独立。
- remote：`https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`
- branch：`codex/cfd-wall-force-numerics-validated-sync-20260830`
- HEAD：`fd21a850a16d0ba17ef3d521123badd53864a3a2`
- 普通 `git status --short` 因缺少 `git-lfs` 失败，原始错误完整保存于 `source_identity_raw.json`。没有安装 LFS 或修改 Git 配置。
- 只读替代：NUL 格式 status 临时禁用 LFS filter 后有 1620 项，335 项经 index LFS 指针的大小+SHA256 核对为已水合相同内容；剩余 1285 项：1114 修改、12 删除、159 未跟踪。完整列表在 `source_git_worktree_assessment.json`，不能把旧工程称为 clean。
- 这些是既有修改，本轮不 reset/restore/clean。旧入口删除与 s1–s5 新文件并存；结论以文件当前内容为准，不以 GitHub/HEAD 内容替代。
- Ultraliser HEAD `3e4b0eee685adbf513e40720a68fd92e66a34b44`，429 行 tracked status；LBPM HEAD `6d686d354e5b8140841d3601e4c8c0e4e4b77e48`，545 行 tracked status。仅保存身份，不判断这些差异均为算法修改。见 `nested_repository_identity.json`。
- 适用 `/home/lzy/.codex/AGENTS.md` 为空；源工程未发现额外 AGENTS.md。

# 3. Limited project map

CONFIRMED_FROM_DATA：总占用约 28 GiB，{stats['counts']['files']:,} 文件、{stats['counts']['directories']:,} 目录；其中 vessel_model 文件内容约 22.71 GB，outputs 约 3.86 GB，.git 约 2.46 GB（这些是文件逻辑大小，不等同 du 占用）。

{table(['目录/模块','实际角色'],[['s1–s5 / configs / utils','主入口、严格配置、几何/ROI/旧 CFD 模块'],['vessel_model','9 个数据集/归档/运行时混合资产族，不能整目录搬运'],['outputs','19 个全局 SWC 批次、19 个采样批次、3 个模型目录、12 个表面目录、42 个 CFD 目录'],['test_data','ROI003274 最小 fixture，而非完整原始数据'],['Ultraliser','几何外部工具源码及 build-wsl'],['external_reference/LBPM','另一个 LBM 参考仓库；不在当前 s1–s5 正式调用链'],['docs / references / outputs/documentation','说明、论文、教学图与历史证据，部分陈旧'],['.git / .codex_tmp / caches / tmp','不迁移的工程元数据与缓存']])}

# 4. Main entrypoints

{table(['入口','实际调用','功能/限制','证据类别'],entryrows)}

代码逐函数行号、imports、调用目标见 `code_structure_index.json`。入口脚本只被读取，连 `--help` 也未执行。C++ 重建入口位于 `Ultraliser/apps/ultraVessMorpho2Mesh/VessMorpho2Mesh.cpp:133`；它是外部几何程序，不是当前血流求解器。

# 5. Configuration chain

CONFIRMED_FROM_CODE：s1 → `utils/swc_roi_yaml_config.py:load_swc_roi_yaml_config`；s2 → `utils/cfd_lumen/model_yaml_config.py:113`；s3/s4/s5 分别调用各子包 `config.py` 的严格 YAML loader。大部分相对路径按 PROJECT_ROOT 解释，s2 的用户 YAML 参数先按 cwd 解析。源码使用 `yaml.safe_load`，校验字段/类型/互斥 ROI selector。

CONFIRMED_FROM_CONFIG：s2 和 s3 指向 `outputs/sampling/20260825_133201_radius_plus_structure_k5`，anchor=3274；s3 的 model_run/rodent_run 为 null，可按兼容性自动解析；s4 和 s5 固定特定旧输出路径。模型目录最新并不等于正式 CFD 输入。`configuration_inventory.json` 保存所有键值、SHA256 和具体路径是否存在。

INFERRED：顶层 `cfd_lumen_config.yaml` 属于旧通用配置入口参考，s2 实际读取完整 `configs/swc_stl_model_generate.yaml`；不能把两个配置当一份。VMTK 配置中的 D:/anaconda3 路径与 ../external/vmtk 在当前 WSL 不成立。当前配置 FRESH_STEADY 在 `utils/cfd_flow/pipeline.py:346` 的实际分支抛出 `CFD_FLOW_FRESH_STEADY_REQUIRES_EXPLICIT_COMPUTE_AUTHORIZATION`；不能声称默认 CLI 可直接跑完整新流场。

# 6. Input data

CONFIRMED_FROM_DATA：`data_file_inventory.json` 对每个匹配路径记录格式、字节、推测角色、代码/配置关联、单位/坐标证据和迁移价值；`data_family_profiles.json` 给出族级证据，族内未读文件不会自动继承代表样本单位。

- fMOST 活动样本：`.../raw_data/analysis_data/analysis_data/{{images,mask,swc}}/fMOST_0_5_6_0_0_6_0001_02_01.*`。image/mask 为 192³ TIFF，SWC 有七列；没有长度单位头。原图与 total_vascular_data 同名图像 SHA 相同，但 mask 和 SWC 不同，不可按名字去重。
- `test_data`：analysis NPZ 有 7419 节点；ROI NPZ 有 109 节点、108 边、3 CUT_PORT、1 TRUE_TERMINAL，保留 local/global IDs 和坐标/半径数组。它是测试输入，不是独立实测流场。
- 当前成功模型 H5：points=(113,4)，structure=(5,2)，connectivity=(4,2)，属性明确 coordinate_unit=um、第四列 diameter_um、radius_scale=0.91。
- 实际格式还包括 STL/VTP/VTU/VTI、CSV/JSON/YAML、TIFF/PNG/JPG、NIfTI、NPZ、H5、MAT、SWC、PKL、GraphML、旧 Lua/LSB/RES；泛用 NPY 示例出现在打包运行时，不要误认活动血管输入。OBJ 来自 Ultraliser 示例。
- `.msh/.ply/.dcm` 和未压缩传统 `.vtk` 在当前格式扫描中 NOT_FOUND；ZIP 内成员另记，不能当已解压输入。AneuX/Brain_arteries/IntrA 等存在压缩归档，中央目录见 `source_dataset_archives.json`，本轮未解压。
- 代表 NIfTI 头 `xyzt_units=0`、qform=0，虽有 sform/pixdim，仍不能猜 mm 或 um。Schmid PKL 有网络和 RBC_trajectories 文件，未反序列化，单位及实验/模拟来源 UNVERIFIED。

# 7. Geometry generation pipeline

CONFIRMED_FROM_CODE：{geometry}

原始 mask 是辅助核查/可视化数据，不会自动重写 canonical SWC 半径和父子边；s2 依赖保存的 ROIRecord，不是直接从任意 TIFF 重建。`utils/cfd_lumen/ultraliser_backend.py:248` 将 H5 radius feed 乘 0.91 并以直径写第四列；`ultraliser_qc.py:40` 导出 um 与 m 副本。ROI 构建保留全局物理位置，不把 local node ID 当局部坐标。

实际存在：lumen surface、STL/VTP surface mesh、centerline/SWC/GraphML、mask/TIFF、voxel/教学 VTI、ROI、branch labels、显式入口/出口 patch。旧 volume mesh 是 Seeder 笛卡尔 LSB 和导出的六面体 VTU；不是 HemoCell 可直接加载的通用四面体网格。原始未封帽、修补与封帽版本均有历史记录，不能仅取“最新 STL”。

逐项存在性、完整路径与用途限制见 `geometry_feature_inventory.json`。计算 ROI 已确认；独立实验 target region 为 UNVERIFIED。CONFIRMED_FROM_CODE / CONFIRMED_FROM_DATA：`docs/CFD_SURFACE_PREPARE.md` 描述旧环形延长/平滑算法，与当前 s4 的 VMTK 调用链不同，归为历史参考。选中历史表面是 cap-only recovery 产物，复用了已有 open surface；不能声称已由本轮当前代码完整重跑。

# 8. Output pipeline

CONFIRMED_FROM_CODE / CONFIRMED_FROM_DATA：各阶段按输出批次写 source config、manifest、geometry/data、QC、figures/report。s1 写 rodent_vasculature 和 sampling；s2 写 input/SWC/H5 与 geometry/STL/VTP、qc、report；s3 写 global_1d、roi/boundary_conditions.json 和 port planes；s4 写 geometry、boundaries、bc、vmtk、qc；s5 的 replay 解析已有 Base/LSB，输出稳态物理场 VTU、metrics CSV、QC JSON 和 production_review.html。

历史 `production_tau1_base_promotion_anchor003274_20260902_013637/qc/run_summary.json` 的状态是 `CFD_FLOW_PRODUCTION_TAU1_INTEGRATION_AND_VISUAL_REGRESSION_PASS`，execution_mode=`VALIDATED_BASE_PROMOTION_REPLAY`。这是读取到的历史记录，**不是本轮验证通过或新鲜求解**。docs/CFD_PRODUCTION_INTEGRATION_CONTEXT.md 的“尚未 promotion”描述已经落后于当前代码/数据。

CONFIRMED_FROM_DATA：当前 `docs/CFD_FLOW.md` 的 scientific scope 记录 Base/Coarse 稳态已接受，但仅为双网格分辨率敏感性；Fine 稳态和正式三网格 GCI 未完成，WSS 为 DEFERRED。本轮未复验，也不据此声称网格无关或壁面剪应力准确。该文档仍使用已删除的旧 CLI 名称，运行入口以当前 s5 源码为准；证据见 `historical_validation_scope.json`。

# 9. Units

{table(['quantity','value','unit','meaning','confidence / source'],[[u['quantity'],u['value'],u['unit'],u['meaning'],u['confidence']+' '+link(u['evidence'])] for u in units])}

UNITS_STATUS：程序内部约定和部分文件单位可追溯；原始 fMOST 物理标定、原始半径、其他 STL/PKL/归档单位 `UNIT_UNVERIFIED`。既不全部宣布未知，也不把变量后缀当原始校准证明。

# 10. Coordinate system and transforms

{table(['item','observed relationship','classification / evidence'],[[k,v,c+' '+link(ev(f,n,c))] for k,v,c,f,n in coords])}

# 11. Boundary / inlet / outlet / wall labels

CONFIRMED_FROM_DATA：正式带标签 VTP 实读 73416 点、67262 cells，cell arrays 含 CellEntityIds、boundary_type_code、boundary_index、boundary_origin_code、boundary_origin、port_id、SurfaceRegionId/SurfaceRegion、RemeshEntityId。实体计数：wall ID1=67071；ID2=49；ID3=44；ID4=56；ID5=42。

boundary_manifest.csv 对应：inlet→ID4（CUT_PORT）；outlet_01→ID3（CUT_PORT）；outlet_02→ID5（CUT_PORT）；outlet_03→ID2（TRUE_TERMINAL）。四个端口 STL 实际存在。旧 CFD partition 从剩余实体推导 wall，并写五份米制 patch。STL 文件自身不能保存这些 VTP 属性，必须成套保留 CSV/JSON/VTP。

CONFIRMED_FROM_CODE / CONFIRMED_FROM_CONFIG：`port_transfer.py:50` 根据 local edge parent/child 方向产生 ASSUMED_INLET/OUTLET；TRUE_TERMINAL 被当作 ASSUMED_OUTLET；`measured_flow_direction=false`。因此这里既不是“只有无标签开口”，也不是“实验确认的入口”。HemoCell 如何使用这些面和 ID 属 UNVERIFIED，不能自动套用旧标签号。完整 CSV 记录见 `boundary_labels.json`。

# 12. Existing fluid / physics parameters

见 `physics_parameters.json`，每项有 old_value、old_unit、old_meaning、source_context、likely_origin、status。该表当前值全部是 `OLD_SOLVER_PARAMETER` 或其数值推导；没有将其标成 EXPERIMENTAL。

主要差别：一维 root=0.7 mm/s、ROI inlet Q=7.693508475538942e-16 m³/s；后期 Tau1 target=2.7369132390905703e-15 m³/s。二者属于不同历史尺度/边界契约，完整重标定来源待评审。rho=1056 kg/m³、nu=3.27e-6 m²/s、mu=0.00345312 Pa·s，不能直接作为 HemoCell plasma/RBC 模型参数。旧 P_ref≈3.3875 MPa 是数值压力偏置；历史 23622.320128 Pa 仅 regression reference。

# 13. Solver-independent components

CONFIRMED_FROM_CODE：SWC/ROI 数据结构、图连通性与 source-edge 映射、中心线/采样特征、网格拓扑、坐标变换、Ultraliser/VMTK 几何处理都没有依赖 HemoCell/Musubi 核心 API。独立于血流 solver 不代表无依赖或可无条件执行：VMTK/Ultraliser 有外部环境，某些 geometry loader 会在内存中合并点，repair/remesh 会生成改动几何。

INFERRED：优先 `utils/mesh/quality.py` 的无量纲拓扑函数；之后评审 `rodent_vasculature/swc_io.py`、`sampling/*`、`cfd_lumen/*`、`cfd_surface_prepare/*`。独立一维网络可以作为边界建模参考，但它是另一个物理模型，不能把结果视为实验真值。

# 14. Solver-specific components

CONFIRMED_FROM_CODE：`utils/cfd_flow/apes.py` 的 Lua 渲染与 WSL 启动、production/pipeline、restart_decode 的 LSB/PDF、Musubi one-step mass replay、D3Q19 link flux、固定 Tau1 数值参数，以及 patches/musubi 和 patches/seeder 均有旧 solver 耦合。

INFERRED：运行入口/输入适配未来应 REWRITE；重启、旧 core 补丁与固定参数 REFERENCE_ONLY；physical_port_flux、steady_state、steady_export 和 doc_visualize_v3 属 mixed，需 ADAPT 场命名/单位/边界契约。不能整套复制进 HemoCell，也不改 Palabos core。

# 15. Experimental / PDMS assets

见 `experimental_assets.json`。fMOST image/mask/SWC 和 NNE2 栈/XML 单列保护；NNE2 67 份采集 XML 有明确 micronsPerPixel。没有在已审阅文本、路径和代表元数据里确认 PDMS 几何、芯片尺寸、用户实测速度/流量/入口波形或 target-region 实验配准。`NOT_FOUND` 仅针对已声明搜索范围，不声称遍查 PDF 或 ZIP 内容。

Schmid 网络/RBC trajectories 文件有科研参考价值，来源和单位待确认，不能与采集图像或当前配置物性混为一类。

Schmid 数据族在清单中为 DATA_ONLY 候选，指后续可选择原始图结构；其压力/流量/RBC 轨迹等历史数值结果为 REFERENCE_ONLY。族级类别不授权整目录迁移。

# 16. Geometry-quality tools

{table(['检查/工具','实际位置与行为','迁移判断'],[['boundary/non-manifold/duplicate faces','utils/mesh/quality.py:97/108/138；VTK/NumPy，counts','ADAPT 最小只读接口'],['degenerate / duplicate vertices','quality.py:80 与 cleanup.py:121 vtkCleanPolyData','ADAPT；必须分开检测和修复'],['normals / connected components / bounds','cleanup.py:284 AutoOrientNormals；quality.py VTK connectivity、bounds/area/volume','ADAPT；面积和体积单位不得硬猜'],['self-intersection','cfd_lumen/ultraliser_qc.py:174；trimesh 空间候选 + VTK triangle intersection；vmtk_qc/guarded_remesh 复用','ADAPT；候选过滤和邻接排除有方法局限'],['radius fidelity / minimum diameter','ultraliser_qc 的截面比较，input/metadata 的 source radius 统计','ADAPT；不是完整最小流道或 RBC passage proof'],['repair / remeshing','pymeshfix cleanup；VMTK TPS、entity remesh、local_cut','ADAPT；未来单独授权几何变更'],['voxelisation / skeleton','utils/voxel VTK image stencil、connectivity、skimage skeleton','ADAPT；不同图像轴/空间单位契约']])}

本轮不执行这些处理或其测试。既有 QC 只作为 CONFIRMED_FROM_DATA 的历史结果，未重新证明网格质量。

# 17. External dependencies

`requirements.txt` 实际列出 NumPy、SciPy、h5py、VTK、PyVista、scikit-image、nibabel、Matplotlib、NetworkX、Jinja2、pymeshfix、psutil、pytest、Pillow、PyYAML、trimesh、Shapely、manifold3d、rtree；约束与 imports 见 `dependency_inventory.json`。静态分析不等于全部依赖已经可用。

CONFIRMED_FROM_DATA：本轮仅探测 `/usr/bin/python3` 的 module spec，numpy/h5py/PIL/yaml/vtk/scipy 可找到，pyvista/nibabel 不可找到；这不代表所有 Python 环境。没有安装依赖。VMTK Windows 环境路径在 WSL 不成立；APES 可执行路径仅 existence 检查，未运行。Ultraliser 是外部几何 C++ 工具，LBPM 是独立参考旧 solver，不能因二者在目录中而称它们都是活动主线。

CONFIRMED_FROM_CODE：Ultraliser 的 CMake 包含 OpenMP、TIFF、HDF5、Eigen3、GLM、FMT、ZLIB、BZip2；要求 CMake≥3.5、C++17，GNU 编译器实际版本检查为≥9.4（注释的 gcc 8.4 已过期）。这里只读构建声明，没有执行 CMake 或编译。

# 18. Migration candidates

`migration_inventory.json` 含 {len(inventory)} 个模块/资产/配置条目，记录 source_path、purpose、entrypoint、solver_dependency、input/output、units、coordinate_system、dependencies、实验相关性、迁移类别、目标候选、证据和 unknowns；均为 DISCOVERED。

{table(['class','候选/范围'],[['KEEP_AS_IS','rodent_vasculature/geometry.py 的独立数学逻辑；sampling/feature_scaling.py 连同 ScalerState'],['ADAPT','SWC/ROI loader、mesh QC、几何/坐标、Ultraliser/VMTK 适配器、边界侧车和混合后处理'],['REWRITE','APES 运行编排/输入渲染的功能接口，未来设计 HemoCell 应用层等价功能'],['REFERENCE_ONLY','旧 solver restart/PDF/core patches/Tau1 参数、历史验证及教学/文档结果'],['DATA_ONLY','经选择的原始图像/标注/ROI、单位明确的几何与标签侧车；原始大归档待评审'],['DO_NOT_MIGRATE','Git 元数据、缓存、build、对象/二进制和打包运行时']])}

目标侧 CONFIRMED_FROM_DATA：当前 HemoCell HEAD `{(A/'target_head_before.txt').read_text().strip()}`；已有 `M cmake/setup_googletest.cmake`、未跟踪验证 build/示例可执行文件及忽略的 build/tmp。它们全部是本轮开始前状态，不清理、不修复。根结构有 scripts/tools/examples/cases/data，没有 py_scripts/test_code。`cases/README.md` 自身称这些 cases 为潜在未验证/未完成案例。INFERRED 候选：只读几何工具可独立 tools/geometry_qc 或未来 py_scripts；真实求解应用才考虑 examples/cases；本轮不确定最终目录、不生成迁移代码。

# 19. DO_NOT_MIGRATE

CONFIRMED_FROM_DATA：.git（含 LFS 对象）、__pycache__/.pytest_cache/.ruff_cache/.codex_tmp/tmp、Ultraliser/build-wsl 的 CMakeFiles/CMakeCache/*.o/二进制、Windows BVLab-Annotation/_internal 运行时和 NNE2 安装器等不应进入新源码。source root 未发现 .venv/venv，但“数据目录”内部确实含打包库和 executable。

不得整目录丢弃 data/results/outputs：accepted Base restart、mesh/qval、物理指标、失败/未完成证据属于 REFERENCE_ONLY 科研凭证；仅区别“是否迁到新工程”，不是建议删除。零碎历史日志根据是否与验收绑定区分；本轮一律不清理。

# 20. Risks

{table(['risk','level','reason'],[[k,v[0],v[1]] for k,v in risks.items()])}

# 21. Unknown / unverified items

'''+''.join('- UNVERIFIED：'+u+'。\n' for u in unknowns)+f'''
未完成事项的含义是“未来迁移的证据缺口”，不是已执行迁移的失败。本轮审计覆盖的文件/方法范围明确记录；最终只读一致性见 `READ_ONLY_VERIFICATION.json`。

# 22. Recommended first minimal migration unit

**FIRST_MINIMAL_MIGRATION_UNIT = READ_ONLY_TRIANGLE_TOPOLOGY_QC**（INFERRED；仅推荐）。

从 `utils/mesh/quality.py:32/97/108` 提取读取三角连接关系和边界边、非流形边、重复面计数的最小接口。输入是明确选择的三角 VTP 或 N×3 face IDs；输出是无量纲计数、输入身份和显式 unknown 单位/坐标元数据。它不改变 mesh，不运行 solver，不依赖实测流量，也不修改 HemoCell/Palabos core。

选择原因：当前代码和输入/输出已可静态确认，拓扑计数不使用物理坐标，因此不受尚未解决的单位/解剖轴影响；容易用四面体、缺面、非流形边、反序重复面验证。完整 `measure_mesh_quality` 的 _um 面积/体积字段应等待显式单位契约，暂不一并照搬。它也不能证明自交、法向、体积或 RBC 可通行性。

候选目标为独立 `tools/geometry_qc/` 或未来 `py_scripts/geometry_qc.py`，尚未选定或创建。详细接口/后续测试建议在 `first_minimal_migration_unit.json`；本轮没有执行这些测试或下一阶段。

NEXT_RECOMMENDED_ACTION = review audit results
'''
(A/'MIGRATION_CONTEXT.md').write_text(md)
print(json.dumps({'migration_inventory_entries':len(inventory),'units':len(units),'physics':len(physics),'dependencies':len(dependencies),'context_bytes':len(md.encode()),'first_minimal_unit':first['name']},indent=2))
