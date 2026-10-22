"""Regenerate only the frozen vascular core and support for P1S / 0.4 / ABS.

Uses the existing s1-11 geometry, native slicing and practical audit kernels.
Historical X1C geometry, casting components and GUI profiles are read-only.
"""
from __future__ import annotations
import copy
import json
import os
import pickle
import shutil
import zipfile
from pathlib import Path
import numpy as np
import trimesh
from . import bambu_manufacturing as native
from .core_support_inputs import default_output as legacy_output, discover_source, freeze, branches
from .core_support_geometry import generate
from .core_support_slicing import evaluate
from .core_support_practical_audit import audit_candidate
from .core_support_practical_geometry import repair_true_missing
from .core_support_practical_topology import topology_candidate
from .core_support_practical_removal import audit_removal, practical_gates
from .support_inputs import save, digest, verify, numeric, profile_values, ROOT

PRESETS = dict(machine='Bambu Lab P1S 0.4 nozzle', process='0.20mm Standard @BBL X1C', filaments=['Bambu ABS Vessel'])

def load(path):
    return json.loads(Path(path).read_text())

def default_output():
    return legacy_output() / 'P1S_04_ABS_20261009'

def validate_profiles(actual):
    m,p,f = (actual[k] for k in ['machine','process','filament'])
    if m.get('printer_model') != 'Bambu Lab P1S' or m.get('name') != PRESETS['machine']:
        raise ValueError('P1S_MACHINE_REQUIRED')
    if numeric(m['nozzle_diameter']) != .4 or len(m['nozzle_diameter']) != 1:
        raise ValueError('P1S_04_NOZZLE_REQUIRED')
    if f.get('filament_type') != ['ABS']:
        raise ValueError('ABS_REQUIRED')
    if p.get('compatible_printers') and m['name'] not in p['compatible_printers']:
        raise ValueError('PROCESS_NOT_COMPATIBLE_WITH_P1S')
    if numeric(p['layer_height']) != .2 or numeric(p['line_width'],.4) != .42:
        raise ValueError('RECORDED_P1S_STANDARD_PROFILE_CHANGED_USE_NEW_RUN_CONFIG')
    return True

def profiles(out,review_only=False):
    if review_only:
        # Reviewing a migrated archive does not require its original Windows
        # executable. Never substitute Linux engine hashes into historical data.
        prior=load(out/'profiles.json')
        validate_profiles(prior['actual'])
        verify({r['path']:r['sha256'] for r in prior['values']['profile_sources']})
        verify({str(out/(kind+'.json')):value for kind,value in prior['hashes'].items()})
        return prior
    old = load(legacy_output()/'final_practical_v3/bambu/profiles/profiles.json')
    roots = [Path(r['path']).parents[3] for r in old['values']['profile_sources'] if '/system/BBL/' in r['path']]
    if os.environ.get('VASCULAR_BAMBU_PROFILE_ROOT'):
        roots.insert(0,Path(os.environ['VASCULAR_BAMBU_PROFILE_ROOT']))
    root = next((p for p in roots if (p/'system/BBL/machine').is_dir()),None)
    if root is None:
        raise FileNotFoundError('BAMBU_PROFILE_ROOT_NOT_FOUND')
    resolved = native.active_profiles(root,PRESETS,'ABS')
    actual = copy.deepcopy(resolved['flattened'])
    validate_profiles(actual)
    original_support = actual['process'].get('enable_support')
    actual['process']['enable_support'] = '0'
    for kind,settings in actual.items():
        settings.update(type=kind,instantiation='true')
    values = profile_values(actual)
    if values['fallback_fields']:
        raise ValueError('P1S_PROFILE_FALLBACK_NOT_ALLOWED')
    p = actual['process']
    values.update(support_line_width=numeric(p['support_line_width'],.4),
        initial_layer_line_width=numeric(p['initial_layer_line_width'],.4),
        roof_gap_mm=numeric(p['support_top_z_distance']),source='EXPLICIT_P1S_04_SAVED_ABS_PROFILE',
        process_name=PRESETS['process'],printer='Bambu Lab P1S',material='ABS',
        profile_sources=resolved['profile_sources'],gap_status='PHYSICAL_REMOVAL_GAP_UNCALIBRATED',
        bridge_limit_mm=2.,bridge_calibration_status='UNCALIBRATED')
    # A byte-identical installed binary has the same version as the recorded
    # successful native runs; no unreliable GUI launch is needed to infer it.
    engine = copy.deepcopy(old['engine'])
    explicit_binary=os.environ.get('VASCULAR_BAMBU_BIN')
    if explicit_binary and not Path(explicit_binary).is_file():
        raise FileNotFoundError('EXPLICIT_BAMBU_BINARY_MISSING: '+explicit_binary)
    if explicit_binary or not Path(engine['binary']).is_file() or digest(engine['binary']) != engine['binary_sha256']:
        from .support_slicer_adapter import discover_engine
        engine = discover_engine('Bambu',dict(binary=explicit_binary or engine.get('binary'),
                                 native_binary=os.environ.get('VASCULAR_BAMBU_NATIVE_BIN')),out/'engine')
    else:
        engine['version_evidence'] = 'Installed binary SHA256 equals previous verified native Bambu binary'
    if engine['status'] != 'ENGINE_FOUND':
        raise RuntimeError('ACTUAL_BAMBU_BINARY_REQUIRED')
    pp=out/'profiles.json'
    if pp.exists():
        prior=load(pp)
        if prior['actual'] != actual or prior['engine']['binary_sha256'] != engine['binary_sha256']:
            raise ValueError('EXISTING_P1S_OUTPUT_PROFILE_CHANGED_CHOOSE_NEW_OUTPUT')
        verify({r['path']:r['sha256'] for r in prior['values']['profile_sources']})
        return prior
    for kind,settings in actual.items():save(out/(kind+'.json'),settings)
    result=dict(values=values,actual=actual,engine=engine,presets=PRESETS,
        profile_root=str(root),user_gui_presets_modified=False,
        overrides={'process.enable_support':dict(before=original_support,after='0',reason='Independent explicitly modeled support')},
        hashes={k:digest(out/(k+'.json')) for k in actual})
    save(pp,result)
    return result

def evaluated(folder,profile_bundle,readonly=False):
    path=folder/'evaluation.json'
    if path.exists() and load(path).get('slice',{}).get('slicer_success'):
        result=load(path)
        signature=dict(project=digest(folder/'bambu/input.3mf'),profiles=profile_bundle['actual'],
                       binary=profile_bundle['engine']['binary_sha256'])
        if load(folder/'bambu/slice_inputs.json')!=signature:
            raise ValueError('CACHED_SLICE_INPUTS_CHANGED_CHOOSE_NEW_OUTPUT')
        with zipfile.ZipFile(result['slice']['archive']) as z:
            if z.testzip() is not None:
                raise ValueError('CACHED_SLICE_ARCHIVE_CORRUPTED')
        return result
    if readonly:
        raise FileNotFoundError('REVIEW_REQUIRES_COMPLETED_EVALUATION: '+str(path))
    for attempt in range(3):
        # Recover an already finished native slice after a Python audit failure.
        cached=folder/'bambu/slice_result.json'
        reuse=cached.exists() and load(cached).get('slicer_success',False)
        result=evaluate(load(folder/'geometry.json'),profile_bundle['actual'],profile_bundle['engine'],reuse_native=reuse)
        if result.get('slice',{}).get('slicer_success'):
            return result
        logs='\n'.join(p.read_text(errors='replace') for p in (folder/'bambu').glob('*.log'))
        if not any(s in logs for s in ['UtilAcceptVsock','UtilConnect','accept4 failed']):
            raise RuntimeError('NATIVE_SLICE_FAILED_SEE_CANDIDATE_LOGS: '+str(folder))
        archive=folder/'runtime_retries'/str(attempt+1)
        shutil.copytree(folder/'bambu',archive,dirs_exist_ok=True)
        print('WSL_INTEROP_RETRY',folder.name,attempt+1,flush=True)
    raise RuntimeError('NATIVE_WSL_INTEROP_FAILED_AFTER_RETRIES')

def score(audit):
    return (audit['structural']['true_missing_support_count'],audit['macro']['macro_fail_count'],
            audit['floating_vascular_regions_unresolved'])

def deliver(out,folder,bundle,history):
    audit,structural,macro=audit_candidate(folder)
    if not (folder/'strategic_cut_plan.json').exists():
        # A rejected topology must not inherit its purported cuts in the
        # selected unmodified parent. Audit that parent with an empty plan.
        save(folder/'strategic_cut_plan.json',dict(weak_links=[],strategic_cuts=[],reason='NO_ACCEPTED_TOPOLOGY_FOR_THIS_CANDIDATE'))
        with (folder/'strategic_cut_solids.pkl').open('wb') as f:pickle.dump([],f)
    removal_path=folder/'practical_removal/removal_audit.json'
    removal=load(removal_path) if removal_path.exists() else audit_removal(folder)
    e=load(folder/'evaluation.json')
    protection=verify(load(out/'validation/protected_before.json'))
    prior=verify(load(out/'validation/historical_protection.json'))
    saved=verify({r['path']:r['sha256'] for r in bundle['values']['profile_sources']})
    g=pickle.load((folder/'geometry.pkl').open('rb'))
    dest=out/'selected';dest.mkdir(exist_ok=True)
    names={'frozen_vascular_core.stl':out/'source/frozen_vascular_core.stl',
        'BG001_RMCA_vascular_core_P1S_04_ABS.stl':folder/'vascular_core.stl',
        'BG001_RMCA_custom_support_P1S_04_ABS.stl':folder/'custom_support.stl',
        'BG001_RMCA_core_support_P1S_04_ABS_project.3mf':folder/'bambu/input.3mf',
        'BG001_RMCA_core_support_P1S_04_ABS_sliced.3mf':Path(e['slice']['archive']),
        'vascular_print_transform.json':folder/'vascular_print_transform.json',
        'inverse_transform.json':folder/'inverse_transform.json'}
    for name,src in names.items():shutil.copy2(src,dest/name)
    archive=dest/'BG001_RMCA_core_support_P1S_04_ABS_sliced.3mf'
    with zipfile.ZipFile(archive) as z:
        crc=z.testzip() is None
        settings=json.loads(z.read('Metadata/project_settings.config'))
        gcode_names=[n for n in z.namelist() if n.endswith('.gcode')]
        for n in gcode_names:(dest/Path(n).name).write_bytes(z.read(n))
        for n in ['Metadata/plate_1.png','Metadata/plate_1_small.png']:
            if n in z.namelist():
                (out/'QC').mkdir(exist_ok=True);(out/'QC'/('native_'+Path(n).name)).write_bytes(z.read(n))
    machine_ok=settings.get('printer_model')=='Bambu Lab P1S' and settings.get('nozzle_diameter')==['0.4'] and settings.get('filament_type')==['ABS']
    archive_valid=crc and len(gcode_names)==1 and machine_ok and e['archive_geometry']['passed'] and e['independent_part_count']==2
    gates=practical_gates(audit,removal,protection['all_unchanged'],archive_valid)
    gates['P1S_04_ABS_PROFILE_CONFIRMED']=machine_ok
    geometry=load(folder/'geometry.json')
    gates['GEOMETRY_ROUTING_ERRORS_ZERO']=not geometry.get('errors',[])
    ready=all(gates.values())
    summary=dict(status='READY_FOR_REAL_PRINT_TRIAL_WITH_PHYSICAL_REMOVAL_UNCALIBRATED' if ready else 'NEEDS_ADJUSTMENT',
        output=str(out),selected_candidate=folder.name,printer=PRESETS['machine'],material='ABS',
        process=PRESETS['process'],layer_height_mm=bundle['values']['layer_height'],line_width_mm=bundle['values']['line_width'],
        nozzle_diameter_mm=.4,scope='VASCULAR_CORE_AND_CUSTOM_SUPPORT_ONLY',
        orientation_source=str(legacy_output()/'final_practical_v3/final_practical_summary.json'),
        anatomy_scaled=False,anatomy_modified=False,support_regenerated=True,box_or_retainers_generated=False,
        profiles=bundle['hashes'],source_protection=protection,historical_protection=prior,saved_profile_protection=saved,
        gates=gates,failures=[k for k,v in gates.items() if not v],practical_audit=audit,
        geometry_routing_errors=geometry.get('errors',[]),
        removal=removal,slice=e['slice'],archive_validation=dict(passed=archive_valid,zip_crc_pass=crc,
        gcode_entries=gcode_names,actual_printer=settings.get('printer_model'),actual_nozzle=settings.get('nozzle_diameter'),
        actual_material=settings.get('filament_type'),part_count=e['independent_part_count'],auto_support_off=e['auto_support_off']),
        history=history,production_frozen=ready,physical_removal_validated=False,no_printer_upload_or_start=True,
        files={str(p.relative_to(out)):digest(p) for p in dest.iterdir() if p.is_file()})
    save(out/'regeneration_summary.json',summary)
    save(out/'validation/source_protection.json',protection)
    save(out/'validation/historical_protection_result.json',prior)
    save(out/'validation/saved_profile_protection.json',saved)
    for p in (folder/'practical').glob('*'):
        if p.is_file():shutil.copy2(p,out/p.name)
    for name in ['safe_removal_sequence.json','removal_audit.json','strategic_cut_seams.csv','final_contact_pieces.csv']:
        p=folder/'practical_removal'/name
        if p.exists():shutil.copy2(p,out/name)
    render(out,g,summary)
    report(out,summary)
    print('P1S_RESULT',summary['status'],summary['failures'],flush=True)
    return summary

def render(out,g,summary):
    import pyvista as pv
    (out/'QC').mkdir(exist_ok=True)
    for name,show_core,show_support in [('assembly',True,True),('vascular_core',True,False),('custom_support',False,True)]:
        plot=pv.Plotter(off_screen=True,window_size=(1400,1050));plot.set_background('white')
        for mesh,color,opacity,show in [(g['vascular'],'#c95949',1.,show_core),(g['support'],'#6699aa',.85,show_support)]:
            if show:plot.add_mesh(pv.wrap(mesh),color=color,opacity=opacity)
        plot.add_axes();plot.show_bounds(xtitle='X (mm)',ytitle='Y (mm)',ztitle='Z (mm)',fmt='%.2f')
        plot.add_text('P1S / 0.4 mm / ABS | '+summary['status'],font_size=12,color='black')
        plot.view_isometric();plot.reset_camera();plot.camera.zoom(.9);plot.screenshot(out/'QC'/(name+'.png'));plot.close()

def report(out,s):
    a=s['practical_audit'];r=s['removal']
    text=f'''# BG001 RMCA：P1S 0.4 mm / ABS 血管芯与支撑

状态：**{s['status']}**。只生成血管芯与显式支撑，不包含外模或 C 插片。

真实切片使用 **Bambu Lab P1S / 0.4 mm / ABS**，当前保存工艺为 `0.20mm Standard @BBL X1C`。
这个官方共享工艺名称包含 X1C，但其兼容列表明确包含 P1S；机器启动/结束 G-code 使用 P1S 机器配置。
层高 {s['layer_height_mm']} mm、常规及支撑线宽 {s['line_width_mm']} mm，自动支撑关闭。
血管原始 STL 逐字节保留，打印姿态沿用已有 O002。支撑根据新工艺重新生成，未放大或增粗血管。

## 文件

`selected/` 包含原始血管、打印坐标血管 STL、独立支撑 STL、可编辑工程 `_project.3mf`、
原生真实切片 `_sliced.3mf`、`plate_1.gcode` 和坐标变换。
在 Bambu Studio 中选择“打开项目”，打开 `_sliced.3mf` 可看真实刀路，打开 `_project.3mf` 可编辑。
不要把两个 STL 分别自动落板；二者共用同一打印坐标，血管通过显式支撑承托。

`QC/` 保存组合、独立血管、独立支撑图和原生切片预览；`candidates/` 保存每次生成与切片日志。
`regeneration_summary.json` 保存配置、门槛、源保护和交付文件 SHA256。

## 本次实际检查

- 真正缺承托：{a['structural']['true_missing_support_count']} 个；最大缺失弧长 {a['structural']['max_true_unsupported_arc_mm']:.6f} mm。
- 支撑宏观失败：{a['macro']['macro_fail_count']}；浮空血管起始区域：{a['floating_vascular_regions_unresolved']}。
- 实际血管刀路特征保留：{a['vascular_features_preserved']}。
- 几何生成阶段遗留路由错误：{len(s['geometry_routing_errors'])} 个，逐项记录于摘要。
- 薄连接验证：{r['thin_web_preserved']}；完整安全拆除路径：{r['sequence']['complete']}。
- 未通过门槛：{', '.join(s['failures']) or '无'}。

切片成功不代表实物制造验收。只有全部软件门槛通过时才标为可首次试印；未通过时这些文件是诊断候选。
ABS 实物表面、拆卸力与血管损伤尚未实测。本次没有向打印机发送文件或启动打印。
历史 X1 Carbon 输出、外模、C 插片及用户 GUI 配置均未改写。
'''
    (out/'P1S_04_ABS使用说明.md').write_text(text)

def run(output=None,stage='all',max_iterations=3):
    out=Path(output or default_output()).resolve()
    if out==legacy_output().resolve() or out== (legacy_output()/'final_practical_v3').resolve():
        raise ValueError('P1S_REQUIRES_SEPARATE_OUTPUT_DIRECTORY')
    if stage=='review':
        # No freeze/copy/render here: migrating the working directory must not
        # rewrite the archived source paths or any original delivery bytes.
        previous_result=load(out/'regeneration_summary.json')
        verify(load(out/'validation/protected_before.json'))
        verify(load(out/'validation/historical_protection.json'))
        verify({str(out/name):value for name,value in previous_result['files'].items()})
        bundle=profiles(out/'bambu/profiles',review_only=True)
        folder=out/'candidates'/previous_result['selected_candidate']
        evaluated(folder,bundle,readonly=True)
        print('P1S_REVIEW',previous_result['status'],previous_result['failures'],flush=True)
        return previous_result
    out.mkdir(parents=True,exist_ok=True)
    source=discover_source();freeze(source,out)
    previous=legacy_output()/'final_practical_v3'
    history_path=out/'validation/historical_protection.json'
    if not history_path.exists():
        protected={str(p):digest(p) for p in (previous/'selected').iterdir() if p.is_file()}
        protected.update({str(previous/n):digest(previous/n) for n in ['final_practical_summary.json','production_candidate_freeze.json']})
        save(history_path,protected)
    verify(load(history_path))
    bundle=profiles(out/'bambu/profiles')
    # This affects this process only and does not change system/GUI configuration.
    if Path('/run/WSL/2_interop').exists():os.environ['WSL_INTEROP']='/run/WSL/2_interop'
    old=load(previous/'final_practical_summary.json');parent=Path(old['candidate_folder'])
    orientation=load(parent/'geometry.json')['orientation']
    bs,top=branches(source,min(bundle['values']['line_width'],bundle['values']['support_line_width'])/2)
    save(out/'source/production_target.json',dict(printer=PRESETS,profile=bundle['values'],orientation=orientation,anatomy_sha256=digest(source['core'])))
    if stage=='prepare':return dict(status='PREPARED',output=str(out))
    base=out/'candidates/P1S_04_BASE'
    if not (base/'geometry.json').exists():
        if stage=='review':raise FileNotFoundError(base/'geometry.json')
        generate(trimesh.load_mesh(source['core']),bs,orientation,bundle['values'],base,
                 gap=bundle['values']['roof_gap_mm'],lines=2,roof_lines=2,reinforce=True,wide_foot=True)
    if stage=='geometry':return dict(status='GEOMETRY_COMPLETE',output=str(out))
    evaluated(base,bundle)
    parent=base;history=[]
    for i in range(1,max_iterations+1):
        a,structural,_=audit_candidate(parent)
        if structural['true_missing_support_count']==0:break
        folder=out/'candidates'/f'P1S_04_REPAIR_R{i}'
        if stage=='review' and not (folder/'evaluation.json').exists():break
        if not (folder/'geometry.json').exists():repair_true_missing(parent,folder,structural,i)
        evaluated(folder,bundle);new=audit_candidate(folder)[0];accepted=score(new)<=score(a)
        history.append(dict(stage='repair',round=i,candidate=folder.name,accepted=accepted,audit=new))
        if accepted:parent=folder
    # Use the current isolated-neck version of the existing topology method;
    # do not rerun the retired decorative R1/R2 neck variants or orientation search.
    folder=out/'candidates/P1S_04_TOPOLOGY_R3'
    if stage!='review' or (folder/'evaluation.json').exists():
        if not (folder/'geometry.json').exists():topology_candidate(parent,folder,3)
        evaluated(folder,bundle);new=audit_candidate(folder)[0];prior=audit_candidate(parent)[0]
        accepted=score(new)<=score(prior)
        history.append(dict(stage='topology',algorithm_revision=3,candidate=folder.name,accepted=accepted,audit=new))
        if accepted:parent=folder
    save(out/'reports/round_history.json',history)
    return deliver(out,parent,bundle,history)
