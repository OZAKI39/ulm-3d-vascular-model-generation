"""Bounded final practical closure. No new orientation or architecture search."""
from pathlib import Path
from datetime import datetime,timezone
import json,pickle,shutil,zipfile
from .core_support_inputs import default_output,discover_source,freeze,current_profiles
from .core_support_practical_audit import audit_candidate
from .core_support_practical_geometry import repair_true_missing
from .core_support_practical_topology import topology_candidate
from .core_support_practical_removal import audit_removal,practical_gates,REMOVAL_AUDIT_VERSION
from .core_support_slicing import evaluate
from .support_inputs import save,table,digest,verify,ROOT

READY='READY_FOR_REAL_PRINT_TRIAL_WITH_PHYSICAL_REMOVAL_UNCALIBRATED'


def load(path):return json.loads(Path(path).read_text())


def baseline_selection(base,out):
    # Read the already ranked top alternatives. Do not launch a new search.
    names=['O002_R25_C2_LAYER_BASE','O000_R20_R3_C3_LAYER','O001_R25_C2'];rows=[]
    for name in names:
        folder=base/'candidates'/name
        if not (folder/'evaluation.json').exists():raise FileNotFoundError(folder/'evaluation.json')
        a,_,_=audit_candidate(folder,out/'baseline_audits'/name)
        oldrem=load(folder/'removal.json') if (folder/'removal.json').exists() else None
        rows.append(dict(candidate=name,folder=str(folder),true_missing=a['structural']['true_missing_support_count'],macro_fail=a['macro']['macro_fail_count'],removal_topology=oldrem,practical=a))
    selected=rows[0]
    for alternative in rows[1:]:
        # Unknown removal evidence cannot be called strictly better. Historical
        # gates on 120 degrees are not used as the new practical acceptance.
        old=selected['removal_topology'];new=alternative['removal_topology']
        better=old and new and new['removable_segment_fraction']>old['removable_segment_fraction'] and new['maximum_wrap_angle_deg']<old['maximum_wrap_angle_deg']
        if alternative['true_missing']<selected['true_missing'] and alternative['macro_fail']<selected['macro_fail'] and better:selected=alternative
    save(out/'reports/baseline_selection.json',dict(selected=selected['candidate'],candidates=rows,reason='Prefer current O002 baseline unless an existing alternative strictly dominates all three measured criteria; no new search'))
    return Path(selected['folder'])


def cached_evaluate(folder,actual,engine,out,stage,round_number):
    """A completed native archive is never sliced again just for a report."""
    folder=Path(folder);result=folder/'evaluation.json'
    if result.exists():return load(result)
    saved=folder/'bambu/slice_result.json'
    if saved.exists() and load(saved).get('slicer_success'):
        return evaluate(load(folder/'geometry.json'),actual,engine,reuse_native=True)
    ledger=out/'validation/native_attempts.json';rows=load(ledger) if ledger.exists() else []
    prior=[r for r in rows if r['stage']==stage and r['round']==round_number]
    if prior:raise RuntimeError('RECOVER_EXISTING_NATIVE_ATTEMPT_WITHOUT_RESLICING: '+str(folder))
    if round_number not in [1,2,3]:raise ValueError('THREE_ROUNDS_MAXIMUM')
    rows.append(dict(stage=stage,round=round_number,folder=str(folder),started=datetime.now(timezone.utc).isoformat()));save(ledger,rows)
    return evaluate(load(folder/'geometry.json'),actual,engine)


def finalise(base,out,folder,profiles,history):
    audit,structural,macro=audit_candidate(folder)
    removal=load(folder/'practical_removal/removal_audit.json')
    e=load(folder/'evaluation.json');protected=verify(load(out/'validation/protected_before.json'));save(out/'validation/source_protection.json',protected)
    with zipfile.ZipFile(e['slice']['archive']) as z:archive_valid=z.testzip() is None
    archive_valid=archive_valid and e['archive_geometry']['passed'] and e['independent_part_count']==2
    gates=practical_gates(audit,removal,protected['all_unchanged'],archive_valid);status=READY if all(gates.values()) else 'NEEDS_ADJUSTMENT'
    selected=out/'selected';selected.mkdir(exist_ok=True)
    # Byte-identical source and the rigid print-frame copy are both delivered;
    # the print-frame copy cannot truthfully have the untransformed file hash.
    files={'frozen_vascular_core.stl':out/'source/frozen_vascular_core.stl','vascular_core_print_frame.stl':folder/'vascular_core.stl','custom_support.stl':folder/'custom_support.stl','final_project.3mf':folder/'bambu/input.3mf','final_sliced.3mf':Path(e['slice']['archive']),'vascular_print_transform.json':folder/'vascular_print_transform.json','inverse_transform.json':folder/'inverse_transform.json'}
    if not files['final_project.3mf'].exists():
        projects=[p for p in (folder/'bambu').glob('*.3mf') if p!=Path(e['slice']['archive'])]
        if len(projects)!=1:raise ValueError('AMBIGUOUS_NATIVE_INPUT_PROJECT')
        files['final_project.3mf']=projects[0]
    for name,source in files.items():shutil.copy2(source,selected/name)
    with zipfile.ZipFile(selected/'final_sliced.3mf') as z:
        for name in z.namelist():
            if name.endswith('.gcode'):
                target=selected/Path(name).name
                with z.open(name) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
    for name in ['support_failure_reclassification.csv','support_self_macro_audit.csv']:
        (out/'tables').mkdir(exist_ok=True);shutil.copy2(folder/'practical'/name,out/'tables'/name);shutil.copy2(folder/'practical'/name,out/name)
    for name in ['final_true_missing_support.json','support_self_macro_audit.json','floating_vascular_regions.json']:shutil.copy2(folder/'practical'/name,out/name)
    for name in ['strategic_cut_seams.csv','major_removal_components.csv','final_contact_pieces.csv','sequential_access.json','safe_removal_sequence.json','bulk_cut_attempts.json']:
        shutil.copy2(folder/'practical_removal'/name,out/name)
        if name.endswith('.csv'):shutil.copy2(folder/'practical_removal'/name,out/'tables'/name)
    save(out/'reports/round_history.json',history)
    inventory=[]
    for candidate in sorted((out/'candidates').iterdir()):
        native=candidate/'bambu/slice_result.json'
        if not native.exists():continue
        result=load(native);archive=result.get('archive')
        inventory.append(dict(candidate=candidate.name,stage='repair' if 'REPAIR' in candidate.name else 'topology',round=int(candidate.name.rsplit('R',1)[1]),slicer_success=result.get('slicer_success'),command_file=str(candidate/'bambu/command.json'),archive=archive,archive_sha256=digest(archive) if archive and Path(archive).exists() else None))
    save(out/'validation/native_run_inventory.json',inventory)
    code=[ROOT/'s1-11_vascular_core_support.py',*sorted((ROOT/'vascular_processing').glob('core_support_*.py'))]
    freeze_manifest=dict(status=status,production_frozen=status==READY,date=datetime.now(timezone.utc).isoformat(),candidate=folder.name,vascular_source_sha256=digest(selected/'frozen_vascular_core.stl'),placed_vascular_sha256=digest(selected/'vascular_core_print_frame.stl'),support_sha256=digest(selected/'custom_support.stl'),project_3mf_sha256=digest(selected/'final_project.3mf'),sliced_3mf_sha256=digest(selected/'final_sliced.3mf'),profile_sha256=profiles['hashes'],bambu_version=profiles['engine']['version'],print_transform=load(selected/'vascular_print_transform.json'),source_protection=protected,files={str(p.relative_to(out)):digest(p) for p in selected.iterdir() if p.is_file()},code_sha256={str(p.relative_to(ROOT)):digest(p) for p in code},physical_removal_force_validated=False,physical_vessel_damage_validated=False,physical_surface_quality_validated=False)
    save(out/'production_candidate_freeze.json',freeze_manifest)
    freeze_manifest['production_candidate']=folder.name if status==READY else None
    freeze_manifest['diagnostic_candidate']=None if status==READY else folder.name
    save(out/'production_candidate_freeze.json',freeze_manifest)
    summary=dict(status=status,selected_candidate=folder.name,candidate_folder=str(folder),output=str(out),gates=gates,failures=[k for k,v in gates.items() if not v],practical_audit=audit,structural={k:v for k,v in structural.items() if k!='rows'},macro={k:v for k,v in macro.items() if k!='rows'},removal={k:v for k,v in removal.items() if k not in ['major_components','final_contact_pieces','cuts','thin_web_checks']},bambu=dict(version=profiles['engine']['version'],slicer_success=e['slice']['slicer_success'],auto_support_off=e['auto_support_off'],profile_hashes=profiles['hashes'],archive=str(selected/'final_sliced.3mf'),metrics=e['slice']),legacy_metrics=dict(coverage=e['coverage'],layer=e['layer_coverage'],support_self=e['support_self']),source_protection=protected,round_history=history,physical_uncertainty=['REAL ABS SUPPORT REMOVAL FORCE NOT VALIDATED','REAL VESSEL DAMAGE RISK NOT VALIDATED','SURFACE QUALITY NOT VALIDATED'],production_frozen=status==READY,repair_round_count=len(list((out/'candidates').glob('PRACTICAL_REPAIR_R*/evaluation.json'))),topology_round_count=len(list((out/'candidates').glob('PRACTICAL_TOPOLOGY_R*/evaluation.json'))))
    save(out/'final_practical_summary.json',summary)
    summary['actual_native_slice_attempt_count']=len(inventory)
    summary['stop_reason']='READY_FREEZE' if status==READY else 'BOUNDED_FINAL_ROUNDS_COMPLETED_WITH_UNRESOLVED_GATES'
    save(out/'final_practical_summary.json',summary)
    from .core_support_practical_review import render,report
    render(base,out,folder,summary,structural,macro,removal);report(out,summary,structural,removal)
    return summary


def run(base=None,max_rounds=3,review_only=False):
    if max_rounds not in [1,2,3]:raise ValueError('THREE_REPAIR_AND_THREE_TOPOLOGY_ROUNDS_MAXIMUM')
    base=Path(base or default_output()).resolve();out=base/'final_practical_v3';out.mkdir(parents=True,exist_ok=True)
    frozen=out/'production_candidate_freeze.json'
    if frozen.exists() and load(frozen).get('production_frozen'):
        manifest=load(frozen);verify({str(out/p):h for p,h in manifest['files'].items()});return load(out/'final_practical_summary.json')
    source=discover_source();freeze(source,out)
    pp=out/'bambu/profiles/profiles.json'
    if not pp.exists():current_profiles(source,pp.parent)
    profiles=load(pp);history=[]
    verify({r['path']:r['sha256'] for r in profiles['values']['profile_sources']})
    parent=baseline_selection(base,out)
    for i in range(1,max_rounds+1):
        a,structural,_=audit_candidate(parent)
        if structural['true_missing_support_count']==0:break
        folder=out/'candidates'/f'PRACTICAL_REPAIR_R{i}'
        if not (folder/'evaluation.json').exists():
            if review_only:break
            if not (folder/'geometry.json').exists():repair_true_missing(parent,folder,structural,i)
            cached_evaluate(folder,profiles['actual'],profiles['engine'],out,'repair',i)
        new,_,_=audit_candidate(folder)
        accepted=new['structural']['true_missing_support_count']<=a['structural']['true_missing_support_count'] and new['macro']['macro_fail_count']<=a['macro']['macro_fail_count']
        history.append(dict(stage='repair',round=i,candidate=folder.name,accepted=accepted,practical=new))
        if accepted:parent=folder
    baseline_audit=audit_candidate(parent)[0];options=[]
    for i in range(1,max_rounds+1):
        folder=out/'candidates'/f'PRACTICAL_TOPOLOGY_R{i}'
        if not (folder/'evaluation.json').exists():
            if review_only:break
            if not (folder/'geometry.json').exists():topology_candidate(parent,folder,i)
            cached_evaluate(folder,profiles['actual'],profiles['engine'],out,'topology',i)
        a,_,_=audit_candidate(folder);rp=folder/'practical_removal/removal_audit.json';rem=load(rp) if rp.exists() else {}
        if rem.get('audit_version')!=REMOVAL_AUDIT_VERSION:rem=audit_removal(folder)
        no_regression=a['structural']['true_missing_support_count']<=baseline_audit['structural']['true_missing_support_count'] and a['macro']['macro_fail_count']<=baseline_audit['macro']['macro_fail_count']
        gates=practical_gates(a,rem)
        history.append(dict(stage='topology',round=i,candidate=folder.name,structural_no_regression=no_regression,practical=a,max_contact_arc_mm=rem['maximum_component_contact_arc_mm'],max_effective_lever_arm_mm=rem['max_effective_lever_arm_mm'],major_components=rem['major_component_count'],sequence_complete=rem['sequence']['complete'],c3_preserved=rem['thin_web_preserved']))
        options.append((not no_regression,a['structural']['true_missing_support_count'],a['macro']['macro_fail_count'],not rem['thin_web_preserved'],not rem['sequence']['complete'],rem['maximum_component_contact_arc_mm'],rem['max_effective_lever_arm_mm'],folder))
        if all(gates.values()):return finalise(base,out,folder,profiles,history)
    if not options:raise RuntimeError('NO_COMPLETED_TOPOLOGY_CANDIDATE_TO_REVIEW')
    chosen=min(options,key=lambda row:row[:-1])[-1]
    return finalise(base,out,chosen,profiles,history)
