"""Explicit, hash-checked cleanup authorized by BCF steer; no source deletion."""
from pathlib import Path
import ast,collections,hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from vascular_processing.support_inputs import discover,digest,save,verify
OUT=ROOT/'research_retrospective';BASE=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner'

def main(execute=False):
    OUT.mkdir(exist_ok=True)
    commands={c:subprocess.run(c,cwd=ROOT,shell=True,capture_output=True,text=True).__dict__ for c in ['git status --short','git rev-parse HEAD','git ls-files']}
    save(OUT/'git_audit.json',commands)
    inputs=discover();refs=collections.defaultdict(list);files=[]
    for directory,dirs,names in os.walk(ROOT,followlinks=False):
        dirs[:]=[d for d in dirs if not (Path(directory)/d).is_symlink() and d not in ('.git','.venv')]
        for n in names:
            p=Path(directory)/n
            if not p.is_symlink():files.append(p)
    for p in files:
        if p.suffix!='.py' or 'third_party' in p.parts:continue
        try:t=ast.parse(p.read_text())
        except (SyntaxError,UnicodeError):continue
        for n in ast.walk(t):
            names=([n.module] if isinstance(n,ast.ImportFrom) and n.module else [a.name for a in n.names] if isinstance(n,ast.Import) else [])
            for name in names:refs[name.split('.')[-1]].append(str(p.relative_to(ROOT)))
    save(OUT/'import_reference_scan.json',refs)
    keepnames={'vascular_support_report.md','vascular_support_summary.json','vascular_support_v2_report.md','vascular_support_v2_summary.json','vascular_support_v3_report.md','vascular_support_v3_summary.json','orientation_then_v3_1_report.md','support_orientation_and_closure_summary.json','closure_summary.json','inverted_study_report.md','robust_full_underbody_report.md','robust_full_support_summary.json','method_comparison.csv','baseline_comparison.csv','v2_vs_v3.csv','protected_before.json','source_protection.json','source_snapshot_manifest.json','pytest.xml','pytest.log','test_results.xml','test_run.log','artifact_verification.json','runtime_environment.json'}
    inventory=[];protected={};corehash=digest(inputs['core'])
    for p in files:
        rel=p.relative_to(ROOT);parts=rel.parts;size=p.stat().st_size;h=digest(p)
        cat='REVIEW_REQUIRED';reason='Outside explicitly retired outputs; preserve without inferring obsolescence.'
        if '__pycache__' in parts or '.pytest_cache' in parts:cat='DELETE_REGENERABLE';reason='Regenerable Python/pytest cache.'
        elif p.is_relative_to(OUT):cat='KEEP_COMPACT_HISTORY';reason='Current cleanup audit.'
        elif p.is_relative_to(BASE):
            sub=p.relative_to(BASE);history=any(k=='history' or k.startswith('iteration_') or k=='candidates' or k=='custom_candidates' for k in sub.parts)
            if not history and (p.name in keepnames or 'profiles' in sub.parts or 'source_snapshot' in sub.parts):
                cat='KEEP_COMPACT_HISTORY';reason='Final report, comparison, tests, source hashes or reproducibility source/profile snapshot.'
            else:cat='DELETE_REGENERABLE';reason='Retired support candidate, intermediate, debug or visual output; final evidence retained.'
        elif p.suffix in ('.py','.sh','.yaml','.yml','.toml') or parts[0] in ('tools','third_party','tests','config','configs','utils','vascular_processing'):
            cat='KEEP_SHARED_CODE';reason='Source/tests/shared helpers retained: no recoverable project Git history; dependency scan recorded.'
        elif parts[0] in ('data','vessel_model','outputs','deliverables','exports','test_data','references'):
            cat='KEEP_PRODUCTION';reason='Outside retired support path; frozen anatomy, production or provenance preserved.'
        if str(p) in inputs['hashes']:cat='KEEP_PRODUCTION';reason='Discovered frozen production input.'
        row=dict(path=str(rel),size_bytes=size,tracked=False,git_state='UNTRACKED_NO_GIT_REPOSITORY',category=cat,referenced_by=refs.get(p.stem,[]),sha256=h,proposed_action='DELETE' if cat.startswith('DELETE_') else 'KEEP',reason=reason,identical_to_production_core=h==corehash)
        inventory.append(row)
        if cat in ('KEEP_PRODUCTION','KEEP_SHARED_CODE') and not p.is_relative_to(OUT):protected[str(p)]=h
    save(ROOT/'cleanup_inventory.json',inventory);save(OUT/'production_protected_before.json',protected)
    deletes=[r for r in inventory if r['proposed_action']=='DELETE'];counts=collections.Counter(r['category'] for r in inventory)
    plan='# 支撑研发输出清理计划（先生成、后执行）\n\n本目录没有 Git 历史，不删除任何算法源码或测试。s1-8 系列降级为研究历史。只删除明确退休支撑目录内可重新生成的输出和 Python/pytest 缓存。所有其他生产、原始数据和共享源码保留；无法确定用途的文件标为 REVIEW_REQUIRED 并保留。\n\n'
    plan+=f'项目文件数：{len(files)}；逻辑字节：{sum(r["size_bytes"] for r in inventory)}。预计删除 {len(deletes)} 文件，释放 {sum(r["size_bytes"] for r in deletes)} 字节。\n\n'
    plan+='\n'.join(f'- {k}: {v}' for k,v in counts.items())+'\n\n逐文件 KEEP/DELETE/REVIEW、导入引用和 SHA256 见 cleanup_inventory.json；生产保护见 research_retrospective/production_protected_before.json。相同 SHA 的候选副本只在明确退休目录删除，权威生产路径永久保留。旧报告内大文件链接属于历史证据，清理后不代表现存交付物。\n'
    (ROOT/'cleanup_plan.md').write_text(plan)
    if not execute:print(plan);return
    if not (OUT/'SUPPORT_RND_RETROSPECTIVE.md').exists():raise RuntimeError('Retrospective required before deletion')
    for r in deletes:
        p=ROOT/r['path']
        if digest(p)!=r['sha256']:raise RuntimeError('Changed since inventory: '+str(p))
        if str(p) in protected:raise RuntimeError('Protected deletion')
        p.unlink()
    for directory,dirs,names in os.walk(BASE,topdown=False):
        p=Path(directory)
        if not any(p.iterdir()):p.rmdir()
    check=verify(protected)
    result=dict(deleted_file_count=len(deletes),deleted_bytes=sum(r['size_bytes'] for r in deletes),kept_file_count=len(inventory)-len(deletes),protected_hash_mismatches=len(check['changed']),protected_files=check['checked_files'],before_file_count=len(inventory),before_bytes=sum(r['size_bytes'] for r in inventory),after_inventory_file_count=len(inventory)-len(deletes),after_inventory_bytes=sum(r['size_bytes'] for r in inventory if r['proposed_action']=='KEEP'),source_core_sha256=corehash,git_commit=None,dead_code_deleted=0,status='CLEANUP_COMPLETED_PENDING_TESTS')
    save(ROOT/'cleanup_result.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main('--execute' in sys.argv)
