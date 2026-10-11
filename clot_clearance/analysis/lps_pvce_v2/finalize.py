"""Read-only final checks, exact source diffs and task-only manifest; no solves."""
from analysis.lps_pvce_v2.common import *
import subprocess,difflib,os,re

def main():
 p=BASE/'provenance';s=json.loads((BASE/'FINAL_STATUS.json').read_text());assert s['DOCUMENT_UPDATE']=='ACTUAL_WSL_UPDATED_AND_READ_BACK';prot=json.loads((p/'PROTECTED_START.json').read_text());changed=[n for n,h in prot.items() if sha(ROOT/n)!=h];assert not changed;save(p/'PROTECTED_FINAL.json',dict(checked=len(prot),changed=changed))
 source=(ROOT/'tools/lammps-contact-onset-v1/src/BPM/bond_bpm_peri.cpp').read_text();assert 'const bool allow_breaks = (update->setupflag == 0) && break_flag' in source;assert 'stretch > s00[type] - alpha[type] * MAX(smin[i1], smin[i2])' in source
 save(p/'NATIVE_CRITERION_AUDIT.json',dict(source='tools/lammps-contact-onset-v1/src/BPM/bond_bpm_peri.cpp',sha256=sha(ROOT/'tools/lammps-contact-onset-v1/src/BPM/bond_bpm_peri.cpp'),exact_criterion='allow_breaks && (stretch > s00[type] - alpha[type] * MAX(smin[i1], smin[i2]))',allow_breaks='update->setupflag == 0 && break_flag',alpha=0,s00=.15,all_reference_bonds_type=1,scope='Native forces checked with task-only break no; offline threshold diagnosis is exact for current alpha0 bulk-only states, not a dynamic event.'))
 nativeold=(ROOT/'analysis/lps_pvce_v1/core.py').read_text().split('class Native:')[1];nativenew=(ROOT/'analysis/lps_pvce_v2/native.py').read_text().split('class Native:')[1];assert nativeold.replace('bond_style bpm/peri\n','bond_style bpm/peri break no\n')==nativenew
 (p/'NATIVE_ADAPTER.diff').write_text(''.join(difflib.unified_diff(nativeold.splitlines(True),nativenew.splitlines(True),fromfile='V1 Native',tofile='V2 Native (task only)')))
 diffs=[]
 for fn in ['fem.py','campaign.py','lps.py']:
  before=(p/('launch_'+fn)).read_text();after=(ROOT/'analysis/lps_pvce_v2'/fn).read_text();diffs.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='launch/'+fn,tofile='final/'+fn))
 (p/'POST_LAUNCH_SOURCE.diff').write_text(''.join(diffs));new=[]
 for f in sorted((ROOT/'analysis/lps_pvce_v2').glob('*.py')):new.extend(difflib.unified_diff([],f.read_text().splitlines(True),fromfile='/dev/null',tofile=str(f.relative_to(ROOT))))
 (p/'NEW_TASK_SOURCE.diff').write_text(''.join(new));save(p/'FINAL_SOURCE_IDENTITY.json',{str(f.relative_to(ROOT)):sha(f) for f in sorted((ROOT/'analysis/lps_pvce_v2').glob('*.py'))})
 checks=json.loads((BASE/'IMPLEMENTATION_CHECKS.json').read_text());checks['native_source_criterion']='allow_breaks && stretch > s00[type] - alpha[type]*MAX(smin_i,smin_j); alpha0 bulk-only gives stretch>.15';save(BASE/'IMPLEMENTATION_CHECKS.json',checks)
 # Capture the isolated installed package identity without copying the environment.
 pkg=Path('/workspace/y_clot/envs/pvce_v2_packages');save(p/'ISOLATED_PACKAGE_FILES.json',{str(f.relative_to(pkg)):sha(f) for f in sorted(pkg.rglob('*')) if f.is_file() and '__pycache__' not in f.parts})
 required=['REVIEW_LPS_PVCE_CONTINUUM_AND_NEAR_ONSET_AUDIT_V2.md','FINAL_STATUS.json','V2_DESIGN_FREEZE.json','CONTINUUM_REFERENCE_AUDIT.json','CONTINUUM_PVCE_COMPARISON.csv','A0_FIRST_FAILURE_SEARCH.json','FIRST_FAILURE_LOCATION_AUDIT.json','ENERGY_AND_EQUILIBRIUM_AUDIT.json','IMPLEMENTATION_CHECKS.json','OPEN_RESULTS.html'];assert all((BASE/f).exists() for f in required);assert not (BASE/'NEAR_ONSET_PVCE_RESULTS.json').exists();assert len(list(FIG.glob('Fig*.png')))==4 and len(list(FIG.glob('Fig*.svg')))==4
 for idx in [BASE/'OPEN_RESULTS.html',FIG/'OPEN_RESULTS.html']:
  for href in re.findall(r'(?:src|href)="([^"]+)"',idx.read_text()):assert (idx.parent/href).exists(),href
 statuses=subprocess.run(['supervisorctl','status','y_clot_pvce_v2_continuum_'+RUN,'y_clot_pvce_v2_onset_'+RUN],capture_output=True,text=True).stdout;assert statuses.count('EXITED')==2
 (BASE/'REPLAY.md').write_text('''# 最小复核与再运行限制

现有服务器：`source /workspace/y_clot/env.sh`，进入`/workspace/y_clot/project`，运行`python -m analysis.lps_pvce_v2.replay --verify-only`。此命令只核验清单和储能差，绝不重新求解。

WSL轻量复核只需Python标准库，无须安装NumPy或加载原生求解库。完整原生重算另需锁定LAMMPS库与energy.so，版本身份在保护清单，不能用其他库悄悄替代。FEM额外安装scikit-fem12.0.2到独立target目录（环境记录中给出），不得升级原环境。

实际数值启动命令与Supervisor配置在provenance/*_launch.sh、*_supervisor.conf；六个FEM ATTEMPT和一个LPS ATTEMPT对应七个独立求解，没有重复计算。旧run目录不可用于覆盖重放。新数值试验需新的路径、冻结方案及授权范围，本轮到停止条件后不自动追加。

三份重复body.data完整文本留服务器；它们是同一个完整未切节点/参考键表，可由旧V1 reference.npz和锁定write_typed_data重建。WSL已保存必要全精度网格、P2/P1自由度、节点状态、力、alive、m及全部数值审计。body.data并非唯一科学数据源。

FILE_MANIFEST.json不收录自己、集合收据及核验回执，以免循环哈希；WSL缺失仅允许上述body.data。最终清单核验同时在服务器和WSL执行。没有迁移其他项目、旧大轨迹、软件环境或视频。
''')
 save(p/'COMPLETION_CHECK.json',dict(required_artifacts=required,high_PVCE_results_absent_by_admission=True,actual_FEM_solves=6,actual_LPS_candidates=1,new_LPS_accepted=0,source_locks_unchanged=True,protected_files=len(prot),supervisor=statuses,document_status=s['DOCUMENT_UPDATE'],visual_review=json.loads((p/'VISUAL_REVIEW.json').read_text())))
 exclude={'FILE_MANIFEST.json','MANIFEST_VERIFICATION.json','WSL_FINAL_VERIFICATION.json','WSL_COMPACT_COLLECTION.json'};files=[f for folder in [BASE,FIG] for f in folder.rglob('*') if f.is_file() and f.name not in exclude];files+=list((ROOT/'analysis/lps_pvce_v2').glob('*.py'));records=[dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=sha(f)) for f in sorted(files)];save(BASE/'FILE_MANIFEST.json',dict(files=records,total_bytes=sum(r['bytes'] for r in records),excluded=list(sorted(exclude)),scientific_scope='Task-scoped source, full precision states, meshes, accepted/failed diagnostics, figures and actual input logs. Production and old results remain read-only.'))
 print('FINALIZED',len(records),'files',sum(r['bytes'] for r in records),'bytes')
if __name__=='__main__':main()
