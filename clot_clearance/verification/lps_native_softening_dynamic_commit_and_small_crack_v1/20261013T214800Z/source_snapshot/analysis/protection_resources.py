from common import *
import subprocess,time,sys
start=time.perf_counter();protected=json.loads((O/'HISTORICAL_PROTECTION_BEFORE.json').read_text());conflicts=[]
for n,v in protected.items():
 p=R/n
 if not p.is_file() or sha(p)!=v['sha256']:conflicts.append(n)
identity=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());asset_conflicts=[]
for n,v in identity['files'].items():
 p=Path(n)
 if not p.is_file() or sha(p)!=v['sha256']:asset_conflicts.append(n)
for n,h in identity['tracked_modified_SHA'].items():
 if sha(R/n)!=h:asset_conflicts.append(n)
gitdiff=subprocess.check_output(['git','diff','--binary'],cwd=R,text=True);(O/'GIT_DIFF_final.patch').write_text(gitdiff);(O/'GIT_STATUS_final.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=R,text=True));same_diff=gitdiff==(O/'GIT_DIFF_before.patch').read_text();assert not conflicts and not asset_conflicts and same_diff
sources={p.name:sha(p) for p in (A).glob('*') if p.suffix in ['.cpp','.h']}
for n,h in sources.items():assert sha(O/'source/src/BPM'/n)==h
newidentity=dict(library=str(LIB),library_sha256=sha(LIB),source_files=sources,patch_sha256=sha(O/'NATIVE_SOURCE_PATCH.diff'),package_metadata=sha(O/'source/src/.gitignore'),production_sha256=sha(R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so'),previous_accepted_library_sha256=sha(PREV/'build/liblammps.so'),binary_used_for_all_formal_steps=True,CPP_not_changed_after_dynamics=True)
dump(O/'FINAL_LIBRARY_IDENTITY.json',newidentity)
dump(O/'HISTORICAL_PROTECTION_AFTER.json',dict(status='PASS',historical_files_checked=len(protected),historical_conflicts=conflicts,protected_asset_files_checked=len(identity['files']),asset_conflicts=asset_conflicts,tracked_diff_identical=same_diff,git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),wall_s=time.perf_counter()-start))
first=json.loads((O/'build_attempt_1/BUILD_RESOURCES.json').read_text());last=json.loads((O/'BUILD_RESOURCES.json').read_text());main=json.loads((O/'DYNAMIC_341_NODE_AUDIT.json').read_text());half=json.loads((O/'dt_half_window_AUDIT.json').read_text());p0=json.loads((O/'P0_COMPLETE.json').read_text());setup=json.loads((O/'setup_attempt_min_syntax/DYNAMIC_341_NODE_AUDIT.json').read_text());p2=json.loads((O/'P2_FIXED_STATE_REGRESSION.json').read_text());fig=json.loads((O/'FIGURE_RESOURCES.json').read_text())
r=dict(P0_wall_s=p0['wall_s'],deleted_code_files=0,compile_jobs=12,build_attempts=2,build_wall_s=first['wall_s']+last['wall_s'],build_stage_records=[first,last],compile_peak_RSS_bytes=max(v['sampled_peak_process_tree_RSS_bytes'] for k in [first,last] for v in k['stages']),compile_peak_sampled_CPU_cores=max(v['sampled_peak_CPU_cores'] for k in [first,last] for v in k['stages']),dynamic_main_cases=1,P4_cases=0,short_dt_half_windows=1,MPI_ranks=1,threads=1,main_wall_s=main['wall_s'],half_window_wall_s=half['wall_s'],formal_dynamic_wall_s=main['wall_s']+half['wall_s'],formal_dynamic_CPU_s=main['CPU_s']+half['CPU_s'],main_physical_time_s=main['simulated_time_s'],main_accepted_steps=244,half_window_physical_time_s=half['simulated_time_s'],half_window_accepted_steps=68,main_wall_per_step_s=main['wall_s']/244,peak_dynamic_RSS_bytes=max(main['peak_RSS_bytes'],half['peak_RSS_bytes']),zero_step_failed_setup_wall_s=setup['wall_s'],P2_fixed_state_wall_s=p2['wall_s'],figure_resources=fig,cgroup_identity=identity['cgroup'],CPU_scope='Actual Linux process CPU and tree-sampled build CPU. GPU not used.',timing_scope='Formal case wall includes initialization, observation, restart writes and FP64 compression. Zero-step syntax failure and plotting are separate. No claim for uninstrumented orchestration/editor elapsed time.',postprocess=json.loads((O/'POSTPROCESS_RESOURCES.json').read_text()))
dump(O/'RUN_RESOURCES.json',r)
print('PROTECTED',len(protected),'ASSETS',len(identity['files']),'BUILD_SECONDS',r['build_wall_s'],'DYNAMICS_SECONDS',r['formal_dynamic_wall_s'],'RAM_MIB',r['peak_dynamic_RSS_bytes']/1024**2,flush=True)
