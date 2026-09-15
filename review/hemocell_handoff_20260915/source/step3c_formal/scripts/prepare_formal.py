from pathlib import Path
import os,shutil,json,hashlib,difflib
R=Path(__file__).resolve().parents[1];S=R.parent/'20260914_stage4_batch_dispatch';B=(S/'frozen_bundle').resolve()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name,target in [('frozen_bundle',B),('frozen_inputs',B/'frozen_inputs'),('frozen_contracts',R/'contracts'),('upstream',S/'upstream'),('nvhpc',S/'nvhpc')]:
 (R/name).symlink_to(target,target_is_directory=True)
(R/'build_gpu/official/libpalabos.a').symlink_to((S/'build_gpu/official/libpalabos.a').resolve())
for p in (S/'source').iterdir():
 if p.is_file():shutil.copyfile(p,R/'source'/p.name)
for p in (B/'frozen_contracts').iterdir():
 if p.is_file():os.link(p,R/'contracts'/p.name)
for name in ['convergence.py','remote_common.py','verify_run.py','test_remote_logic.py','export_remote_fields.py']:
 shutil.copyfile(R/'reference'/name,R/'scripts'/name)
shutil.copyfile(S/'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json',R/'contracts/GPU_NUMERICAL_COMPARISON_CONTRACT.json')
old=(B/'frozen_contracts/step3_solver_parameters.txt').read_text().splitlines();lines=list(old)
nums=lines[2].split();assert nums[-1]=='1000';nums[-1]='700000';lines[2]=' '.join(nums)
lines[0]='frozen_inputs/step1_geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl';lines[1]='frozen_inputs/step2'
assert lines[2].split()[:-1]==old[2].split()[:-1] and lines[3:]==old[3:]
(R/'contracts/solver_parameters.txt').write_text('\n'.join(lines)+'\n')
(R/'contracts/monitor_parameters.txt').write_text('100 5000 119872 240000 180543\n')
x=(R/'source/vascularPoC.cpp').read_text();before=x
def replace(a,b):
 global x
 assert x.count(a)==1,(a[:100],x.count(a))
 x=x.replace(a,b)
replace('test="GPU_POC_STAGE_4"','test="STEP3C_FORMAL_GPU_LONG_VALIDATION"')
replace('require(requested==200||requested==1000||requested==5000,"Stage4 frozen horizon mismatch");','require(requested==700000,"Formal frozen horizon mismatch");')
replace('const int warmup=requested==200?100:requested/10;','const int warmup=500; // Original Stage4 5000-case timing marker; full-run time also reported.')
replace('maxSteps=requested; // Independent Step3C execution cap; frozen Step3B file stays byte-identical.\n        require(maxSteps<=5000,"Stage4 hard cap 5000");','require(maxSteps==requested,"Formal max_steps must equal original 700000 contract");')
replace('std::set<int> milestones={0,100,200,1000,5000};\n        if(std::getenv("STAGE4_STRONG_CHECKPOINTS")){milestones.insert(1);milestones.insert(10);}','std::set<int> milestones={0,1000,5000,10000,50000,100000,200000,240000,300000,400000,500000,600000};')
replace('bool safe=!s.nan','bool safe=std::isfinite(s.rsum)&&std::isfinite(s.cvsum)&&std::isfinite(s.speedSum)&&s.rsum>0&&s.cvsum>0&&!s.nan')
replace('{Stage2Scope p("OUTPUT_FLOW_HISTORY");history<<step','if(step%scalarEvery==0||!safe) {Stage2Scope p("OUTPUT_FLOW_HISTORY");history<<step')
replace('require(bool(history),"History write failure");}\n                previousSample=step;previousCV=cvMass;for(int p=0;p<4;++p) previousMassFlux[p]=mflux[6*p+2];','require(bool(history),"History write failure");\n                previousSample=step;previousCV=cvMass;for(int p=0;p<4;++p) previousMassFlux[p]=mflux[6*p+2];}')
a='''            if(milestones.count(step)){Stage2Scope p("OUTPUT_SNAPSHOTS");
                snapshot(run,step,nodes);sampledSnapshot(run,step,sampleNodes);
                sparse<<step<<",poc_same_output_snapshot\\n";sparse.flush();

            }'''
b='''            bool eval=step>0&&step%evaluateEvery==0;
            bool witness=step+fieldDelta>=firstEvaluation&&step+fieldDelta<=maxSteps&&(step+fieldDelta)%evaluateEvery==0;
            if(milestones.count(step)||eval||witness||step==maxSteps){
                Stage2Scope p("OUTPUT_SNAPSHOTS");snapshot(run,step,nodes);
                if(milestones.count(step)||eval||step==maxSteps)sampledSnapshot(run,step,sampleNodes);
            }
            if(milestones.count(step)){sparse<<step<<",milestone\\n";sparse.flush();}
            if(eval){
                evaluate(root,run,step); // Synchronous; step 5000 identity gate must return PASS.
                std::ifstream decision(run+"/diagnostics/decision.txt");int evaluated,stop,candidate;
                decision>>evaluated>>stop>>candidate;
                require(bool(decision)&&evaluated==step&&(stop==0||stop==1)&&(candidate==0||candidate==1),"Invalid frozen evaluator decision");
                if(candidate){sparse<<step<<",first_candidate\\n";sparse.flush();}
                if(stop){sparse<<step<<",confirmation\\n";sparse.flush();autoConverged=true;}
            }'''
replace(a,b)
replace('if(step>0&&step%100==0){','if(step>0&&step%1000==0){')
replace("<<step-99<<','<<step<<\",100,\"<<elapsed<<','<<100/elapsed", "<<step-999<<','<<step<<\",1000,\"<<elapsed<<','<<1000/elapsed")
# Both terminal checks should recognize confirmed early stop; no other stopping trigger.
assert x.count('if(step==maxSteps)')==2
x=x.replace('if(step==maxSteps)','if(autoConverged||step==maxSteps)')
replace('<<maxSteps-warmup','<<completed-warmup')
replace('<<(maxSteps-warmup)/timedElapsed','<<(completed-warmup)/timedElapsed')
replace('<<maxSteps<<",\\"full_iteration_seconds\\"','<<completed<<",\\"full_iteration_seconds\\"')
replace('<<maxSteps/fullElapsed','<<completed/fullElapsed')
replace('<<"STAGE3_HORIZON_COMPLETE"','<<(autoConverged?"AUTO_PASS_HUMAN_PENDING":"NOT_CONVERGED_WITHIN_VALIDATION_HORIZON")')
replace('GPU_POC_STAGE3_COMPLETE','STEP3C_FORMAL_GPU_COMPLETE')
(R/'source/vascularPoC.cpp').write_text(x)
(R/'provenance/FORMAL_EXECUTION_HARNESS.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),x.splitlines(True),fromfile='sealed_stage4/vascularPoC.cpp',tofile='formal/source/vascularPoC.cpp')))
proof={}
for name in ['gpu_adapter.hpp','stage4_lbm_batch.hpp','stage2_profile.hpp','mpi_support.hpp','CMakeLists.txt']:
 assert sha(S/'source'/name)==sha(R/'source'/name);proof[name]=sha(R/'source'/name)
(R/'provenance/STAGE4_COMPUTE_IDENTITY.json').write_text(json.dumps(dict(status='PASS',unchanged_files=proof,changes='Execution horizon/cadence, original frozen evaluation and snapshot schedule, synchronous 5000 identity, finite scalar safety assertion only. No compute/BC/halo/monitor reduction/layout changes.',convergence_algebra_sha256=sha(R/'scripts/convergence.py')),indent=2)+'\n')
print('GENERATED',R)
