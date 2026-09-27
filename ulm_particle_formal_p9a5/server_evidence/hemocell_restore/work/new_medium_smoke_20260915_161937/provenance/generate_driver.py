from pathlib import Path
import hashlib,json,difflib
C=Path(__file__).resolve().parent
P=Path('/home/lzy/projects/server_migrations/hemocell_20260915_094331/source_snapshot/hemocell_gpu_poc/20260914_stage4_batch_dispatch/source/vascularPoC.cpp')
old=P.read_text();assert hashlib.sha256(P.read_bytes()).hexdigest()=='e23363ce5cedd9278e39bedec749ae90a38ff595a263d16e8e2e7a7dbb97e508'
s=old
helper=r'''
// Host input/identity adapter only: consume generated values; never regenerate dt or boundary density.
std::string contractText(std::string const& path) {
    std::ifstream f(path);require(bool(f),"Missing frozen contract: "+path);
    std::ostringstream out;out<<f.rdbuf();require(!out.str().empty(),"Empty frozen contract");return out.str();
}
T contractNumber(std::string const& body,std::string const& key) {
    std::string token="\""+key+"\"";size_t pos=0;
    while((pos=body.find(token,pos))!=std::string::npos){
        size_t end=pos+token.size();while(end<body.size()&&std::isspace(static_cast<unsigned char>(body[end])))++end;
        if(end<body.size()&&body[end]==':'){
            size_t used=0;T v=std::stod(body.substr(end+1),&used);require(used>0&&std::isfinite(v),"Nonfinite contract value: "+key);return v;
        }
        pos+=token.size();
    }
    throw std::runtime_error("Missing numerical contract key: "+key);
}
'''
s=s.replace('#include <chrono>','#include <chrono>\n#include <cctype>')
s=s.replace('int main(int argc,char** argv)',helper+'\nint main(int argc,char** argv)')
s=s.replace('test="GPU_POC_STAGE_4"','test="PURE_FLUID_NEW_MEDIUM_SMOKE"')
s=s.replace('require(requested==200||requested==1000||requested==5000,"Stage4 frozen horizon mismatch");','require(requested==500||requested==5000,"New-medium smoke horizon must be500 or5000");')
s=s.replace('const int warmup=requested==200?100:requested/10;','const int warmup=requested/10;')
s=s.replace('require(maxSteps==700000&&rampSteps==10&&tau==1.,"Frozen execution contract mismatch");','require(maxSteps==5000&&rampSteps==10&&tau==1.,"Smoke execution contract mismatch");')
s=s.replace('maxSteps=requested; // Independent Step3C execution cap; frozen Step3B file stays byte-identical.','maxSteps=requested; // New independent smoke cap; original Stage4/Step3C contracts remain unchanged.')
a=s.index('        int scalarEvery,evaluateEvery,fieldDelta,firstEvaluation,cvExpected;');b=s.index('        std::vector<Port> ports(4);',a)
s=s[:a]+r'''        const int scalarEvery=100,cvExpected=180543;
        std::string numericalText=contractText(run+"/contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json");
        std::string safetyText=contractText(run+"/contracts/SMOKE_SAFETY_CONTRACT.json");
        require(dt==contractNumber(numericalText,"dt_s")&&expectedDx==contractNumber(numericalText,"dx_m")&&tau==contractNumber(numericalText,"tau")&&rhoPhys==contractNumber(numericalText,"rho_kg_m3"),"Actual input differs from generated numerics contract");
        require(qTarget==contractNumber(numericalText,"physical_Qtarget_m3_s")&&multiplier==contractNumber(numericalText,"inlet_numerical_multiplier"),"Physical inlet mismatch");
        require(rhoPhys==1000.&&contractNumber(numericalText,"nu_m2_s")==1e-6&&tau==1.&&dt>0,"Provisional medium input mismatch");
        require(rhoLow==.99&&rhoHigh==1.01&&maLimit==.05,"Existing project safety bounds changed");
        const T maxDrift=contractNumber(safetyText,"max_total_mass_drift_fraction");
        const T maxJump=contractNumber(safetyText,"max_mass_jump_per50_fraction");
        const T maxFlux=contractNumber(safetyText,"max_abs_flux_over_Qtarget");
''' +s[b:]
needle='        require(bool(input),"Incomplete parameter file");'
new=r'''        require(bool(input),"Incomplete parameter file");
        for(int k=1;k<=3;++k)require(ports[k].density==contractNumber(numericalText,"outlet_0"+std::to_string(k)+"_rho_lu"),"Actual outlet density differs from generated contract");
        const char* digest=std::getenv("NEW_MEDIUM_NUMERICS_SHA256");require(digest&&std::string(digest).size()==64,"Missing frozen numerics hash from verified launcher");
        {RankFile audit(run+"/diagnostics/ACTUAL_NUMERICS_READBACK.json");audit<<std::setprecision(17)
            <<"{\"status\":\"PASS\",\"dt_s\":"<<dt<<",\"dx_m\":"<<expectedDx<<",\"tau\":"<<tau<<",\"rho_kg_m3\":"<<rhoPhys
            <<",\"nu_m2_s\":"<<contractNumber(numericalText,"nu_m2_s")<<",\"pressure_unit_pa_read\":"<<contractNumber(numericalText,"pressure_unit_pa")
            <<",\"outlet_01_rho_lu\":"<<ports[1].density<<",\"outlet_02_rho_lu\":"<<ports[2].density<<",\"outlet_03_rho_lu\":"<<ports[3].density
            <<",\"Qtarget\":"<<qTarget<<",\"multiplier\":"<<multiplier<<",\"numerics_contract_sha256\":\""<<digest
            <<"\",\"dt_or_outlet_density_recomputed_in_solver\":false,\"new_medium_inlet_multiplier_validation\":\"NOT_PERFORMED\"}\n";}
'''
assert s.count(needle)==1;s=s.replace(needle,new)
a=s.index('        T initialMass=inspect(nodes).rsum*massUnit;');b=s.index('    } catch(std::exception const& e)',a)
s=s[:a]+r'''        T initialMass=inspect(nodes).rsum*massUnit;
        require(std::isfinite(initialMass)&&initialMass>0,"Invalid initial physical mass");
        RankFile history(run+"/diagnostics/NEW_MEDIUM_RUNTIME_HISTORY.csv");history<<std::setprecision(17)
            <<"iteration,time_s,rho_min,rho_max,rho_mean,Mach_max,ghost_rho_min,ghost_rho_max,ghost_Mach_max,total_mass_kg,total_mass_drift\n";
        RankFile flows(run+"/diagnostics/NEW_MEDIUM_FLUX_HISTORY.csv");flows<<std::setprecision(17)<<"iteration,time_s";
        for(int p=0;p<4;++p)for(int g=0;g<6;++g)flows<<','<<(p==0?"Qin":p==1?"Qout01":p==2?"Qout02":"Qout03")<<"_g"<<g;
        for(int g=0;g<24;++g)flows<<",mass_outward_g"<<g;flows<<'\n';
        RankFile mass(run+"/diagnostics/NEW_MEDIUM_MASS_HISTORY.csv");mass<<std::setprecision(17)
            <<"iteration,time_s,total_mass_kg,total_mass_drift,CV_mass_kg,CV_mass_change_kg,integrated_net_outward_mass_kg,CV_balance_error_kg,CV_balance_error_over_initial_CV_mass\n";
        RankFile safety(run+"/diagnostics/safety_counts.csv");safety<<"iteration,fluid_nan,fluid_inf,ghost_nan,ghost_inf\n";
        RankFile sparse(run+"/diagnostics/snapshot_events.csv");sparse<<"iteration,reason\n";
        const std::set<int> milestones={0,5000};
        uint64_t safetyChecks=0;T allRmin=2,allRmax=0,allMach=0,allGhostMin=2,allGhostMax=0,maxAbsDrift=0;
        T previousMass50=initialMass,initialCV=0,previousCV=0,previousNetMass=0,integratedNetMass=0;
        int previousFluxStep=-1;bool backflow=false;
        double timedStart=0,timedElapsed=0,fullStart=0,fullElapsed=0,segmentStart=0;
        RankFile windows(run+"/diagnostics/iteration_windows.csv");windows<<"first_step,last_step,steps,seconds,steps_per_second\n";
        stage2Profiler.initialize(run);
        for(int step=0;step<=maxSteps;++step){
            stage2Profiler.step=step;
            if(step){state.ramp=std::min(T(1),T(step)/rampSteps);
#ifdef GPU_POC
                gpu.step(state.ramp);
#else
                lattice->collideAndStream();
#endif
                completed=step;
            }
            // Original Stage4 GPU reductions and mathematics, retained every step for safety.
            uint64_t localCounts[3]={uint64_t(nodes.size()),uint64_t(foundCV),uint64_t(sampleNodes.size())},globalCounts[3];
            MPI_Allreduce(localCounts,globalCounts,3,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
            require(globalCounts[0]==actualFluid&&globalCounts[1]==uint64_t(cvExpected)&&globalCounts[2]==sampleIds.size(),"Runtime ownership changed");
            Stats s=inspect(nodes),g=inspect(ghosts);++safetyChecks;
            safety<<step<<','<<s.nan<<','<<s.inf<<','<<g.nan<<','<<g.inf<<'\n';safety.flush();require(bool(safety),"Safety output failed");
            T m=s.rsum*massUnit,cv=s.cvsum*massUnit,drift=(m-initialMass)/initialMass;
            bool safe=!s.nan&&!s.inf&&!g.nan&&!g.inf&&s.rmin>=rhoLow&&s.rmax<=rhoHigh&&g.rmin>0
                &&s.umax/std::sqrt(DESCRIPTOR<T>::cs2)<=maLimit&&g.umax/std::sqrt(DESCRIPTOR<T>::cs2)<=maLimit
                &&std::isfinite(m)&&std::isfinite(cv)&&m>0&&cv>0&&std::abs(drift)<=maxDrift;
            allRmin=std::min(allRmin,s.rmin);allRmax=std::max(allRmax,s.rmax);
            allMach=std::max(allMach,std::max(s.umax,g.umax)/std::sqrt(DESCRIPTOR<T>::cs2));
            allGhostMin=std::min(allGhostMin,g.rmin);allGhostMax=std::max(allGhostMax,g.rmax);maxAbsDrift=std::max(maxAbsDrift,std::abs(drift));
            if(step%50==0||step==maxSteps||!safe){
                safe=safe&&std::abs(m-previousMass50)/initialMass<=maxJump;previousMass50=m;
                history<<step<<','<<step*dt<<','<<s.rmin<<','<<s.rmax<<','<<s.rsum/actualFluid<<','<<s.umax/std::sqrt(DESCRIPTOR<T>::cs2)
                    <<','<<g.rmin<<','<<g.rmax<<','<<g.umax/std::sqrt(DESCRIPTOR<T>::cs2)<<','<<m<<','<<drift<<'\n';history.flush();require(bool(history),"Runtime CSV write failed");
            }
            if(step%100==0||step==maxSteps||!safe){
                checkQuadratureOwnership(samples,run,false);
                T flux[24]={0},mflux[24]={0};
#ifdef GPU_POC
                gpu.flux(flux,mflux,qUnit,rhoPhys);
#else
                throw std::runtime_error("This independent smoke requires the validated GPU implementation");
#endif
                MPI_Allreduce(MPI_IN_PLACE,flux,24,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,mflux,24,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                flows<<step<<','<<step*dt;
                for(int j=0;j<24;++j){safe=safe&&std::isfinite(flux[j])&&std::isfinite(mflux[j])&&std::abs(flux[j])<=maxFlux*qTarget;flows<<','<<(j<6?-1:1)*flux[j];}
                for(int j=0;j<24;++j)flows<<','<<mflux[j];flows<<'\n';flows.flush();require(bool(flows),"Flux CSV write failed");
                T net=mflux[2]+mflux[8]+mflux[14]+mflux[20];
                if(previousFluxStep>=0)integratedNetMass+=.5*(previousNetMass+net)*(step-previousFluxStep)*dt;
                previousFluxStep=step;previousNetMass=net;
                for(int p=1;p<4;++p)backflow=backflow||flux[6*p+2]<0;
                if(step%500==0||step==maxSteps||!safe){
                    if(step==0)initialCV=previousCV=cv;
                    T change=cv-previousCV,error=change+integratedNetMass;
                    mass<<step<<','<<step*dt<<','<<m<<','<<drift<<','<<cv<<','<<change<<','<<integratedNetMass<<','<<error<<','<<error/initialCV<<'\n';mass.flush();require(bool(mass),"Mass CSV write failed");
                    previousCV=cv;integratedNetMass=0;
                }
                if(step%500==0||!safe)pcout<<"SMOKE_STEP "<<step<<" safe="<<safe<<" Qin="<<-flux[2]<<" Qout="<<flux[8]<<","<<flux[14]<<","<<flux[20]<<" Mach="<<std::max(s.umax,g.umax)/std::sqrt(DESCRIPTOR<T>::cs2)<<std::endl;
            }
            if(milestones.count(step)||!safe){snapshot(run,step,nodes);sampledSnapshot(run,step,sampleNodes);sparse<<step<<','<<(safe?"permitted_smoke_snapshot":"failure_snapshot")<<'\n';sparse.flush();}
            require(safe,"Immediate smoke safety failure: nonfinite/density/Mach/mass/flux");
            if(step==0){
                for(int k=0;k<5;++k)require(sumCount(state.calls[k])>0,"A native profile was never used");
                auto ready=std::chrono::steady_clock::now();fullStart=segmentStart=MPI_Wtime();
                RankFile init(run+"/diagnostics/initialization_timing.json");init<<std::setprecision(17)
                    <<"{\"main_entry_monotonic_ns\":"<<std::chrono::duration_cast<std::chrono::nanoseconds>(mainStart.time_since_epoch()).count()
                    <<",\"first_iteration_monotonic_ns\":"<<std::chrono::duration_cast<std::chrono::nanoseconds>(ready.time_since_epoch()).count()
                    <<",\"main_entry_to_ready_seconds\":"<<std::chrono::duration<double>(ready-mainStart).count()
                    <<",\"MPI_initialization_seconds\":"<<std::chrono::duration<double>(mpiReady-mainStart).count()<<",\"includes_step0_required_output\":true}\n";
            }
            if(step==warmup){gpu.sync();MPI_Barrier(MPI_COMM_WORLD);timedStart=MPI_Wtime();}
            stage2Profiler.afterStep(step);
            if(step>0&&step%100==0){double now=MPI_Wtime(),elapsed=now-segmentStart;windows<<std::setprecision(17)<<step-99<<','<<step<<",100,"<<elapsed<<','<<100/elapsed<<'\n';windows.flush();segmentStart=now;}
            if(step==maxSteps){
                gpu.sync();MPI_Barrier(MPI_COMM_WORLD);double end=MPI_Wtime();timedElapsed=end-timedStart;fullElapsed=end-fullStart;
                RankFile timing(run+"/diagnostics/solver_timing.json");timing<<std::setprecision(17)
                    <<"{\"warmup_steps\":"<<warmup<<",\"timed_steps\":"<<maxSteps-warmup<<",\"timed_seconds\":"<<timedElapsed
                    <<",\"steps_per_second\":"<<(maxSteps-warmup)/timedElapsed<<",\"full_iteration_steps\":"<<maxSteps
                    <<",\"full_iteration_seconds\":"<<fullElapsed<<",\"full_iteration_steps_per_second\":"<<maxSteps/fullElapsed<<"}\n";
            }
        }
        stage2Profiler.write();
        RankFile status(run+"/diagnostics/solver_status.json");status<<std::setprecision(17)
            <<"{\"runtime_safety\":\"PASS\",\"status\":\"SMOKE_HORIZON_COMPLETE\",\"timesteps\":"<<completed
            <<",\"auto_converged\":false,\"cell_count\":0,\"mpi_ranks\":1,\"mpi_runtime_correctness\":\"PASS\",\"quadrature_ownership_checks\":"<<ownershipChecks
            <<",\"safety_checks\":"<<safetyChecks<<",\"all_step_fluid_rho_min\":"<<allRmin<<",\"all_step_fluid_rho_max\":"<<allRmax
            <<",\"all_step_max_mach_including_ghosts\":"<<allMach<<",\"all_step_ghost_rho_min\":"<<allGhostMin<<",\"all_step_ghost_rho_max\":"<<allGhostMax
            <<",\"max_abs_mass_drift\":"<<maxAbsDrift<<",\"transient_outlet_backflow\":"<<(backflow?"true":"false")<<"}\n";
        require(bool(status),"Final status write failure");pcout<<"NEW_MEDIUM_SMOKE_COMPLETE steps="<<completed<<" RBC=0"<<std::endl;return 0;
''' +s[b:]
(C/'vascularPoC.cpp').write_text(s)
(C/'DRIVER_ADAPTATION.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='frozen_stage4/vascularPoC.cpp',tofile='new_medium_host_driver/vascularPoC.cpp')))
# Exact identity of all fluid/Guo setup between completed parameter read and GPU state construction.
a=old.index('        // Step2 closed helper');b=old.index('        T initialMass=inspect(nodes)')
assert old[a:b] in s
assert old[old.index('class RampedNativeProfile'):old.index('int main(')] in s or old[old.index('class RampedNativeProfile'):old.index('int main(')].rstrip() in s
(C/'DRIVER_STATIC_PROOF.json').write_text(json.dumps(dict(status='PASS',baseline_source_sha256=hashlib.sha256(old.encode()).hexdigest(),new_driver_sha256=hashlib.sha256(s.encode()).hexdigest(),geometry_BC_Guo_lattice_GPU_initialization_block_byte_identical=True,RampedNativeProfile_byte_identical=True,changed_scope=['Generated contract read/identity audit','500/5000 task horizons','Host monitoring/output cadence and smoke-only safety','Only0/5000/failure snapshots'],gpu_core_headers_changed=False,solver_regenerates_dt_or_outlet_rho=False),indent=2)+'\n')
