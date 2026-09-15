// Stage 2 scratch instrumentation. No additional CUDA synchronization per stage.
#pragma once
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <sys/mman.h>
#ifdef GPU_POC
#include <nvtx3/nvToolsExt.h>
#include <cuda_profiler_api.h>
#endif
struct Stage2Record { const char* name; int step; long long start,end; };
struct Stage2Profiler {
    Stage2Record* rows=nullptr; size_t count=0; static constexpr size_t capacity=20000;
    int step=0; bool enabled=false,trace=false,traceOpen=false;
    std::string run;
    static long long now(){return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();}
    void initialize(std::string const& path) {
        run=path;const char* m=std::getenv("STAGE2_PROFILE_MODE");
        require(m&&(std::strcmp(m,"off")==0||std::strcmp(m,"timers")==0||std::strcmp(m,"trace")==0),"Explicit Stage2 instrumentation mode required");
        enabled=std::strcmp(m,"off")!=0;trace=std::strcmp(m,"trace")==0;
        if(enabled){void* p=mmap(nullptr,capacity*sizeof(Stage2Record),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);require(p!=MAP_FAILED,"Profiling buffer allocation failure");rows=static_cast<Stage2Record*>(p);}
    }
    bool active()const{return enabled&&step>100&&step<=200;}
    void afterStep(int s) {
#ifdef GPU_POC
        if(trace&&s==100){require(cudaProfilerStart()==cudaSuccess,"Profiler start failed");nvtxRangePushA("STAGE3_TRACE_WINDOW");traceOpen=true;}
        if(trace&&s==200){nvtxRangePop();traceOpen=false;require(cudaProfilerStop()==cudaSuccess,"Profiler stop failed");}
#endif
    }
    void write() {
        if(global::mpi().getRank()!=0)return;
        std::ofstream f(run+"/diagnostics/host_stage_ranges.csv");f<<"iteration,component,start_ns,end_ns,duration_ns\n";
        for(size_t i=0;i<count;++i){auto const&r=rows[i];f<<r.step<<','<<r.name<<','<<r.start<<','<<r.end<<','<<r.end-r.start<<'\n';}
        require(bool(f),"Profiling evidence write failure");
    }
    ~Stage2Profiler(){if(rows)munmap(rows,capacity*sizeof(Stage2Record));}
} stage2Profiler;
struct Stage2Scope {
    const char* name; bool active; long long start=0;
    explicit Stage2Scope(const char* n):name(n),active(stage2Profiler.active()){
        if(!active)return;
#ifdef GPU_POC
        nvtxRangePushA(name);
#endif
        start=Stage2Profiler::now();
    }
    ~Stage2Scope(){
        if(!active)return;auto end=Stage2Profiler::now();
#ifdef GPU_POC
        nvtxRangePop();
#endif
        require(stage2Profiler.count<Stage2Profiler::capacity,"Profiling buffer overflow");
        stage2Profiler.rows[stage2Profiler.count++]={name,stage2Profiler.step,start,end};
    }
};
