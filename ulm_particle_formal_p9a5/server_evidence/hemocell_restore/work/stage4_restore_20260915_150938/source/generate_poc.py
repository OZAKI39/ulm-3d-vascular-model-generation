from pathlib import Path
import hashlib,json,difflib
D=Path(__file__).resolve().parent
s=(D/'frozen_reference/inletCalibrationMPI.cpp').read_text();base=s
s=s.replace('#include <vector>','#include <vector>\n#include <chrono>')
s=s.replace('#include "mpi_support.hpp"','#ifdef GPU_POC\n#include "gpu_adapter.hpp"\n#endif\n#include "mpi_support.hpp"')
a=s.index('        require(argc==6');b=s.index('        run=argv[1];',a)
s=s[:a]+'''        require(argc==4,"Usage: vascularPoC RUN_DIR BUNDLE_ROOT MULTIPLIER");
#ifdef GPU_POC
        require(global::mpi().getSize()==1,"GPU MPI1 only");
        defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);
#else
        require(global::mpi().getSize()==12,"CPU MPI12 only");
#endif
        std::string root=argv[2],test="GPU_POC_STAGE_1";
        T multiplier=std::stod(argv[3]);require(multiplier==T(1.1197286861799598),"Frozen multiplier mismatch");
        const int requested=1000;
'''+s[b:]
s=s.replace('maxSteps<=10000','maxSteps==1000').replace('Long validation forbidden in calibration','PoC hard cap 1000')
s=s.replace('auto* model=new GuoOffLatticeModel3D<T,DESCRIPTOR>(new TriangleFlowShape3D<T,Vec>(domain.getBoundary(),profiles),voxelFlag::inside,true);','''#ifdef GPU_POC
        Capture capture{&state,{}};
        auto* model=new CapturingGuo(new TriangleFlowShape3D<T,Vec>(domain.getBoundary(),profiles),&capture);
#else
        auto* model=new GuoOffLatticeModel3D<T,DESCRIPTOR>(new TriangleFlowShape3D<T,Vec>(domain.getBoundary(),profiles),voxelFlag::inside,true);
#endif''')
s=s.replace('bc(model,domain,*lattice); bc.insert();','''bc(model,domain,*lattice);
#ifndef GPU_POC
        bc.insert();
#endif''')
s=s.replace('lattice->initialize();','''lattice->initialize();
#ifdef GPU_POC
        bc.apply(); // Native initial Guo completion; no CPU Guo processor copied into GPU lattice.
#endif''')
s=s.replace('        T initialMass=inspect(nodes).rsum*massUnit;','''#ifdef GPU_POC
        GPUState gpu(*lattice,capture,nx,ny,nz,nodes,ghosts,samples,sampleIds,run);activeGPU=&gpu;gpu.sync();
#endif
        T initialMass=inspect(nodes).rsum*massUnit;''')
s=s.replace('std::set<int> milestones={0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};','std::set<int> milestones={0,100,1000};')
s=s.replace('if(step) {state.ramp=std::min(T(1),T(step)/rampSteps);lattice->collideAndStream();completed=step;}','''if(step) {state.ramp=std::min(T(1),T(step)/rampSteps);
#ifdef GPU_POC
                gpu.step(state.ramp);
#else
                lattice->collideAndStream();
#endif
                completed=step;
            }''')
a=s.index('            if(step==10){');b=s.index('            uint64_t localCounts',a)
s=s[:a]+'''            // Initialization verifies exact ownership on both backends; no CPU-only distributed proxy diagnostic in the fair PoC.
'''+s[b:]
a=s.index('                auto macros=sampledMacros(sampleIds,sampleNodes);');b=s.index('                MPI_Allreduce(MPI_IN_PLACE,flux',a)
s=s[:a]+'''#ifdef GPU_POC
                gpu.flux(flux,mflux,qUnit,rhoPhys);
#else
'''+s[a:b]+'''#endif
'''+s[b:]
s=s.replace('                if(step>=1000)mpiFieldResidual(run,step,step-1000,nodes,dx,dt,rhoPhys);','')
a=s.index('            if(test=="A"&&step==5000){');b=s.index('            if(step==maxSteps)',a)
s=s[:a]+'''            if(step==100){
#ifdef GPU_POC
                gpu.sync();
#endif
                MPI_Barrier(MPI_COMM_WORLD);timedStart=MPI_Wtime();
            }
            if(step==1000){
#ifdef GPU_POC
                gpu.sync();
#endif
                MPI_Barrier(MPI_COMM_WORLD);timedElapsed=MPI_Wtime()-timedStart;
                RankFile timing(run+"/diagnostics/solver_timing.json");timing<<std::setprecision(17)<<"{\\"warmup_steps\\":100,\\"timed_steps\\":900,\\"timed_seconds\\":"<<timedElapsed<<",\\"steps_per_second\\":"<<900/timedElapsed<<"}\\n";
            }
'''+s[b:]
s=s.replace('SHORT_CALIBRATION_HORIZON_COMPLETE','POC_1000_COMPLETE').replace('STEP3C_SHORT_COMPLETE','GPU_POC_SHORT_COMPLETE').replace('short_calibration_snapshot','poc_same_output_snapshot')
# Numerical source stays literal except backend dispatch and bounded run/IO/timing.
(D/'vascularPoC.cpp').write_text(s)
h=(D/'frozen_reference/mpi_support.hpp').read_text();hb=h
h=h.replace('    Stats s;\n    for(auto const& n:nodes)', '''#ifdef GPU_POC
    if(activeGPU)return activeGPU->inspect(!nodes.empty()&& !nodes.front().inCV && nodes.front().index==activeGPU->ghost.front().id);
#endif
    Stats s;
    for(auto const& n:nodes)''',1)
# Avoid relying on node order/CV flags to identify fluid vs ghost sets.
h=h.replace('!nodes.empty()&& !nodes.front().inCV && nodes.front().index==activeGPU->ghost.front().id','!nodes.empty()&& !activeGPU->ghost.empty() && nodes.front().index==activeGPU->ghost.front().id')
h=h.replace('    std::vector<PackedNode> local(nodes.size());','''#ifdef GPU_POC
    if(activeGPU){std::vector<uint64_t> ids;for(auto const& n:nodes)ids.push_back(n.index);auto v=activeGPU->pack(ids);std::vector<PackedNode> out(v.size());std::memcpy(out.data(),v.data(),out.size()*sizeof(PackedNode));std::sort(out.begin(),out.end(),[](PackedNode const&a,PackedNode const&b){return a.index<b.index;});return out;}
#endif
    std::vector<PackedNode> local(nodes.size());''',1)
(D/'mpi_support.hpp').write_text(h)
p=(D/'gpu_adapter.hpp');g=p.read_text().replace('BGKdynamics<T,DESCRIPTOR>::id','BGKdynamics<T,DESCRIPTOR>(1.).getId()').replace('index(d+2*v)','index(Dot3D(d.x+2*v.x,d.y+2*v.y,d.z+2*v.z))');p.write_text(g)
patch=''.join(difflib.unified_diff(base.splitlines(True),s.splitlines(True),fromfile='frozen_reference/inletCalibrationMPI.cpp',tofile='source/vascularPoC.cpp'))+''.join(difflib.unified_diff(hb.splitlines(True),h.splitlines(True),fromfile='frozen_reference/mpi_support.hpp',tofile='source/mpi_support.hpp'))+''.join(difflib.unified_diff([],g.splitlines(True),fromfile='/dev/null',tofile='source/gpu_adapter.hpp'))
(D/'vascular_adapter.diff').write_text(patch)
print('GENERATED',len(s.splitlines()),'driver lines;',len(g.splitlines()),'adapter lines')
