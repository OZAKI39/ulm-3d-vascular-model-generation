// Pure fluid only. No HemoCell/Palabos core changes or population reconstruction.
// Closed native geometry stages reuse Step2 closedVascularVoxelizer.cpp;
// source SHA256: 0cdba2795f2b435abbdaa1ded5dab1bc9342588f0ff8fa6f6b618cdcd100aaad.
// Native Guo setup follows bundled Palabos examples/showCases/aneurysm.
#include "palabos3D.h"
#include "palabos3D.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <sys/wait.h>
#include <unistd.h>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
#include <chrono>
using namespace plb;
typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor
typedef Array<T,3> Vec;
typedef Cell<T,DESCRIPTOR> LCell;

void require(bool ok, std::string const& why) { if(!ok) throw std::runtime_error(why); }
struct Port { int label; Vec normal; T density; std::vector<int> ids; std::vector<std::array<Vec,3> > faces; };
struct State { T ramp=0, inletSpeed=0; uint64_t calls[5]={0,0,0,0,0}; };

// Public boundary DATA composition only. Native profiles and native Guo own all
// density/velocity-to-population reconstruction. Clone shares startup state.
class RampedNativeProfile : public BoundaryProfile3D<T,Vec> {
    int label; Vec normal; T density; State* state;
public:
    RampedNativeProfile(int label_,Vec n,T rho,State* s):label(label_),normal(n),density(rho),state(s) {}
    void setNormal(Vec const& n) override { normal=n; }
    void defineCircularShape(Vec const&,T) override {} // no circular cap construction
    void getData(Vec const& pos,plint id,AtomicBlock3D const* arg,Vec& data,OffBoundary::Type& type) const override {
        ++state->calls[label];
        if(label==0) { NoSlipProfile3D<T> p; p.getData(pos,id,arg,data,type); }
        else if(label==1) {
            VelocityPlugProfile3D<T> p(state->ramp*state->inletSpeed); p.setNormal(normal);
            p.getData(pos,id,arg,data,type);
        } else {
            DensityNeumannBoundaryProfile3D<T> p(1+state->ramp*(density-1));
            p.getData(pos,id,arg,data,type);
        }
    }
    RampedNativeProfile* clone() const override { return new RampedNativeProfile(*this); }
};
struct ExactVertexSet {
    std::set<std::array<T,3> > vertices;
    bool operator()(Vec const& p) const { return vertices.count({{p[0],p[1],p[2]}})>0; }
};
std::vector<uint8_t> raw(std::string const& path,size_t n) {
    std::ifstream f(path,std::ios::binary); std::vector<uint8_t> v(n);
    f.read(reinterpret_cast<char*>(v.data()),n); require(bool(f)&&f.peek()==EOF,"Raw field size/read failed: "+path); return v;
}
struct FluidNode { uint64_t index; LCell* cell; bool inCV=false; FluidNode(uint64_t i,LCell* c):index(i),cell(c){} };
struct Sample { int port,owner; T area; Vec normal; std::array<uint64_t,8> indices; std::array<T,8> weights; };
struct Stats { T rmin=2,rmax=0,rsum=0,umax=0,cvsum=0,speedSum=0; uint64_t nan=0,inf=0; };
#include "stage2_profile.hpp"
#ifdef GPU_POC
#include "gpu_adapter.hpp"
#endif
#include "mpi_support.hpp"

int main(int argc,char** argv) {
    const auto mainStart=std::chrono::steady_clock::now();
    plbInit(&argc,&argv);
    const auto mpiReady=std::chrono::steady_clock::now();
    std::string run; int completed=0;
    try {
        require(argc==5,"Usage: vascularPoC RUN_DIR BUNDLE_ROOT MULTIPLIER STEPS");
#ifdef GPU_POC
        require(global::mpi().getSize()==1,"GPU MPI1 only");
        defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);
#else
        require(global::mpi().getSize()==12,"CPU MPI12 only");
#endif
        std::string root=argv[2],test="GPU_POC_STAGE_4";
        T multiplier=std::stod(argv[3]);require(multiplier==T(1.1197286861799598),"Frozen multiplier mismatch");
        const int requested=std::stoi(argv[4]);
        require(requested==200||requested==1000||requested==5000,"Stage4 frozen horizon mismatch");
        const int warmup=requested==200?100:requested/10;
        run=argv[1]; global::directories().setOutputDir(run+"/");
        pcout<<std::setprecision(17);
        std::ifstream input(run+"/contracts/solver_parameters.txt"); std::string stl,s2;
        input>>stl>>s2;
        int nx,ny,nz,rampSteps,maxSteps; Vec expectedOrigin; T expectedDx,dt,tau,rhoPhys,qTarget,maLimit,rhoLow,rhoHigh;
        input>>nx>>ny>>nz>>expectedOrigin[0]>>expectedOrigin[1]>>expectedOrigin[2]>>expectedDx>>dt>>tau>>rhoPhys>>qTarget>>maLimit>>rhoLow>>rhoHigh>>rampSteps>>maxSteps;
        require(maxSteps==700000&&rampSteps==10&&tau==1.,"Frozen execution contract mismatch");
        maxSteps=requested; // Independent Step3C execution cap; frozen Step3B file stays byte-identical.
        require(maxSteps<=5000,"Stage4 hard cap 5000");
        int scalarEvery,evaluateEvery,fieldDelta,firstEvaluation,cvExpected;
        std::ifstream monitorInput(run+"/contracts/monitor_parameters.txt");
        monitorInput>>scalarEvery>>evaluateEvery>>fieldDelta>>firstEvaluation>>cvExpected;
        require(bool(monitorInput)&&scalarEvery==100&&evaluateEvery==5000&&fieldDelta==119872&&firstEvaluation==240000&&cvExpected==180543,"Monitor contract mismatch");
        std::vector<Port> ports(4); std::map<int,int> expectedTags;
        for(auto& p:ports) {
            int count; input>>p.label>>p.normal[0]>>p.normal[1]>>p.normal[2]>>p.density>>count;
            for(int i=0;i<count;++i) {
                int id; input>>id; std::array<Vec,3> face;
                for(auto& v:face) input>>v[0]>>v[1]>>v[2];
                require(!expectedTags.count(id),"Overlapping cap IDs"); expectedTags[id]=p.label;
                p.ids.push_back(id); p.faces.push_back(face);
            }
        }
        require(bool(input),"Incomplete parameter file");
        // Step2 closed helper stages verbatim in numerical meaning; envelope 2
        // accommodates native Guo's second neighbor and is not a mesh change.
        TriangleSet<T> triangles(stl,DBL);
        DEFscaledMesh<T> scaled(triangles,494,1,1,0);
        const Vec origin=scaled.getPhysicalLocation(); const T dx=scaled.getDx();
        require(std::abs(dx-expectedDx)<1e-20 && norm(origin-expectedOrigin)<1e-17,"Step2 transform mismatch");
        TriangleBoundary3D<T> boundary(scaled,false);
        auto& mesh=boundary.getMesh();
        std::vector<int> tags;
        for(auto const& p:ports) {
            ExactVertexSet domain;
            for(size_t k=0;k<p.ids.size();++k) {
                int id=p.ids[k]; require(id<mesh.getNumTriangles(),"Cap ID outside mesh");
                for(int j=0;j<3;++j) {
                    Vec v=mesh.getVertex(id,j); T nearest=std::numeric_limits<T>::max();
                    for(int c=0;c<3;++c) nearest=std::min(nearest,norm(v-(p.faces[k][c]-origin)/dx));
                    require(nearest<1e-8,"Step2 cap triangle coordinates/IDs differ in native mesh");
                    domain.vertices.insert({{v[0],v[1],v[2]}});
                }
            }
            int tag=boundary.tagDomain(domain); require(tag==p.label,"Native cap tag mismatch"); tags.push_back(tag);
        }
        for(int i=0;i<mesh.getNumTriangles();++i) {
            int expected=expectedTags.count(i)?expectedTags[i]:0;
            require(boundary.getTag(i)==expected,"TagDomain affected a non-cap triangle or missed a cap");
        }
        boundary.getMesh().inflate(); // same official 0.001 LU used and audited in Step2
        State state;
        T inletProjectedArea=0,unmodifiedArea=0;
        RankFile capAudit(run+"/diagnostics/native_cap_area.csv"); capAudit<<std::setprecision(17)<<"label,projected_native_area_lu2,frozen_area_m2\n";
        for(auto const& p:ports) {
            T projected=0,areaM=0;
            for(size_t k=0;k<p.ids.size();++k) {
                int i=p.ids[k]; Vec a=mesh.getVertex(i,0),b=mesh.getVertex(i,1),c=mesh.getVertex(i,2);
                projected+=dot(crossProduct(b-a,c-a),p.normal)/2;
                auto const& face=p.faces[k]; areaM+=norm(crossProduct(face[1]-face[0],face[2]-face[0]))/2;
            }
            require(projected>0,"Cap orientation disagrees with outward normal");
            capAudit<<p.label<<','<<projected<<','<<areaM<<'\n';
            if(p.label==1) { inletProjectedArea=projected; unmodifiedArea=areaM; }
        }
        T qUnit=dx*dx*dx/dt,velUnit=dx/dt,massUnit=rhoPhys*dx*dx*dx;
        state.inletSpeed=(qTarget/qUnit)/inletProjectedArea;
        T integralError=std::abs(state.inletSpeed*inletProjectedArea*qUnit-qTarget)/qTarget;
        require(integralError<=1e-12,"Inlet profile does not integrate to Qtarget");
        // Only physical change in command application: fixed numerical inlet compensation.
        // qTarget and nominal integral above remain the original physical target.
        T nominalInletSpeed=state.inletSpeed;
        state.inletSpeed*=multiplier;
        {RankFile f(run+"/diagnostics/inlet_command.json");f<<std::setprecision(17)
          <<"{\"test\":\""<<test<<"\",\"physical_Qtarget\":"<<qTarget
          <<",\"multiplier\":"<<multiplier<<",\"nominal_velocity_LU\":"<<nominalInletSpeed
          <<",\"numerical_velocity_LU\":"<<state.inletSpeed<<",\"nominal_integral_error\":"<<integralError
          <<",\"numerical_profile_integral_m3_s\":"<<multiplier*qTarget<<"}\n";}
        BoundaryProfiles3D<T,Vec> profiles;
        profiles.setWallProfile(new RampedNativeProfile(0,Vec(0.,0.,0.),1.,&state));
        for(auto const& p:ports) profiles.defineProfile(p.label,new RampedNativeProfile(p.label,p.normal,p.density,&state));
        pcout<<"NATIVE_CAP_TAGS_PASS triangles="<<expectedTags.size()<<" profile_integral_error="<<integralError<<std::endl;
        VoxelizedDomain3D<T> domain(boundary,voxelFlag::inside,0,1,2,16);
        auto& voxels=domain.getVoxelMatrix(); Box3D box=voxels.getBoundingBox();
        require(box.x0==0&&box.y0==0&&box.z0==0&&box.getNx()==nx&&box.getNy()==ny&&box.getNz()==nz,"Step2 lattice bbox changed");
        size_t count=size_t(nx)*ny*nz;
        auto native=raw(s2+"/diagnostics/palabos_native_flags.u8",count);
        auto closed=raw(s2+"/diagnostics/closed_flag_matrix.u8",count);
        auto labels=raw(s2+"/diagnostics/port_label_field.u8",count);
        auto opened=raw(s2+"/diagnostics/opened_flag_matrix.u8",count);
        auto const& bulks=voxels.getMultiBlockManagement().getSparseBlockStructure().getBulks();
        uint64_t actualFluid=0,labelCounts[5]={0,0,0,0,0},allocated=0;
        for(auto const& entry:bulks) {
            if(!voxels.getMultiBlockManagement().getThreadAttribution().isLocal(entry.first)) continue;
            auto const& v=voxels.getComponent(entry.first); Dot3D loc=v.getLocation(); auto b=entry.second;
            for(plint z=b.z0;z<=b.z1;++z) for(plint y=b.y0;y<=b.y1;++y) for(plint x=b.x0;x<=b.x1;++x) {
                size_t i=(z*ny+y)*nx+x; int flag=v.get(x-loc.x,y-loc.y,z-loc.z);
                require(flag==native[i],"Native voxel flag changed from Step2");
                bool fluid=flag==voxelFlag::inside||flag==voxelFlag::innerBorder;
                require(fluid==bool(closed[i]),"Closed lumen flag mismatch"); actualFluid+=fluid; ++allocated;
            }
        }
        actualFluid=sumCount(actualFluid);allocated=sumCount(allocated);
        for(size_t i=0;i<count;++i) {
            require(labels[i]<=4,"Invalid frozen port label");
            require(opened[i]==uint8_t(closed[i]||labels[i]),"Frozen opened/closed/port relation invalid");
            if(labels[i]) { require(native[i]==voxelFlag::outerBorder,"Port diagnostic node is not native boundary support"); ++labelCounts[labels[i]]; }
        }
        require(actualFluid==182694&&labelCounts[1]==194&&labelCounts[2]==85&&labelCounts[3]==93&&labelCounts[4]==125,"Step2 counts changed");
        RankFile geometryAudit(run+"/diagnostics/runtime_geometry_check.json"); geometryAudit<<std::setprecision(17);
        geometryAudit<<"{\"status\":\"PASS\",\"native_bulk_flags_equal_step2\":true,\"closed_lumen_equal_step2\":true,\"port_labels_unchanged\":true,\"cap_triangle_tags_exact\":true,\"physical_fluid_voxels\":"<<actualFluid<<",\"allocated_bulk_voxels\":"<<allocated<<",\"native_envelope\":2,\"ghost_nodes_are_physical_fluid\":false,\"profile_integral_error\":"<<integralError<<",\"inlet_native_projected_area_m2\":"<<inletProjectedArea*dx*dx<<",\"inlet_frozen_area_m2\":"<<unmodifiedArea<<"}\n"; geometryAudit.close();
        closed.clear();closed.shrink_to_fit();opened.clear();opened.shrink_to_fit();labels.clear();labels.shrink_to_fit();
        auto lattice=generateMultiBlockLattice<T,DESCRIPTOR>(voxels,2,new BGKdynamics<T,DESCRIPTOR>(1/tau));
        lattice->periodicity().toggleAll(false); lattice->toggleInternalStatistics(false);
        defineDynamics(*lattice,voxels,lattice->getBoundingBox(),new NoDynamics<T,DESCRIPTOR>,voxelFlag::outside);
        #ifdef GPU_POC
        Capture capture{&state,{}};
        auto* model=new CapturingGuo(new TriangleFlowShape3D<T,Vec>(domain.getBoundary(),profiles),&capture);
#else
        auto* model=new GuoOffLatticeModel3D<T,DESCRIPTOR>(new TriangleFlowShape3D<T,Vec>(domain.getBoundary(),profiles),voxelFlag::inside,true);
#endif
        model->setVelIsJ(false);model->selectSecondOrder(true);model->selectUseRegularizedModel(true);model->selectComputeStat(false);
        OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Vec> bc(model,domain,*lattice);
#ifndef GPU_POC
        bc.insert();
#endif
        initializeAtEquilibrium(*lattice,lattice->getBoundingBox(),1.,Vec(0.,0.,0.)); lattice->initialize();
#ifdef GPU_POC
        bc.apply(); // Native initial Guo completion; no CPU Guo processor copied into GPU lattice.
#endif
        std::vector<FluidNode> nodes,ghosts;
        for(auto const& entry:bulks) {
            if(!voxels.getMultiBlockManagement().getThreadAttribution().isLocal(entry.first)) continue;
            auto& block=lattice->getComponent(entry.first); Dot3D loc=block.getLocation(); auto b=entry.second;
            for(plint z=b.z0;z<=b.z1;++z) for(plint y=b.y0;y<=b.y1;++y) for(plint x=b.x0;x<=b.x1;++x) {
                size_t i=(z*ny+y)*nx+x; auto* cell=&block.get(x-loc.x,y-loc.y,z-loc.z);
                if(native[i]==3||native[i]==4) nodes.push_back({i,cell});
                else if(native[i]==2) ghosts.push_back({i,cell});
            }
        }
        require(sumCount(nodes.size())==actualFluid,"FAIL_MPI_RUNTIME_CORRECTNESS: fluid cache global count mismatch");
        std::map<uint64_t,LCell*> sampleNodes;
        std::map<LCell*,uint64_t> pointerOwners;
        std::set<uint64_t> sampleIdSet;
        std::vector<Sample> samples; std::ifstream quad(run+"/contracts/multiplane_quadrature.tsv"); int label; T x,y,z,a,n0,n1,n2;
        auto const& management=lattice->getMultiBlockManagement();auto const& attribution=management.getThreadAttribution();
        while(quad>>label>>x>>y>>z>>a>>n0>>n1>>n2){
            require(label>=0&&label<24,"Quadrature group outside frozen24 groups");Sample sample;sample.port=label;sample.area=a;sample.normal=Vec(n0,n1,n2);
            int base[3]={int(std::floor(x)),int(std::floor(y)),int(std::floor(z))};T frac[3]={x-base[0],y-base[1],z-base[2]};
            plint pointBlock=management.getSparseBlockStructure().locate(base[0],base[1],base[2]);require(pointBlock>=0,"Quadrature point missing bulk owner");
            sample.owner=attribution.getMpiProcess(pointBlock);require(sample.owner>=0&&sample.owner<rankCount(),"FAIL_MPI_RUNTIME_CORRECTNESS: invalid quadrature owner");int j=0;
            for(int ox=0;ox<2;++ox)for(int oy=0;oy<2;++oy)for(int oz=0;oz<2;++oz){
                int ix=base[0]+ox,iy=base[1]+oy,iz=base[2]+oz;require(ix>=0&&ix<nx&&iy>=0&&iy<ny&&iz>=0&&iz<nz,"Flux sample outside grid");
                uint64_t id=(uint64_t(iz)*ny+iy)*nx+ix;require(native[id]>=1&&native[id]<=4,"Flux sample unallocated");
                sample.indices[j]=id;sampleIdSet.insert(id);sample.weights[j]=(ox?frac[0]:1-frac[0])*(oy?frac[1]:1-frac[1])*(oz?frac[2]:1-frac[2]);++j;
            }
            samples.push_back(sample);
        }
        require(quad.eof()&&!samples.empty(),"Incomplete flux quadrature");
        std::vector<uint64_t> sampleIds(sampleIdSet.begin(),sampleIdSet.end());
        std::vector<int> localSampleOwnership(sampleIds.size(),0),globalSampleOwnership(sampleIds.size());
        for(size_t slot=0;slot<sampleIds.size();++slot){
            uint64_t id=sampleIds[slot];int ix=id%nx,iy=(id/nx)%ny,iz=id/(uint64_t(nx)*ny);
            plint blockId=management.getSparseBlockStructure().locate(ix,iy,iz);require(blockId>=0,"Sample cell missing bulk owner");
            if(!attribution.isLocal(blockId))continue;
            auto const& box=management.getSparseBlockStructure().getBulks().at(blockId);
            require(ix>=box.x0&&ix<=box.x1&&iy>=box.y0&&iy<=box.y1&&iz>=box.z0&&iz<=box.z1,"FAIL_MPI_RUNTIME_CORRECTNESS: nonbulk pointer");
            auto& block=lattice->getComponent(blockId);Dot3D loc=block.getLocation();LCell* cell=&block.get(ix-loc.x,iy-loc.y,iz-loc.z);
            require(!pointerOwners.count(cell)||pointerOwners[cell]==id,"FAIL_MPI_RUNTIME_CORRECTNESS: persistent atomic pointer alias");
            pointerOwners[cell]=id;sampleNodes[id]=cell;localSampleOwnership[slot]=1;
        }
        MPI_Allreduce(localSampleOwnership.data(),globalSampleOwnership.data(),int(sampleIds.size()),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        for(int owners:globalSampleOwnership)require(owners==1,"FAIL_MPI_RUNTIME_CORRECTNESS: sample cell ownership");
        require(sumCount(pointerOwners.size())==sampleIds.size(),"FAIL_MPI_RUNTIME_CORRECTNESS: global pointer alias count");
        for(auto& sample:samples)for(auto& id:sample.indices)id=std::lower_bound(sampleIds.begin(),sampleIds.end(),id)-sampleIds.begin();
        checkQuadratureOwnership(samples,run,true);native.clear();native.shrink_to_fit();

        std::vector<uint64_t> cvIds;uint64_t id;
        std::ifstream cvFile(run+"/contracts/control_volume_indices.txt");
        while(cvFile>>id) cvIds.push_back(id);
        require(cvFile.eof()&&cvIds.size()==size_t(cvExpected)&&std::is_sorted(cvIds.begin(),cvIds.end())&&std::adjacent_find(cvIds.begin(),cvIds.end())==cvIds.end(),"Invalid frozen CV indices");
        size_t foundCV=0;
        for(auto& node:nodes) {node.inCV=std::binary_search(cvIds.begin(),cvIds.end(),node.index);foundCV+=node.inCV;}
        require(sumCount(foundCV)==cvIds.size(),"FAIL_MPI_RUNTIME_CORRECTNESS: CV global count mismatch");
        {RankFile check(run+"/diagnostics/control_volume_ownership_check.json");check<<"{\"status\":\"PASS\",\"CV_CELL_GLOBAL_COUNT\":"<<cvExpected<<",\"physical_fluid_global_count\":"<<actualFluid<<"}\n";}
#ifdef GPU_POC
        GPUState gpu(*lattice,capture,nx,ny,nz,nodes,ghosts,samples,sampleIds,run,expectedTags);activeGPU=&gpu;gpu.sync();
#endif
        T initialMass=inspect(nodes).rsum*massUnit;
        RankFile history(run+"/diagnostics/flow_history.csv");history<<std::setprecision(17);
        history<<"iteration,time_s,rho_min,rho_max,rho_mean,u_max_m_s,Mach_max,mean_speed_m_s,total_mass,control_volume_mass,relative_mass_drift,BOUNDARY_PROFILE_Q_TARGET";
        const char* portNames[4]={"Qin","Qout01","Qout02","Qout03"};
        for(int p=0;p<4;++p) for(int k=0;k<3;++k) for(int r=0;r<2;++r)
            history<<','<<portNames[p]<<'_'<<(2+2*k)<<"dx"<<(r?"_dx2":"");
        for(int g=0;g<24;++g) history<<",mass_outward_g"<<g;
        history<<",instantaneous_flow_closure,control_volume_mass_balance,ghost_rho_min,ghost_rho_max,ghost_Mach_max\n";
        RankFile safety(run+"/diagnostics/safety_counts.csv");safety<<"iteration,fluid_nan,fluid_inf,ghost_nan,ghost_inf\n";
        RankFile sparse(run+"/diagnostics/snapshot_events.csv");sparse<<"iteration,reason\n";
        std::set<int> milestones={0,100,200,1000,5000};
        if(std::getenv("STAGE4_STRONG_CHECKPOINTS")){milestones.insert(1);milestones.insert(10);}
        uint64_t safetyChecks=0;T allRmin=2,allRmax=0,allMach=0,allGhostMin=2,allGhostMax=0;
        T previousCV=0,previousMassFlux[4]={0,0,0,0};int previousSample=-1;
        bool autoConverged=false;T window3sum=0,window4sum=0;int window3count=0,window4count=0;
        double timedStart=0,timedElapsed=0,fullStart=0,fullElapsed=0,segmentStart=0;
        long long firstIterationNs=0;
        RankFile windows(run+"/diagnostics/iteration_windows.csv");
        windows<<"first_step,last_step,steps,seconds,steps_per_second\n";
        stage2Profiler.initialize(run);
        for(int step=0;step<=maxSteps;++step) {
            stage2Profiler.step=step;
            { Stage2Scope wholeStep("TIMESTEP");
            if(step) {state.ramp=std::min(T(1),T(step)/rampSteps);
#ifdef GPU_POC
                gpu.step(state.ramp);
#else
                lattice->collideAndStream();
#endif
                completed=step;
            }
            // Initialization verifies exact ownership on both backends; no CPU-only distributed proxy diagnostic in the fair PoC.
            {Stage2Scope p("MPI_OWNERSHIP_SCALARS");
            uint64_t localCounts[3]={uint64_t(nodes.size()),uint64_t(foundCV),uint64_t(sampleNodes.size())},globalCounts[3];MPI_Allreduce(localCounts,globalCounts,3,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
            require(globalCounts[0]==actualFluid&&globalCounts[1]==uint64_t(cvExpected)&&globalCounts[2]==sampleIds.size(),"FAIL_MPI_RUNTIME_CORRECTNESS: runtime owned cell count changed");}
            Stats s=inspect(nodes),g=inspect(ghosts);++safetyChecks;
            {Stage2Scope p("OUTPUT_SAFETY");safety<<step<<','<<s.nan<<','<<s.inf<<','<<g.nan<<','<<g.inf<<'\n';safety.flush();require(bool(safety),"Safety evidence write failure");}
            bool safe=!s.nan&&!s.inf&&!g.nan&&!g.inf&&s.rmin>=rhoLow&&s.rmax<=rhoHigh&&g.rmin>0&&s.umax/std::sqrt(DESCRIPTOR<T>::cs2)<=maLimit&&g.umax/std::sqrt(DESCRIPTOR<T>::cs2)<=maLimit;
            allRmin=std::min(allRmin,s.rmin);allRmax=std::max(allRmax,s.rmax);
            allMach=std::max(allMach,std::max(s.umax,g.umax)/std::sqrt(DESCRIPTOR<T>::cs2));
            allGhostMin=std::min(allGhostMin,g.rmin);allGhostMax=std::max(allGhostMax,g.rmax);
            { // Exact discrete arithmetic means: measure every completed step, including all24groups.
                T flux[24]={0},mflux[24]={0};
                if(step%scalarEvery==0||!safe){Stage2Scope p("MPI_QUADRATURE_OWNERSHIP");checkQuadratureOwnership(samples,run,false);}
#ifdef GPU_POC
                gpu.flux(flux,mflux,qUnit,rhoPhys);
#else
                auto macros=sampledMacros(sampleIds,sampleNodes);
                for(auto const& sample:samples){
                    if(sample.owner!=rankId())continue;
                    Vec u(0.,0.,0.);T rho=0;
                    for(int j=0;j<8;++j){size_t k=4*sample.indices[j];Vec v(macros[k+1],macros[k+2],macros[k+3]);u+=sample.weights[j]*v;rho+=sample.weights[j]*macros[k];}
                    T dq=sample.area*dot(u,sample.normal)*qUnit;flux[sample.port]+=dq;mflux[sample.port]+=rhoPhys*rho*dq;
                }
#endif
                {Stage2Scope p("MPI_FLUX_SCALARS");MPI_Allreduce(MPI_IN_PLACE,flux,24,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,mflux,24,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);}
                for(int g=0;g<24;++g)require(std::isfinite(flux[g])&&std::isfinite(mflux[g]),"FAIL_MPI_RUNTIME_CORRECTNESS: flux reduction nonfinite");
                T qi=-flux[2],net=flux[2]+flux[8]+flux[14]+flux[20];
                if(step>=3001&&step<=4000){window3sum+=qi;++window3count;}
                if(step>=4001&&step<=5000){window4sum+=qi;++window4count;}
                T mass=s.rsum*massUnit,cvMass=s.cvsum*massUnit;
                {Stage2Scope p("OUTPUT_FLOW_HISTORY");history<<step<<','<<step*dt<<','<<s.rmin<<','<<s.rmax<<','<<s.rsum/actualFluid<<','<<s.umax*velUnit<<','<<s.umax/std::sqrt(DESCRIPTOR<T>::cs2)<<','<<s.speedSum/actualFluid*velUnit<<','<<mass<<','<<cvMass<<','<<(mass-initialMass)/initialMass<<','<<qTarget*state.ramp;
                for(int group=0;group<24;++group) history<<','<<(group<6?-1:1)*flux[group];
                for(int group=0;group<24;++group) history<<','<<mflux[group];
                history<<','<<std::abs(net)/std::max(std::abs(qi),qTarget*1e-12)<<',';
                if(previousSample>=0) {
                    T inward=0,outwardSum=0,netMass=0;
                    for(int p=0;p<4;++p) {T mean=.5*(previousMassFlux[p]+mflux[6*p+2]);netMass+=mean;if(p==0) inward=std::abs(mean);else outwardSum+=std::abs(mean);}
                    T residual=std::abs((cvMass-previousCV)/((step-previousSample)*dt)+netMass)/std::max(std::max(inward,outwardSum),rhoPhys*qTarget*1e-12);
                    history<<residual;
                }
                history<<','<<g.rmin<<','<<g.rmax<<','<<g.umax/std::sqrt(DESCRIPTOR<T>::cs2)<<'\n';history.flush();require(bool(history),"History write failure");}
                previousSample=step;previousCV=cvMass;for(int p=0;p<4;++p) previousMassFlux[p]=mflux[6*p+2];
                if(step%1000==0||!safe) pcout<<"STAGE "<<step<<" safe="<<safe<<" Qin="<<qi<<" CV_mass="<<cvMass<<std::endl;
            }
            if(!safe) throw std::runtime_error("Immediate runtime safety gate failure: NaN/Inf/density/Mach");
            if(milestones.count(step)){Stage2Scope p("OUTPUT_SNAPSHOTS");
                snapshot(run,step,nodes);sampledSnapshot(run,step,sampleNodes);
                sparse<<step<<",poc_same_output_snapshot\n";sparse.flush();

            }
            if(step==0)for(int k=0;k<5;++k)require(sumCount(state.calls[k])>0,"A native surface profile was never used globally");
            } // End measured step; excludes profiler control and timing-file serialization.
            if(step==0){
                // The existing initial GPU synchronization, scalar checks and snapshots have completed.
                firstIterationNs=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
                fullStart=segmentStart=MPI_Wtime();
                RankFile init(run+"/diagnostics/initialization_timing.json");init<<std::setprecision(17)
                    <<"{\"main_entry_monotonic_ns\":"<<std::chrono::duration_cast<std::chrono::nanoseconds>(mainStart.time_since_epoch()).count()
                    <<",\"first_iteration_monotonic_ns\":"<<firstIterationNs
                    <<",\"main_entry_to_ready_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-mainStart).count()
                    <<",\"MPI_initialization_seconds\":"<<std::chrono::duration<double>(mpiReady-mainStart).count()
                    <<",\"includes_step0_required_output\":true}\n";
            }
            if(step==warmup){
#ifdef GPU_POC
                gpu.sync();
#endif
                MPI_Barrier(MPI_COMM_WORLD);timedStart=MPI_Wtime();
            }
            stage2Profiler.afterStep(step);
            if(step>0&&step%100==0){
                double now=MPI_Wtime(),elapsed=now-segmentStart;
                windows<<std::setprecision(17)<<step-99<<','<<step<<",100,"<<elapsed<<','<<100/elapsed<<'\n';windows.flush();
                require(bool(windows),"Window timing evidence write failure");segmentStart=now;
            }
            if(step==maxSteps){
#ifdef GPU_POC
                gpu.sync();
#endif
                MPI_Barrier(MPI_COMM_WORLD);const double end=MPI_Wtime();timedElapsed=end-timedStart;fullElapsed=end-fullStart;
                RankFile timing(run+"/diagnostics/solver_timing.json");timing<<std::setprecision(17)
                    <<"{\"warmup_steps\":"<<warmup<<",\"timed_steps\":"<<maxSteps-warmup
                    <<",\"timed_seconds\":"<<timedElapsed<<",\"steps_per_second\":"<<(maxSteps-warmup)/timedElapsed
                    <<",\"full_iteration_steps\":"<<maxSteps<<",\"full_iteration_seconds\":"<<fullElapsed
                    <<",\"full_iteration_steps_per_second\":"<<maxSteps/fullElapsed<<"}\n";
            }
            if(step==maxSteps){sparse<<step<<",final\n";sparse.flush();break;}
        }
        stage2Profiler.write();
        for(int k=0;k<5;++k)state.calls[k]=sumCount(state.calls[k]);
        RankFile status(run+"/diagnostics/solver_status.json");status<<std::setprecision(17);
        status<<"{\"runtime_safety\":\"PASS\",\"status\":\""<<"STAGE3_HORIZON_COMPLETE"<<"\",\"timesteps\":"<<completed<<",\"auto_converged\":"<<(autoConverged?"true":"false")<<",\"cell_count\":0,\"mpi_ranks\":"<<rankCount()<<",\"mpi_runtime_correctness\":\"PASS\",\"quadrature_ownership_checks\":"<<ownershipChecks<<",\"safety_checks\":"<<safetyChecks<<",\"all_step_fluid_rho_min\":"<<allRmin<<",\"all_step_fluid_rho_max\":"<<allRmax<<",\"all_step_max_mach_including_ghosts\":"<<allMach<<",\"all_step_ghost_rho_min\":"<<allGhostMin<<",\"all_step_ghost_rho_max\":"<<allGhostMax<<",\"native_profile_calls\":[";
        for(int k=0;k<5;++k) status<<(k?",":"")<<state.calls[k];status<<"]}\n";
        require(bool(status),"Final status write failed");
        pcout<<"GPU_POC_STAGE3_COMPLETE timesteps="<<completed<<" auto_converged="<<autoConverged<<" CELL_COUNT=0"<<std::endl;return 0;
    } catch(std::exception const& e) {
        pcout<<"STEP3C_ERROR completed_timesteps="<<completed<<" "<<e.what()<<std::endl;
        if(!run.empty()) {RankFile f(run+"/diagnostics/solver_failure.txt");f<<"completed_timesteps="<<completed<<"\n"<<e.what()<<'\n';}
        std::cerr<<"RANK "<<rankId()<<" STEP3C_MPI_ERROR "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;
    }
}
