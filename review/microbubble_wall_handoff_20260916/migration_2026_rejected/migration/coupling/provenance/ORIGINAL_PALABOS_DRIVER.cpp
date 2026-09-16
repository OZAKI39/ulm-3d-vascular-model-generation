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
#include <cctype>
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
        std::string root=argv[2],test="PURE_FLUID_NEW_MEDIUM_SMOKE";
        T multiplier=std::stod(argv[3]);require(multiplier==T(1.1197286861799598),"Frozen multiplier mismatch");
        const int requested=std::stoi(argv[4]);
        require(requested==500||requested==5000,"New-medium smoke horizon must be500 or5000");
        const int warmup=requested/10;
        run=argv[1]; global::directories().setOutputDir(run+"/");
        pcout<<std::setprecision(17);
        std::ifstream input(run+"/contracts/solver_parameters.txt"); std::string stl,s2;
        input>>stl>>s2;
        int nx,ny,nz,rampSteps,maxSteps; Vec expectedOrigin; T expectedDx,dt,tau,rhoPhys,qTarget,maLimit,rhoLow,rhoHigh;
        input>>nx>>ny>>nz>>expectedOrigin[0]>>expectedOrigin[1]>>expectedOrigin[2]>>expectedDx>>dt>>tau>>rhoPhys>>qTarget>>maLimit>>rhoLow>>rhoHigh>>rampSteps>>maxSteps;
        require(maxSteps==5000&&rampSteps==10&&tau==1.,"Smoke execution contract mismatch");
        maxSteps=requested; // New independent smoke cap; original Stage4/Step3C contracts remain unchanged.
        require(maxSteps<=5000,"Stage4 hard cap 5000");
        const int scalarEvery=100,cvExpected=180543;
        std::string numericalText=contractText(run+"/contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json");
        std::string safetyText=contractText(run+"/contracts/SMOKE_SAFETY_CONTRACT.json");
        require(dt==contractNumber(numericalText,"dt_s")&&expectedDx==contractNumber(numericalText,"dx_m")&&tau==contractNumber(numericalText,"tau")&&rhoPhys==contractNumber(numericalText,"rho_kg_m3"),"Actual input differs from generated numerics contract");
        require(qTarget==contractNumber(numericalText,"physical_Qtarget_m3_s")&&multiplier==contractNumber(numericalText,"inlet_numerical_multiplier"),"Physical inlet mismatch");
        require(rhoPhys==1000.&&contractNumber(numericalText,"nu_m2_s")==1e-6&&tau==1.&&dt>0,"Provisional medium input mismatch");
        require(rhoLow==.99&&rhoHigh==1.01&&maLimit==.05,"Existing project safety bounds changed");
        const T maxDrift=contractNumber(safetyText,"max_total_mass_drift_fraction");
        const T maxJump=contractNumber(safetyText,"max_mass_jump_per50_fraction");
        const T maxFlux=contractNumber(safetyText,"max_abs_flux_over_Qtarget");
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
        for(int k=1;k<=3;++k)require(ports[k].density==contractNumber(numericalText,"outlet_0"+std::to_string(k)+"_rho_lu"),"Actual outlet density differs from generated contract");
        const char* digest=std::getenv("NEW_MEDIUM_NUMERICS_SHA256");require(digest&&std::string(digest).size()==64,"Missing frozen numerics hash from verified launcher");
        {RankFile audit(run+"/diagnostics/ACTUAL_NUMERICS_READBACK.json");audit<<std::setprecision(17)
            <<"{\"status\":\"PASS\",\"dt_s\":"<<dt<<",\"dx_m\":"<<expectedDx<<",\"tau\":"<<tau<<",\"rho_kg_m3\":"<<rhoPhys
            <<",\"nu_m2_s\":"<<contractNumber(numericalText,"nu_m2_s")<<",\"pressure_unit_pa_read\":"<<contractNumber(numericalText,"pressure_unit_pa")
            <<",\"outlet_01_rho_lu\":"<<ports[1].density<<",\"outlet_02_rho_lu\":"<<ports[2].density<<",\"outlet_03_rho_lu\":"<<ports[3].density
            <<",\"Qtarget\":"<<qTarget<<",\"multiplier\":"<<multiplier<<",\"numerics_contract_sha256\":\""<<digest
            <<"\",\"dt_or_outlet_density_recomputed_in_solver\":false,\"new_medium_inlet_multiplier_validation\":\"NOT_PERFORMED\"}\n";}

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
    } catch(std::exception const& e) {
        pcout<<"STEP3C_ERROR completed_timesteps="<<completed<<" "<<e.what()<<std::endl;
        if(!run.empty()) {RankFile f(run+"/diagnostics/solver_failure.txt");f<<"completed_timesteps="<<completed<<"\n"<<e.what()<<'\n';}
        std::cerr<<"RANK "<<rankId()<<" STEP3C_MPI_ERROR "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;
    }
}
