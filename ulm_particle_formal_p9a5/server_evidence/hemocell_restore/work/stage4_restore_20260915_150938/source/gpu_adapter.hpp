// Scratch-only adapter for the frozen, unforced D3Q19 BGK / GuoPiNeq case.
// Algorithm source: Palabos guoOffLatticeModel3D.hh lines 296-562 (AGPL-3.0-or-later).
// Original geometry preparation and native CPU initialization are retained.
#include <execution>
#include <numeric>
#include <cuda_runtime.h>
using AL=AtomicAcceleratedLattice3D<T,DESCRIPTOR>;
using ML=AcceleratedLattice3D<T,DESCRIPTOR>;
using DD=DESCRIPTOR<T>;
using MM=momentTemplatesImpl<T,typename DD::BaseDescriptor>;
using DT=dynamicsTemplatesImpl<T,typename DD::BaseDescriptor>;
using Pi=Array<T,6>;
// Scratch MPI1 full-node halo plan. Native communicator defines every overlap.
// A single device kernel executes the verified plan at each of the two halo stages.
struct FrozenHaloJob {AL const* source;AL* destination;plint sourceIndex,destinationIndex;plint sourceId,destinationId;};
struct FrozenHaloBlock {plint id;Box3D bulkLocal;};
struct FrozenHaloPlan {
    std::vector<FrozenHaloJob> jobs;
    std::map<AL const*,FrozenHaloBlock> blocks;
    bool recording=true;
    size_t overlaps=0;
    void record(AL* dest,AL const* src,Box3D to,plint dx,plint dy,plint dz) {
        if(!recording)return;
        require(src&&dest,"Invalid accelerated halo blocks");
        auto sb=src->getBoundingBox();auto db=dest->getBoundingBox();
        auto srcMeta=blocks.at(src),dstMeta=blocks.at(dest);
        require(contained(to.x0,to.y0,to.z0,db)&&contained(to.x1,to.y1,to.z1,db)
             &&contained(to.x0+dx,to.y0+dy,to.z0+dz,sb)&&contained(to.x1+dx,to.y1+dy,to.z1+dz,sb),"Halo overlap exceeds original block bounds");
        require(contained(to.x0+dx,to.y0+dy,to.z0+dz,srcMeta.bulkLocal)
             &&contained(to.x1+dx,to.y1+dy,to.z1+dz,srcMeta.bulkLocal),"Halo source must belong to original source bulk");
        plint sy=src->getNy(),sz=src->getNz(),dySize=dest->getNy(),dzSize=dest->getNz();
        for(plint x=to.x0;x<=to.x1;++x)for(plint y=to.y0;y<=to.y1;++y)for(plint z=to.z0;z<=to.z1;++z){
            require(!contained(x,y,z,dstMeta.bulkLocal),"Halo destination must be in original envelope");
            jobs.push_back({src,dest,(z+dz)+sz*((y+dy)+sy*(x+dx)),z+dzSize*(y+dySize*x),srcMeta.id,dstMeta.id});
        }
        ++overlaps;
    }
    void seal(std::string const& run) {
        auto less=[](FrozenHaloJob const& a,FrozenHaloJob const& b){return a.destinationId<b.destinationId||(a.destinationId==b.destinationId&&a.destinationIndex<b.destinationIndex);};
        std::sort(jobs.begin(),jobs.end(),less);
        size_t before=jobs.size();
        for(size_t i=1;i<jobs.size();++i){
            auto const& a=jobs[i-1];auto const& b=jobs[i];
            if(a.destinationId==b.destinationId&&a.destinationIndex==b.destinationIndex)
                require(a.sourceId==b.sourceId&&a.sourceIndex==b.sourceIndex,"Conflicting native halo destinations; cannot batch");
        }
        jobs.erase(std::unique(jobs.begin(),jobs.end(),[](FrozenHaloJob const& a,FrozenHaloJob const& b){return a.destinationId==b.destinationId&&a.destinationIndex==b.destinationIndex;}),jobs.end());
        require(!jobs.empty()&&overlaps>0,"Missing original halo plan");
        recording=false;
        std::ofstream f(run+"/diagnostics/halo_plan_inventory.json");
        f<<"{\"status\":\"PASS\",\"native_overlaps\":"<<overlaps<<",\"original_copy_nodes\":"<<before<<",\"unique_copy_nodes\":"<<jobs.size()<<",\"idempotent_duplicate_copies\":"<<(before-jobs.size())<<",\"source_bulk_verified\":true,\"destination_envelope_verified\":true,\"conflicting_writes\":0,\"populations_per_node\":19,\"halo_kernels_per_step\":2}\n";
        require(bool(f),"Halo plan inventory write failure");
    }
    void apply() {
        Stage2Scope p("HALO_BATCH");
        require(!recording,"Unsealed halo plan");
        std::for_each(std::execution::par_unseq,jobs.begin(),jobs.end(),[](FrozenHaloJob const& n){
            Array<T,19> f;n.source->pullPop(n.sourceIndex,f);n.destination->pushPop(n.destinationIndex,f);
        });
    }
};
FrozenHaloPlan* activeHaloPlan=nullptr;
class FrozenFullHaloTransfer : public AcceleratedLatticeDataTransfer3D<T,DESCRIPTOR> {
    using Base=AcceleratedLatticeDataTransfer3D<T,DESCRIPTOR>;
    AL* destination=nullptr;
public:
    using Base::attribute;
    FrozenFullHaloTransfer* clone() const override {return new FrozenFullHaloTransfer(*this);}
    void setBlock(AtomicBlock3D& block) override {
        Base::setBlock(block);destination=dynamic_cast<AL*>(&block);
        require(destination,"Halo destination must be accelerated lattice");
    }
    void attribute(Box3D to,plint dx,plint dy,plint dz,AtomicBlock3D const& block,modif::ModifT kind) override {
        require(activeHaloPlan&&kind==modif::staticVariables,"Frozen halo adapter supports MPI1 static populations only");
        static_assert(DD::ExternalField::numScalars==0,"Frozen unforced descriptor only");
        if(activeHaloPlan->recording)activeHaloPlan->record(destination,dynamic_cast<AL const*>(&block),to,dx,dy,dz);
    }
};

struct CapturedLink {plint one,two;int depth,type,triangle;T delta,weight;Vec normal,fullData;};
struct CapturedDry {BlockLattice3D<T,DESCRIPTOR>* cpu;plint index;std::vector<CapturedLink> links;};
struct Capture {State* state;std::vector<CapturedDry> dry;};
class CapturingGuo: public GuoOffLatticeModel3D<T,DESCRIPTOR> {
    Capture* capture;
public:
    CapturingGuo(BoundaryShape3D<T,Vec>* shape,Capture* c):GuoOffLatticeModel3D<T,DESCRIPTOR>(shape,voxelFlag::inside,true),capture(c){}
    CapturingGuo* clone() const override{return new CapturingGuo(*this);}
    void boundaryCompletion(AtomicBlock3D& a,AtomicContainerBlock3D& container,std::vector<AtomicBlock3D*> const& args) override {
        require(args.empty()&&usesRegularizedModel()&&usesSecondOrder()&&!computesStat()&&!velIsJ()&&!getPartialReplace(),"Unsupported Guo mode; no semantic substitution");
        GuoOffLatticeModel3D<T,DESCRIPTOR>::boundaryCompletion(a,container,args);
        auto& block=dynamic_cast<BlockLattice3D<T,DESCRIPTOR>&>(a);
        auto* info=dynamic_cast<GuoOffLatticeInfo3D*>(container.getData());require(info,"Missing original Guo links");
        auto offset=container.getLocation();T oldRamp=capture->state->ramp;capture->state->ramp=1;
        auto index=[&](Dot3D const& p){require(contained(p.x,p.y,p.z,block.getBoundingBox()),"Guo link exceeds original two-cell envelope");return p.z+block.getNz()*(p.y+block.getNy()*p.x);};
        for(size_t n=0;n<info->getDryNodes().size();++n){
            auto d=info->getDryNodes()[n];CapturedDry cd{&block,index(d),{}};
            auto const& dirs=info->getDryNodeFluidDirections()[n];auto const& ids=info->getDryNodeIds()[n];
            require(block.get(d.x,d.y,d.z).getDynamics().getId()==BGKdynamics<T,DESCRIPTOR>(1.).getId(),"Dry regularization requires original BGK dynamics");
            for(size_t k=0;k<dirs.size();++k){
                auto dir=dirs[k];auto c=NextNeighbor<T>::c[dir.first];Dot3D v(c[0],c[1],c[2]);
                Vec wall,normal,data;T distance;OffBoundary::Type type;plint tri=ids[k];
                require(pointOnSurface(d+offset,v,wall,distance,normal,data,type,tri),"Original Guo surface intersection failed");
                require(type==OffBoundary::dirichlet||type==OffBoundary::densityNeumann,"Unexpected BC in frozen vascular case");
                T inv=NextNeighbor<T>::invD[dir.first];require(distance<=NextNeighbor<T>::d[dir.first],"Invalid original Guo distance");
                Vec nv{T(c[0]),T(c[1]),T(c[2])};nv*=inv;
                require(block.get(d.x+c[0],d.y+c[1],d.z+c[2]).getDynamics().getId()==BGKdynamics<T,DESCRIPTOR>(1.).getId(),"First Guo neighbor is not original fluid BGK");
                cd.links.push_back({index(d+v),index(Dot3D(d.x+2*v.x,d.y+2*v.y,d.z+2*v.z)),dir.second,int(type),int(tri),1-distance*inv,std::fabs(dot(nv,normal)),normal,data});
            }
            require(!cd.links.empty()&&cd.links.size()<=26,"Invalid original direction count");capture->dry.push_back(cd);
        }
        capture->state->ramp=oldRamp;
    }
};
struct Handle {AL* block;plint index;uint64_t id;bool inCV;};
struct DeviceDry {AL* block;plint index;size_t begin,end;};
struct DevicePacked {uint64_t index;T rho;T u[3];};
static_assert(sizeof(DevicePacked)==40,"Snapshot layout");
struct FluxPair {T q=0,m=0;};
struct DeviceSample {int group;T area;Vec normal;Handle cells[8];T weights[8];};
inline void moments(Handle const& h,T& rb,Vec& j){
    if(h.block->getCollisionMatrix()[h.index]==CollisionModel::NoDynamics){rb=0;j.resetToZero();return;}
    Array<T,19> f;h.block->pullPop(h.index,f);MM::get_rhoBar_j(f,rb,j);
}
inline DevicePacked macro(Handle const& h){T rb;Vec j;moments(h,rb,j);T rho=DD::fullRho(rb),inv=DD::invRho(rb);return {h.id,rho,{j[0]*inv,j[1]*inv,j[2]*inv}};}
struct ReduceStats {Stats operator()(Stats const& a,Stats const& b)const{Stats r;r.rmin=std::min(a.rmin,b.rmin);r.rmax=std::max(a.rmax,b.rmax);r.rsum=a.rsum+b.rsum;r.umax=std::max(a.umax,b.umax);r.cvsum=a.cvsum+b.cvsum;r.speedSum=a.speedSum+b.speedSum;r.nan=a.nan+b.nan;r.inf=a.inf+b.inf;return r;}};
#include "stage4_lbm_batch.hpp"
struct GPUState {
    Stage4LBMBatch lbmBatch;bool useBatch=false;
    FrozenHaloPlan halo;
    std::unique_ptr<ML> lattice;std::vector<DeviceDry> dry;std::vector<CapturedLink> links;std::vector<Handle> fluid,ghost;std::vector<DeviceSample> quadrature;std::array<size_t,25> groupOffset{};
    int nx,ny,nz;uint64_t boundaryKernels=0,monitorCalls=0,snapshotCalls=0;std::string run;
    Handle handle(uint64_t id,bool cv=false){
        int x=id%nx,y=(id/nx)%ny,z=id/(uint64_t(nx)*ny);plint b=lattice->getMultiBlockManagement().getSparseBlockStructure().locate(x,y,z);require(b>=0,"Unallocated GPU sample");auto& a=lattice->getComponent(b);auto o=a.getLocation();return {&a,(z-o.z)+a.getNz()*((y-o.y)+a.getNy()*(x-o.x)),id,cv};
    }
    GPUState(MultiBlockLattice3D<T,DESCRIPTOR>& cpu,Capture const& capture,int x,int y,int z,std::vector<FluidNode> const& nodes,std::vector<FluidNode> const& ghosts,std::vector<Sample> const& samples,std::vector<uint64_t> const& sampleIds,std::string path,std::map<int,int> const& capTags):nx(x),ny(y),nz(z),run(path){
        require(global::mpi().getSize()==1,"Single GPU/MPI1 only");lattice.reset(new ML(cpu));
        std::map<BlockLattice3D<T,DESCRIPTOR>*,AL*> mapping;
        activeHaloPlan=&halo;
        auto const& originalBulks=lattice->getMultiBlockManagement().getSparseBlockStructure().getBulks();
        for(auto const& b:cpu.getBlockLattices()){
            auto& a=lattice->getComponent(b.first);auto o=a.getLocation();auto g=originalBulks.at(b.first);
            halo.blocks[&a]={b.first,Box3D(g.x0-o.x,g.x1-o.x,g.y0-o.y,g.y1-o.y,g.z0-o.z,g.z1-o.z)};
            a.setDataTransfer(new FrozenFullHaloTransfer);mapping[b.second]=&a;
        }
        lattice->duplicateOverlaps(modif::staticVariables);halo.seal(run);halo.apply();
        for(auto const& d:capture.dry){size_t begin=links.size();links.insert(links.end(),d.links.begin(),d.links.end());dry.push_back({mapping.at(d.cpu),d.index,begin,links.size()});}
        for(auto const& n:nodes)fluid.push_back(handle(n.index,n.inCV));
        for(auto const& n:ghosts)ghost.push_back(handle(n.index));
        for(int g=0;g<24;++g){groupOffset[g]=quadrature.size();for(auto const& s:samples)if(s.port==g){DeviceSample q;q.group=g;q.area=s.area;q.normal=s.normal;for(int k=0;k<8;++k){q.cells[k]=handle(sampleIds.at(s.indices[k]));q.weights[k]=s.weights[k];}quadrature.push_back(q);}}groupOffset[24]=quadrature.size();
        require(fluid.size()==182694&&quadrature.size()==samples.size()&&!dry.empty(),"GPU static inventory mismatch");
        {
            std::map<AL*,size_t> fluidCounts;
            for(auto const& h:fluid)++fluidCounts[h.block];
            std::ofstream out(run+"/diagnostics/gpu_block_distribution.csv");
            out<<"block_id,nx,ny,nz,allocated_cells,physical_fluid_owned_cells,bulk_cells\n";
            for(auto const& b:halo.blocks){auto* a=b.first;out<<b.second.id<<','<<a->getNx()<<','<<a->getNy()<<','<<a->getNz()<<','<<a->getBoundingBox().nCells()<<','<<fluidCounts[const_cast<AL*>(a)]<<','<<b.second.bulkLocal.nCells()<<'\n';}
            require(bool(out),"Block inventory write failure");
            size_t linkCounts[5]={0},nodeCounts[5]={0},pureNodes[5]={0},mixedNodes=0;
            for(auto const& d:capture.dry){unsigned mask=0;for(auto const& l:d.links){auto it=capTags.find(l.triangle);int label=it==capTags.end()?0:it->second;require(label>=0&&label<5,"Invalid captured boundary label");++linkCounts[label];mask|=1u<<label;}
                for(int k=0;k<5;++k){nodeCounts[k]+=bool(mask&(1u<<k));pureNodes[k]+=mask==(1u<<k);}mixedNodes+=(mask&(mask-1))!=0;
            }
            std::ofstream bc(run+"/diagnostics/guo_boundary_inventory.csv");
            bc<<"label,boundary,links,nodes_touching_category,pure_category_nodes,mixed_nodes_total\n";
            const char* names[5]={"Guo_wall","inlet","outlet_01","outlet_02","outlet_03"};
            for(int k=0;k<5;++k)bc<<k<<','<<names[k]<<','<<linkCounts[k]<<','<<nodeCounts[k]<<','<<pureNodes[k]<<','<<mixedNodes<<'\n';
            require(bool(bc),"Guo inventory write failure");
        }
        {
            std::ofstream out(run+"/diagnostics/gpu_address_inventory.csv");
            out<<"block_id,region,address,bytes,cuda_pointer_type\n";
            auto emit=[&](plint id,const char* region,const void* pointer,size_t bytes){
                cudaPointerAttributes attr{};auto rc=cudaPointerGetAttributes(&attr,pointer);
                int type=rc==cudaSuccess?int(attr.type):-1;if(rc!=cudaSuccess)cudaGetLastError();
                out<<id<<','<<region<<','<<reinterpret_cast<uintptr_t>(pointer)<<','<<bytes<<','<<type<<'\n';
            };
            for(auto const& b:halo.blocks){auto* a=b.first;
                emit(b.second.id,"AL_METADATA",a,sizeof(AL));
                emit(b.second.id,"BACKGROUND_DYNAMICS",&a->getBackgroundDynamics(),sizeof(BGKdynamics<T,DESCRIPTOR>));
                emit(b.second.id,"COLLISION_MATRIX",a->getCollisionMatrix(),a->getN()*sizeof(int));
            }
            emit(-1,"HALO_JOBS",halo.jobs.data(),halo.jobs.size()*sizeof(FrozenHaloJob));
            emit(-1,"GUO_DRY",dry.data(),dry.size()*sizeof(DeviceDry));
            emit(-1,"GUO_LINKS",links.data(),links.size()*sizeof(CapturedLink));
            emit(-1,"FLUID_HANDLES",fluid.data(),fluid.size()*sizeof(Handle));
            emit(-1,"GHOST_HANDLES",ghost.data(),ghost.size()*sizeof(Handle));
            emit(-1,"QUADRATURE",quadrature.data(),quadrature.size()*sizeof(DeviceSample));
            require(bool(out),"Address inventory write failure");
        }
        std::ofstream f(run+"/diagnostics/gpu_static_inventory.json");f<<"{\"dry_nodes\":"<<dry.size()<<",\"boundary_links\":"<<links.size()<<",\"fluid_nodes\":"<<fluid.size()<<",\"ghost_nodes\":"<<ghost.size()<<",\"quadrature_points\":"<<quadrature.size()<<",\"blocks\":"<<mapping.size()<<",\"full_field_writeBack_calls\":0,\"core_source_changes\":0}\n";
        // Micro 1: residency hints for small shared AL objects only.
        // Frozen arrays, collision/Guo mathematics, execution order and reductions are unchanged.
        {
            cudaMemLocation hostLocation{};hostLocation.type=cudaMemLocationTypeHost;
            cudaMemLocation deviceLocation{};deviceLocation.type=cudaMemLocationTypeDevice;deviceLocation.id=0;
            std::ofstream evidence(run+"/diagnostics/metadata_residency_advice.csv");
            evidence<<"block_id,address,bytes,preferred_location,accessed_by,return_code\n";
            for(auto const& entry:halo.blocks){
                auto* object=entry.first;cudaPointerAttributes attr{};
                require(cudaPointerGetAttributes(&attr,object)==cudaSuccess&&attr.type==cudaMemoryTypeManaged,"Metadata advice requires observed managed allocation");
                auto r1=cudaMemAdvise(object,sizeof(AL),cudaMemAdviseSetPreferredLocation,hostLocation);
                auto r2=cudaMemAdvise(object,sizeof(AL),cudaMemAdviseSetAccessedBy,deviceLocation);
                require(r1==cudaSuccess&&r2==cudaSuccess,"Metadata residency advice failed");
                evidence<<entry.second.id<<','<<reinterpret_cast<uintptr_t>(object)<<','<<sizeof(AL)<<",CPU,GPU0,0\n";
            }
            require(bool(evidence),"Residency advice evidence write failure");
        }
        const char* variant=std::getenv("STAGE4_LBM");
        require(variant&&(std::string(variant)=="baseline"||std::string(variant)=="batch"),"Explicit Stage4 LBM variant required");
        useBatch=std::string(variant)=="batch";
        lbmBatch.initialize(*lattice,run,useBatch);
    }
    void sync(){require(cudaDeviceSynchronize()==cudaSuccess,"CUDA synchronize failed");}
    void step(T ramp){
        {Stage2Scope p("LBM_COLLIDE_STREAM_FRAMEWORK");if(useBatch)lattice->stage4CollideAndStream(CollisionKernel<T,DESCRIPTOR,CollisionModel::BGK,CollisionModel::NoDynamics>(),lbmBatch);else lattice->collideAndStream(CollisionKernel<T,DESCRIPTOR,CollisionModel::BGK,CollisionModel::NoDynamics>());}
        halo.apply();
        {Stage2Scope p("GUO_FUSED_WALL_INLET_OUTLETS");
        auto* ls=links.data();
        std::for_each(std::execution::par_unseq,dry.begin(),dry.end(),[ls,ramp](DeviceDry const& n){
            T rb=0,weightSum=0;Vec j(0.,0.,0.);Pi pi;pi.resetToZero();
            for(size_t k=n.begin;k<n.end;++k){
                auto const& l=ls[k];Array<T,19> f;T rb1;Vec j1,j2;Pi p1;
                n.block->pullPop(l.one,f);MM::compute_rhoBar_j_PiNeq(f,rb1,j1,p1);
                T tmp;moments({n.block,l.two,0,false},tmp,j2);
                Vec wd=l.fullData;if(l.type==OffBoundary::densityNeumann){wd[0]=1+ramp*(wd[0]-1);}else wd*=ramp;
                T r=l.type==OffBoundary::densityNeumann?DD::rhoBar(wd[0]):rb1;Vec wj=DD::fullRho(r)*wd,jd;
                if(l.depth<2){jd=wj;} // Both frozen depth<2 branches use wall_j (native complex-geometry fix).
                else if(l.delta<T(.75)){jd=wj+(l.delta-T(1))*j1+(T(1)-l.delta)/(T(1)+l.delta)*(T(2)*wj+(l.delta-T(1))*j2);}
                else{jd=T(1)/l.delta*(wj+(l.delta-T(1))*j1);}
                if(l.type==OffBoundary::densityNeumann)jd=dot(j1,l.normal)*l.normal;
                rb+=r*l.weight;j+=jd*l.weight;pi+=p1*l.weight;weightSum+=l.weight;
            }
            rb/=weightSum;j/=weightSum;pi/=weightSum;T js=normSqr(j),inv=DD::invRho(rb);Array<T,19> f;
            f[0]=DT::bgk_ma2_equilibrium(0,rb,inv,j,js)+offEquilibriumTemplates<T,DESCRIPTOR>::fromPiToFneq(0,pi);
            for(int i=1;i<=9;++i){f[i]=DT::bgk_ma2_equilibrium(i,rb,inv,j,js);f[i+9]=DT::bgk_ma2_equilibrium(i+9,rb,inv,j,js);T neq=offEquilibriumTemplates<T,DESCRIPTOR>::fromPiToFneq(i,pi);f[i]+=neq;f[i+9]+=neq;}
            n.block->pushPop(n.index,f);
        });
        }
        {Stage2Scope p("HALO_FRAMEWORK_POST_GUO");lattice->duplicateOverlaps(modif::staticVariables);}
        halo.apply();++boundaryKernels;
    }
    Stats inspect(bool ghosts){Stage2Scope p(ghosts?"MONITOR_GHOST_RHO_MACH_MASS":"MONITOR_FLUID_RHO_MACH_MASS");auto const& list=ghosts?ghost:fluid;++monitorCalls;return std::transform_reduce(std::execution::par_unseq,list.begin(),list.end(),Stats{},ReduceStats{},[](Handle const& h){auto m=macro(h);T speed=std::sqrt(m.u[0]*m.u[0]+m.u[1]*m.u[1]+m.u[2]*m.u[2]);Stats s;s.rmin=m.rho;s.rmax=m.rho;s.rsum=m.rho;s.cvsum=h.inCV?m.rho:0;s.umax=speed;s.speedSum=speed;s.nan=std::isnan(m.rho)||std::isnan(speed);s.inf=std::isinf(m.rho)||std::isinf(speed);return s;});}
    std::vector<DevicePacked> pack(std::vector<uint64_t> const& ids){std::vector<Handle> hs;for(auto id:ids)hs.push_back(handle(id));std::vector<DevicePacked> out(ids.size());std::transform(std::execution::par_unseq,hs.begin(),hs.end(),out.begin(),[](Handle const& h){return macro(h);});sync();++snapshotCalls;return out;}
    void flux(T* flux,T* mflux,T qunit,T rhoPhys){
        static const char* names[24]={"MONITOR_FLUX_G00","MONITOR_FLUX_G01","MONITOR_FLUX_G02","MONITOR_FLUX_G03","MONITOR_FLUX_G04","MONITOR_FLUX_G05","MONITOR_FLUX_G06","MONITOR_FLUX_G07","MONITOR_FLUX_G08","MONITOR_FLUX_G09","MONITOR_FLUX_G10","MONITOR_FLUX_G11","MONITOR_FLUX_G12","MONITOR_FLUX_G13","MONITOR_FLUX_G14","MONITOR_FLUX_G15","MONITOR_FLUX_G16","MONITOR_FLUX_G17","MONITOR_FLUX_G18","MONITOR_FLUX_G19","MONITOR_FLUX_G20","MONITOR_FLUX_G21","MONITOR_FLUX_G22","MONITOR_FLUX_G23"};
        for(int g=0;g<24;++g){Stage2Scope p(names[g]);auto r=std::transform_reduce(std::execution::par_unseq,quadrature.begin()+groupOffset[g],quadrature.begin()+groupOffset[g+1],FluxPair{},[](FluxPair const& a,FluxPair const& b){return FluxPair{a.q+b.q,a.m+b.m};},[qunit,rhoPhys](DeviceSample const& s){Vec u(0.,0.,0.);T rho=0;for(int k=0;k<8;++k){auto m=macro(s.cells[k]);u+=s.weights[k]*Vec(m.u[0],m.u[1],m.u[2]);rho+=s.weights[k]*m.rho;}T q=s.area*dot(u,s.normal)*qunit;return FluxPair{q,rhoPhys*rho*q};});flux[g]=r.q;mflux[g]=r.m;}
    }
};
GPUState* activeGPU=nullptr;
