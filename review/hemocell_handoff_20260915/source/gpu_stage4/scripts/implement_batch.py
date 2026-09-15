from pathlib import Path
import hashlib,json
R=Path.cwd()
def put(p,t):
 q=p.with_name(p.name+'.stage4.tmp');q.write_text(t);q.replace(p)
native=R/'upstream/palabos/src';base=R/'baseline/native/palabos/src'
hh=(base/'atomicBlock/atomicAcceleratedLattice3D.hh').read_text()
start=hh.index('void AtomicAcceleratedLattice3D<T, Descriptor>::collideAndStream(CollFun const &collFun)')
tail=hh[start:hh.index('template <typename T, template <typename U> class Descriptor>',start+5)]
body=tail.split('            size_t i = &f0 - populations;\n',1)[1].split('\n        });',1)[0]
h=(base/'atomicBlock/atomicAcceleratedLattice3D.h').read_text()
marker='    template <class CollFun>\n    void collideAndStream(CollFun const &collFun);'
assert h.count(marker)==1
decl="""
    // Stage4 scratch-only immutable batch view. No object layout/vtable changes.
    struct Stage4BatchView {
        static inline bool readBit(const uint32_t &mask, uint8_t position)
        {
            return (mask >> position) & 1ULL;
        }
        plint nx, ny, nz, N;
        Box3D fullDomain;
        T* buffers[2];
        ExternalFieldArray<T, typename Descriptor<T>::ExternalField>* externalScalars;
        int* collisionMatrix;
        uint32_t* hw_bb_links;
        plint* dynamicScalarIndex;
        T* dynamicScalarsPtr;
        bool hasDynamicScalars;
        int backgroundCollisionModel;
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars;
        template <class CollFun>
        void collideCell(size_t i, CollFun const& collFun, int parity) const;
    };
    Stage4BatchView stage4BatchView();
    void stage4SwapBuffers();
"""
put(native/'atomicBlock/atomicAcceleratedLattice3D.h',h.replace(marker,marker+decl))
definitions="""
template <typename T, template <typename U> class Descriptor>
typename AtomicAcceleratedLattice3D<T, Descriptor>::Stage4BatchView
AtomicAcceleratedLattice3D<T, Descriptor>::stage4BatchView()
{
    Stage4BatchView v;
    v.nx=this->getNx();v.ny=this->getNy();v.nz=this->getNz();v.N=N;
    v.fullDomain=this->getBoundingBox();
    v.buffers[0]=populations;v.buffers[1]=tmpPopulations;
    v.externalScalars=externalScalars;v.collisionMatrix=collisionMatrix;
    v.hw_bb_links=hw_bb_links;v.dynamicScalarIndex=dynamicScalarIndex;
    v.hasDynamicScalars=!dynamicScalars.empty();v.dynamicScalarsPtr=dynamicScalars.data();
    v.backgroundCollisionModel=toCollisionModel(this->getBackgroundDynamics());
    getStaticScalars(this->getBackgroundDynamics(),v.backgroundCollisionModel,v.staticScalars);
    return v;
}
template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::stage4SwapBuffers()
{
    std::swap(populations,tmpPopulations);
    std::swap(populationGrid,tmpPopulationGrid);
}
template <typename T, template <typename U> class Descriptor>
template <class CollFun>
void AtomicAcceleratedLattice3D<T, Descriptor>::Stage4BatchView::collideCell(
    size_t i, CollFun const& collFun, int parity) const
{
    T* populations=buffers[parity];
    T* tmpPopulations=buffers[1-parity];
"""+body+"""
}
"""
assert body in definitions
put(native/'atomicBlock/atomicAcceleratedLattice3D.hh',hh.replace('\n}  // namespace plb','\n'+definitions+'\n}  // namespace plb'))
assert definitions in (native/'atomicBlock/atomicAcceleratedLattice3D.hh').read_text()
mh=(base/'multiBlock/acceleratedLattice3D.h').read_text()
put(native/'multiBlock/acceleratedLattice3D.h',mh.replace(marker,marker+'\n    template <class CollFun, class BatchPlan>\n    void stage4CollideAndStream(CollFun const& collFun, BatchPlan& plan);'))
mhh=(base/'multiBlock/acceleratedLattice3D.hh').read_text()
a=mhh.index('template <typename T, template <typename U> class Descriptor>\ntemplate <class CollFun>\nvoid AcceleratedLattice3D<T, Descriptor>::collideAndStream(CollFun const &collFun)')
b=mhh.index('\ntemplate <typename T, template <typename U> class Descriptor>',a+10)
wrapper=mhh[a:b].replace('template <class CollFun>','template <class CollFun, class BatchPlan>').replace('::collideAndStream(CollFun const &collFun)','::stage4CollideAndStream(CollFun const &collFun, BatchPlan& plan)').replace('collideAndStreamImplementation<CollFun>(collFun);','plan.apply(collFun);')
put(native/'multiBlock/acceleratedLattice3D.hh',mhh[:b]+'\n'+wrapper+'\n'+mhh[b:])
plan="""// Stage4 scratch: batching only. Original collision/streaming body lives in
// the scratch native view, mechanically copied by scripts/implement_batch.py.
#pragma once
struct Stage4LBMBatch {
    struct Job { uint32_t block,cell; };
    std::vector<AL::Stage4BatchView> views;
    std::vector<Job> jobs;
    std::vector<AL*> batched,independent;
    int parity=0;
    void initialize(ML& lattice,std::string const& run,bool enabled) {
        static_assert(sizeof(AL)==504,"Stage3 atomic ABI changed");
        std::ofstream f(run+"/diagnostics/lbm_runtime_compatibility.csv");
        f<<"block_id,N,nx,ny,nz,background_collision,collision_models,external_scalars,dynamic_scalars,batch_compatible,selected_path\\n";
        for(auto const& entry:lattice.getBlockLattices()) {
            auto* a=entry.second;auto v=a->stage4BatchView();
            std::set<int> models;for(plint i=0;i<v.N;++i)models.insert(v.collisionMatrix[i]);
            bool compatible=DD::ExternalField::numScalars==0 &&
                v.backgroundCollisionModel==CollisionModel::BGK && !v.hasDynamicScalars;
            for(int m:models)compatible=compatible&&(m==CollisionModel::BGK||m==CollisionModel::NoDynamics);
            f<<entry.first<<','<<v.N<<','<<v.nx<<','<<v.ny<<','<<v.nz<<','<<v.backgroundCollisionModel<<',';
            for(int m:models)f<<m<<';';
            f<<','<<DD::ExternalField::numScalars<<','<<v.hasDynamicScalars<<','<<compatible<<','<<(enabled&&compatible?"BATCH":"ORIGINAL")<<'\\n';
            if(enabled&&compatible) {
                uint32_t id=views.size();views.push_back(v);batched.push_back(a);
                for(plint i=0;i<v.N;++i)jobs.push_back({id,uint32_t(i)});
            } else if(enabled) independent.push_back(a);
        }
        require(bool(f),"LBM compatibility inventory write failure");
        if(!enabled)return;
        require(!jobs.empty(),"No safely compatible LBM batch");
        // Advice/prefetch only once, after host initialization of persistent arrays.
        cudaMemLocation dev{};dev.type=cudaMemLocationTypeDevice;dev.id=0;
        auto persist=[&](void* ptr,size_t bytes){
            require(cudaMemAdvise(ptr,bytes,cudaMemAdviseSetPreferredLocation,dev)==cudaSuccess,"Batch metadata preferred device failed");
            require(cudaMemPrefetchAsync(ptr,bytes,dev,0,0)==cudaSuccess,"Batch metadata prefetch failed");
        };
        persist(views.data(),views.size()*sizeof(views[0]));
        persist(jobs.data(),jobs.size()*sizeof(jobs[0]));
        require(cudaDeviceSynchronize()==cudaSuccess,"Initial batch metadata prefetch failed");
        std::ofstream meta(run+"/diagnostics/batch_metadata.json");
        meta<<"{\\"BATCH_METADATA_CREATED_PER_STEP\\":false,\\"descriptor_count\\":"<<views.size()<<",\\"job_count\\":"<<jobs.size()
            <<",\\"metadata_bytes\\":"<<views.size()*sizeof(views[0])+jobs.size()*sizeof(jobs[0])
            <<",\\"fallback_blocks\\":"<<independent.size()<<",\\"atomic_ABI_bytes\\":"<<sizeof(AL)<<"}\\n";
        require(bool(meta),"Batch metadata evidence failure");
    }
    template<class CollFun> void apply(CollFun const& collFun) {
        auto const* descriptor=views.data();int phase=parity;
        std::for_each(std::execution::par_unseq,jobs.begin(),jobs.end(),
            [descriptor,phase,collFun](Job const& job){
                descriptor[job.block].collideCell(size_t(job.cell),collFun,phase);
            });
        for(AL* a:batched)a->stage4SwapBuffers();
        for(AL* a:independent)a->collideAndStream(collFun);
        parity=1-parity;
    }
};
"""
put(R/'source/stage4_lbm_batch.hpp',plan)
g=(R/'baseline/source/gpu_adapter.hpp').read_text()
g=g.replace('struct GPUState {','#include "stage4_lbm_batch.hpp"\nstruct GPUState {\n    Stage4LBMBatch lbmBatch;bool useBatch=false;')
g=g.replace('    void sync(){','    void sync(){',1)
needle='    }\n    void sync(){'
assert g.count(needle)==1
g=g.replace(needle,'''        const char* variant=std::getenv("STAGE4_LBM");
        require(variant&&(std::string(variant)=="baseline"||std::string(variant)=="batch"),"Explicit Stage4 LBM variant required");
        useBatch=std::string(variant)=="batch";
        lbmBatch.initialize(*lattice,run,useBatch);
    }
    void sync(){''')
old='lattice->collideAndStream(CollisionKernel<T,DESCRIPTOR,CollisionModel::BGK,CollisionModel::NoDynamics>());'
new='if(useBatch)lattice->stage4CollideAndStream(CollisionKernel<T,DESCRIPTOR,CollisionModel::BGK,CollisionModel::NoDynamics>(),lbmBatch);else '+old
assert g.count(old)==1;g=g.replace(old,new);put(R/'source/gpu_adapter.hpp',g)
c=(R/'baseline/source/vascularPoC.cpp').read_text().replace('GPU_POC_STAGE_3','GPU_POC_STAGE_4').replace('||requested==10000','').replace('Stage3 frozen horizon','Stage4 frozen horizon').replace('maxSteps<=10000,"Stage3 hard cap 10000"','maxSteps<=5000,"Stage4 hard cap 5000"').replace('std::set<int> milestones={0,100,200,1000,5000,10000};','std::set<int> milestones={0,100,200,1000,5000};\n        if(std::getenv("STAGE4_STRONG_CHECKPOINTS")){milestones.insert(1);milestones.insert(10);}')
put(R/'source/vascularPoC.cpp',c)
proof=dict(status='PASS',body_extraction_sha256=hashlib.sha256(body.encode()).hexdigest(),original_per_cell_body_verbatim=True,wrapper_changes='dispatch call only plus method signature',object_layout_changed=False,halo_guo_monitor_bodies_changed=False,created_before_first_compile=True)
(R/'provenance/STATIC_MATH_REUSE_PROOF.json').write_text(json.dumps(proof,indent=2))
print('IMPLEMENTED',proof)