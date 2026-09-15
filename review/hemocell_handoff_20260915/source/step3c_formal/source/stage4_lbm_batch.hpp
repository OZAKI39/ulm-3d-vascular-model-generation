// Stage4 scratch: batching only. Original collision/streaming body lives in
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
        f<<"block_id,N,nx,ny,nz,background_collision,collision_models,external_scalars,dynamic_scalars,batch_compatible,selected_path\n";
        for(auto const& entry:lattice.getBlockLattices()) {
            auto* a=entry.second;auto v=a->stage4BatchView();
            std::set<int> models;for(plint i=0;i<v.N;++i)models.insert(v.collisionMatrix[i]);
            bool compatible=DD::ExternalField::numScalars==0 &&
                v.backgroundCollisionModel==CollisionModel::BGK && !v.hasDynamicScalars;
            for(int m:models)compatible=compatible&&(m==CollisionModel::BGK||m==CollisionModel::NoDynamics);
            f<<entry.first<<','<<v.N<<','<<v.nx<<','<<v.ny<<','<<v.nz<<','<<v.backgroundCollisionModel<<',';
            for(int m:models)f<<m<<';';
            f<<','<<DD::ExternalField::numScalars<<','<<v.hasDynamicScalars<<','<<compatible<<','<<(enabled&&compatible?"BATCH":"ORIGINAL")<<'\n';
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
        meta<<"{\"BATCH_METADATA_CREATED_PER_STEP\":false,\"descriptor_count\":"<<views.size()<<",\"job_count\":"<<jobs.size()
            <<",\"metadata_bytes\":"<<views.size()*sizeof(views[0])+jobs.size()*sizeof(jobs[0])
            <<",\"fallback_blocks\":"<<independent.size()<<",\"atomic_ABI_bytes\":"<<sizeof(AL)<<"}\n";
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
