// Adaptation for the bounded RBC repair campaign. No change to bounce physics.
// Included only by an independently built, explicitly selected Mirheo library.
#pragma once
#include <mirheo/core/mirheo_state.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>
#include <vector>

namespace mirheo {
namespace rbc_repair {

template<class T>
std::vector<T> download(const T* ptr, int n, cudaStream_t stream)
{
    std::vector<T> host(n);
    if (n > 0) CUDA_Check(cudaMemcpyAsync(host.data(), ptr, sizeof(T)*n, cudaMemcpyDeviceToHost, stream));
    CUDA_Check(cudaStreamSynchronize(stream));
    return host;
}

inline bool enabled()
{
    const char* value = std::getenv("RBC_REPAIR_NATIVE_TRACE");
    return value && std::string(value) == "1";
}

inline void dumpPV(const std::string& filename, const real4* x, const real4* old,
                   const real4* vel, const real4* force, int n, cudaStream_t stream)
{
    auto xx=download(x,n,stream), oo=download(old,n,stream);
    auto vv=download(vel,n,stream), ff=download(force,n,stream);
    std::ofstream out(filename);
    if (!out) die("Cannot open repair evidence %s", filename.c_str());
    out << std::setprecision(17) << "local_index,id,x,y,z,old_x,old_y,old_z,vx,vy,vz,fx,fy,fz\n";
    for (int i=0;i<n;++i) {
        Particle particle(xx[i],vv[i]);
        out << i << ',' << particle.getId() << ',' << xx[i].x << ',' << xx[i].y << ',' << xx[i].z
            << ',' << oo[i].x << ',' << oo[i].y << ',' << oo[i].z
            << ',' << vv[i].x << ',' << vv[i].y << ',' << vv[i].z
            << ',' << ff[i].x << ',' << ff[i].y << ',' << ff[i].z << '\n';
    }
}

// Called before resolution. Existing coarse/fine downloads have completed.
// Full point states are saved only at step 0 or on an overflow. No trajectory
// restart, candidate deletion, capacity increase or membership correction.
inline void trace(const MirState* state, const std::string& pvName, const std::string& locality,
                  int nObjects, const OVviewWithNewOldVertices& ov, const PVviewWithOldParticles& pv,
                  int coarse, int coarseCap, const int2* coarsePairs,
                  int fine, int fineCap, const int2* finePairs, cudaStream_t stream)
{
    if (!enabled()) return;
    const bool overflow = coarse > coarseCap || fine > fineCap;
    if (!overflow && state->currentStep > 1 && state->currentStep % 250 != 0) return;
    const char* stageEnv = std::getenv("RBC_REPAIR_PHASE");
    const std::string stage = stageEnv ? stageEnv : "unknown";
    const std::string prefix = "bounce_" + stage + "_" + pvName + "_" + locality + "_" + std::to_string(state->currentStep);
    const int nc = std::min(coarse,coarseCap), nf = std::max(0,std::min(fine,fineCap));
    const auto cc = download(coarsePairs,nc,stream);
    const auto fc = download(finePairs,nf,stream);
    std::set<std::pair<int,int>> unique;
    // The full pair files retain exact native indices even for non-sphere diagnostics.
    for (auto named : {std::make_pair("coarse",&cc),std::make_pair("fine",&fc)}) {
        std::ofstream out(prefix+"_"+named.first+".csv");
        if (!out) die("Cannot write repair pair evidence");
        out << "particle_local_index,global_triangle_index\n";
        for (const auto& pair : *named.second) out << pair.x << ',' << pair.y << '\n';
    }
    for (const auto& pair : cc) unique.insert({pair.x,pair.y});
    std::ofstream out("bounce_counts_"+stage+".csv",std::ios::app);
    if (!out) die("Cannot write repair counter evidence");
    if (out.tellp()==0) out << "step,time,dt,pv,locality,n_objects,n_particles,coarse_count,coarse_capacity,fine_count,fine_capacity,stored_coarse,stored_fine,stored_coarse_duplicates\n";
    out << std::setprecision(17) << state->currentStep << ',' << state->currentTime << ',' << state->getDt()
        << ',' << pvName << ',' << locality << ',' << nObjects << ',' << pv.size
        << ',' << coarse << ',' << coarseCap << ',' << fine << ',' << fineCap
        << ',' << nc << ',' << nf << ',' << nc-unique.size() << '\n';
    if (overflow || state->currentStep==0) {
        dumpPV(prefix+"_fluid.csv",pv.positions,pv.oldPositions,pv.velocities,pv.forces,pv.size,stream);
        dumpPV(prefix+"_membrane.csv",ov.vertices,ov.old_vertices,ov.velocities,ov.vertexForces,ov.nvertices*ov.nObjects,stream);
    }
}
} // namespace rbc_repair
} // namespace mirheo
