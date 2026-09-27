// Application-level MPI mechanics only. Native physics and library code are unchanged.
#include <mpi.h>
#include <cstring>
#include <streambuf>
int rankId(){return global::mpi().getRank();}
int rankCount(){return global::mpi().getSize();}
struct DiscardBuffer:std::streambuf {int overflow(int c) override{return traits_type::not_eof(c);}std::streamsize xsputn(const char*,std::streamsize n) override{return n;}};
class RankFile:public std::ostream {
    std::ofstream file;DiscardBuffer discard;
public:
    RankFile(std::string const& path,std::ios::openmode mode=std::ios::out):std::ostream(nullptr){
        if(rankId()==0){file.open(path,mode);rdbuf(file.rdbuf());if(!file)setstate(std::ios::badbit);}
        else rdbuf(&discard);
    }
    void close(){if(rankId()==0)file.close();}
};
uint64_t sumCount(uint64_t value){uint64_t result=0;MPI_Allreduce(&value,&result,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);return result;}
T maxValue(T value){T result=0;MPI_Allreduce(&value,&result,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);return result;}
Stats inspect(std::vector<FluidNode> const& nodes) {
#ifdef GPU_POC
    if(activeGPU)return activeGPU->inspect(!nodes.empty()&& !activeGPU->ghost.empty() && nodes.front().index==activeGPU->ghost.front().id);
#endif
    Stats s;
    for(auto const& n:nodes){T r=n.cell->computeDensity();Vec u;n.cell->computeVelocity(u);T speed=norm(u);
        if(std::isnan(r)||std::isnan(speed))++s.nan;
        if(std::isinf(r)||std::isinf(speed))++s.inf;
        s.rmin=std::min(s.rmin,r);s.rmax=std::max(s.rmax,r);s.rsum+=r;if(n.inCV)s.cvsum+=r;s.speedSum+=speed;s.umax=std::max(s.umax,speed);
    }
    Stats g;T localSums[3]={s.rsum,s.cvsum,s.speedSum},globalSums[3];T localMax[2]={s.rmax,s.umax},globalMax[2];uint64_t flags[2]={s.nan,s.inf},globalFlags[2];
    MPI_Allreduce(&s.rmin,&g.rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(localMax,globalMax,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(localSums,globalSums,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(flags,globalFlags,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    g.rmax=globalMax[0];g.umax=globalMax[1];g.rsum=globalSums[0];g.cvsum=globalSums[1];g.speedSum=globalSums[2];g.nan=globalFlags[0];g.inf=globalFlags[1];return g;
}
struct PackedNode {uint64_t index;T rho;T u[3];};
static_assert(sizeof(PackedNode)==40,"Frozen snapshot serialization requires 40 bytes per cell");
std::vector<PackedNode> gatherNodes(std::vector<FluidNode> const& nodes){
#ifdef GPU_POC
    if(activeGPU){std::vector<uint64_t> ids;for(auto const& n:nodes)ids.push_back(n.index);auto v=activeGPU->pack(ids);std::vector<PackedNode> out(v.size());std::memcpy(out.data(),v.data(),out.size()*sizeof(PackedNode));std::sort(out.begin(),out.end(),[](PackedNode const&a,PackedNode const&b){return a.index<b.index;});return out;}
#endif
    std::vector<PackedNode> local(nodes.size());
    for(size_t i=0;i<nodes.size();++i){local[i].index=nodes[i].index;local[i].rho=nodes[i].cell->computeDensity();Vec u;nodes[i].cell->computeVelocity(u);for(int d=0;d<3;++d)local[i].u[d]=u[d];}
    require(local.size()*sizeof(PackedNode)<size_t(std::numeric_limits<int>::max()),"MPI field count overflow");
    int count=int(local.size()*sizeof(PackedNode));std::vector<int> counts(rankCount()),offsets(rankCount());
    MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);int total=0;
    if(rankId()==0)for(int i=0;i<rankCount();++i){offsets[i]=total;total+=counts[i];}
    std::vector<PackedNode> global(size_t(total)/sizeof(PackedNode));
    MPI_Gatherv(local.data(),count,MPI_BYTE,global.data(),counts.data(),offsets.data(),MPI_BYTE,0,MPI_COMM_WORLD);
    if(rankId()==0){std::sort(global.begin(),global.end(),[](PackedNode const&a,PackedNode const&b){return a.index<b.index;});for(size_t i=1;i<global.size();++i)require(global[i-1].index<global[i].index,"FAIL_MPI_RUNTIME_CORRECTNESS: duplicate global field owner");}
    return global;
}
void saveGathered(std::string const& path,std::vector<FluidNode> const& nodes){
    auto global=gatherNodes(nodes);
    if(rankId()==0){std::ofstream f(path,std::ios::binary);uint64_t count=global.size();f.write(reinterpret_cast<char*>(&count),8);f.write(reinterpret_cast<char*>(global.data()),global.size()*sizeof(PackedNode));require(bool(f),"Global snapshot write failure");}
    MPI_Barrier(MPI_COMM_WORLD);
}
void snapshot(std::string const& run,int step,std::vector<FluidNode> const& nodes){saveGathered(run+"/diagnostics/field_samples/fields_"+std::to_string(step)+".bin",nodes);}
void sampledSnapshot(std::string const& run,int step,std::map<uint64_t,LCell*> const& samples){
    std::vector<FluidNode> nodes;for(auto const&e:samples)nodes.push_back(FluidNode(e.first,e.second));
    saveGathered(run+"/diagnostics/field_samples/samples_"+std::to_string(step)+".bin",nodes);
}
std::vector<T> sampledMacros(std::vector<uint64_t> const& ids,std::map<uint64_t,LCell*> const& owned){
    std::vector<T> local(4*ids.size(),0),global(local.size());
    for(auto const&e:owned){size_t slot=std::lower_bound(ids.begin(),ids.end(),e.first)-ids.begin();require(slot<ids.size()&&ids[slot]==e.first,"Sample coordinate missing");local[4*slot]=e.second->computeDensity();Vec u;e.second->computeVelocity(u);for(int d=0;d<3;++d)local[4*slot+d+1]=u[d];}
    MPI_Allreduce(local.data(),global.data(),int(global.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(T value:global)require(std::isfinite(value),"FAIL_MPI_RUNTIME_CORRECTNESS: nonfinite reduced sample macro");return global;
}
uint64_t ownershipChecks=0;
void checkQuadratureOwnership(std::vector<Sample> const& samples,std::string const& run,bool emit){
    std::vector<int> local(samples.size()),global(samples.size());
    for(size_t i=0;i<samples.size();++i)local[i]=samples[i].owner==rankId();
    MPI_Allreduce(local.data(),global.data(),int(global.size()),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    uint64_t duplicate=0,unowned=0;for(int count:global){duplicate+=count>1;unowned+=count==0;}
    require(!duplicate&&!unowned,"FAIL_MPI_RUNTIME_CORRECTNESS: quadrature ownership invariant");++ownershipChecks;
    if(emit){RankFile f(run+"/diagnostics/quadrature_ownership_check.json");f<<"{\"status\":\"PASS\",\"quadrature_points\":"<<samples.size()<<",\"duplicate_owned_points\":0,\"unowned_points\":0,\"owner_definition\":\"rank owning the floor-coordinate bulk cell\",\"runtime_recheck_interval_steps\":100}\n";}
}
void mpiFieldResidual(std::string const& run,int step,int previous,std::vector<FluidNode> const& nodes,T dx,T dt,T rhoPhys){
    std::ifstream f(run+"/diagnostics/field_samples/fields_"+std::to_string(previous)+".bin",std::ios::binary);uint64_t count=0;f.read(reinterpret_cast<char*>(&count),8);require(count==182694,"MPI residual previous field size");
    std::vector<PackedNode> old(count);f.read(reinterpret_cast<char*>(old.data()),count*sizeof(PackedNode));require(bool(f),"MPI residual previous field read");
    require(std::is_sorted(old.begin(),old.end(),[](PackedNode const&a,PackedNode const&b){return a.index<b.index;}),"MPI residual previous global field order");
    T local[4]={0,0,0,0},global[4],scale=dx/dt,punit=rhoPhys*scale*scale/3;
    for(auto const& node:nodes){auto it=std::lower_bound(old.begin(),old.end(),node.index,[](PackedNode const&a,uint64_t id){return a.index<id;});require(it!=old.end()&&it->index==node.index,"MPI residual ownership index");Vec u;node.cell->computeVelocity(u);
        for(int d=0;d<3;++d){T now=u[d]*scale,before=it->u[d]*scale;local[0]+=(now-before)*(now-before);local[1]+=now*now;}
        T now=(node.cell->computeDensity()-1)*punit,before=(it->rho-1)*punit;local[2]+=(now-before)*(now-before);local[3]+=now*now;
    }
    MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(T value:global)require(std::isfinite(value)&&value>=0,"FAIL_MPI_RUNTIME_CORRECTNESS: residual sum invalid");
    RankFile out(run+"/diagnostics/mpi_field_residual_"+std::to_string(step)+".json");out<<std::setprecision(17)<<"{\"iteration\":"<<step<<",\"previous_iteration\":"<<previous<<",\"R_velocity\":"<<std::sqrt(global[0]/std::max(global[1],count*1e-24))<<",\"R_pressure\":"<<std::sqrt(global[2]/std::max(global[3],T(count)))<<",\"global_u_difference_squared_sum\":"<<global[0]<<",\"global_u_squared_sum\":"<<global[1]<<",\"global_p_difference_squared_sum\":"<<global[2]<<",\"global_p_squared_sum\":"<<global[3]<<",\"method\":\"local numerator/denominator SUM then MPI_SUM before ratio\"}\n";
}
void evaluate(std::string const& root,std::string const& run,int step){
    int ok=1;
    if(rankId()==0){std::string script=root+"/scripts/monitor_online.py",number=std::to_string(step);pid_t pid=fork();require(pid>=0,"Cannot fork read-only evaluator");if(pid==0){execl("/usr/bin/python3","/usr/bin/python3",script.c_str(),root.c_str(),run.c_str(),number.c_str(),(char*)nullptr);_exit(127);}int status=0;ok=waitpid(pid,&status,0)==pid&&WIFEXITED(status)&&WEXITSTATUS(status)==0;}
    MPI_Bcast(&ok,1,MPI_INT,0,MPI_COMM_WORLD);require(ok,"Convergence evaluator failed; stop without changing BC");
}
