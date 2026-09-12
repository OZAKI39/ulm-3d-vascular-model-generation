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
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
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
struct FluidNode { uint64_t index; LCell* cell; };
struct Sample { int port; T area; Vec normal; std::array<LCell*,8> cells; std::array<T,8> weights; };
struct Stats { T rmin=2,rmax=0,rsum=0,umax=0; uint64_t nan=0,inf=0; };
Stats inspect(std::vector<FluidNode> const& nodes) {
    Stats s;
    for(auto const& n:nodes) {
        T r=n.cell->computeDensity(); Vec u; n.cell->computeVelocity(u); T speed=norm(u);
        if(std::isnan(r)||std::isnan(speed)) ++s.nan;
        if(std::isinf(r)||std::isinf(speed)) ++s.inf;
        s.rmin=std::min(s.rmin,r); s.rmax=std::max(s.rmax,r); s.rsum+=r; s.umax=std::max(s.umax,speed);
    }
    return s;
}
void snapshot(std::string const& run,int step,std::vector<FluidNode> const& nodes) {
    std::ofstream f(run+"/diagnostics/fields_"+std::to_string(step)+".bin",std::ios::binary);
    uint64_t count=nodes.size(); f.write(reinterpret_cast<char*>(&count),8);
    for(auto const& n:nodes) {
        f.write(reinterpret_cast<char const*>(&n.index),8);
        T r=n.cell->computeDensity(); Vec u; n.cell->computeVelocity(u);
        f.write(reinterpret_cast<char*>(&r),8);
        for(int d=0;d<3;++d) f.write(reinterpret_cast<char*>(&u[d]),8);
    }
    require(bool(f),"Snapshot write failed");
}

int main(int argc,char** argv) {
    plbInit(&argc,&argv);
    std::string run; int completed=0;
    try {
        require(argc==2 && global::mpi().getSize()==1,"Usage: vascular_pure_fluid RUN_DIR; MPI1 only");
        run=argv[1]; global::directories().setOutputDir(run+"/");
        pcout<<std::setprecision(17);
        std::ifstream input(run+"/contracts/solver_parameters.txt"); std::string stl,s2;
        input>>stl>>s2;
        int nx,ny,nz,rampSteps,maxSteps; Vec expectedOrigin; T expectedDx,dt,tau,rhoPhys,qTarget,maLimit,rhoLow,rhoHigh;
        input>>nx>>ny>>nz>>expectedOrigin[0]>>expectedOrigin[1]>>expectedOrigin[2]>>expectedDx>>dt>>tau>>rhoPhys>>qTarget>>maLimit>>rhoLow>>rhoHigh>>rampSteps>>maxSteps;
        require(maxSteps<=1000&&rampSteps>0&&tau>.5,"Unsafe execution contract");
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
        std::ofstream capAudit(run+"/diagnostics/native_cap_area.csv"); capAudit<<std::setprecision(17)<<"label,projected_native_area_lu2,frozen_area_m2\n";
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
            auto const& v=voxels.getComponent(entry.first); Dot3D loc=v.getLocation(); auto b=entry.second;
            for(plint z=b.z0;z<=b.z1;++z) for(plint y=b.y0;y<=b.y1;++y) for(plint x=b.x0;x<=b.x1;++x) {
                size_t i=(z*ny+y)*nx+x; int flag=v.get(x-loc.x,y-loc.y,z-loc.z);
                require(flag==native[i],"Native voxel flag changed from Step2");
                bool fluid=flag==voxelFlag::inside||flag==voxelFlag::innerBorder;
                require(fluid==bool(closed[i]),"Closed lumen flag mismatch"); actualFluid+=fluid; ++allocated;
            }
        }
        for(size_t i=0;i<count;++i) {
            require(labels[i]<=4,"Invalid frozen port label");
            require(opened[i]==uint8_t(closed[i]||labels[i]),"Frozen opened/closed/port relation invalid");
            if(labels[i]) { require(native[i]==voxelFlag::outerBorder,"Port diagnostic node is not native boundary support"); ++labelCounts[labels[i]]; }
        }
        require(actualFluid==182694&&labelCounts[1]==194&&labelCounts[2]==85&&labelCounts[3]==93&&labelCounts[4]==125,"Step2 counts changed");
        std::ofstream geometryAudit(run+"/diagnostics/runtime_geometry_check.json"); geometryAudit<<std::setprecision(17);
        geometryAudit<<"{\"status\":\"PASS\",\"native_bulk_flags_equal_step2\":true,\"closed_lumen_equal_step2\":true,\"port_labels_unchanged\":true,\"cap_triangle_tags_exact\":true,\"physical_fluid_voxels\":"<<actualFluid<<",\"allocated_bulk_voxels\":"<<allocated<<",\"native_envelope\":2,\"ghost_nodes_are_physical_fluid\":false,\"profile_integral_error\":"<<integralError<<",\"inlet_native_projected_area_m2\":"<<inletProjectedArea*dx*dx<<",\"inlet_frozen_area_m2\":"<<unmodifiedArea<<"}\n"; geometryAudit.close();
        closed.clear();closed.shrink_to_fit();opened.clear();opened.shrink_to_fit();labels.clear();labels.shrink_to_fit();
        auto lattice=generateMultiBlockLattice<T,DESCRIPTOR>(voxels,2,new BGKdynamics<T,DESCRIPTOR>(1/tau));
        lattice->periodicity().toggleAll(false); lattice->toggleInternalStatistics(false);
        defineDynamics(*lattice,voxels,lattice->getBoundingBox(),new NoDynamics<T,DESCRIPTOR>,voxelFlag::outside);
        auto* model=new GuoOffLatticeModel3D<T,DESCRIPTOR>(new TriangleFlowShape3D<T,Vec>(domain.getBoundary(),profiles),voxelFlag::inside,true);
        model->setVelIsJ(false);model->selectSecondOrder(true);model->selectUseRegularizedModel(true);model->selectComputeStat(false);
        OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Vec> bc(model,domain,*lattice); bc.insert();
        initializeAtEquilibrium(*lattice,lattice->getBoundingBox(),1.,Vec(0.,0.,0.)); lattice->initialize();
        std::vector<FluidNode> nodes,ghosts;
        for(auto const& entry:bulks) {
            auto& block=lattice->getComponent(entry.first); Dot3D loc=block.getLocation(); auto b=entry.second;
            for(plint z=b.z0;z<=b.z1;++z) for(plint y=b.y0;y<=b.y1;++y) for(plint x=b.x0;x<=b.x1;++x) {
                size_t i=(z*ny+y)*nx+x; auto* cell=&block.get(x-loc.x,y-loc.y,z-loc.z);
                if(native[i]==3||native[i]==4) nodes.push_back({i,cell});
                else if(native[i]==2) ghosts.push_back({i,cell});
            }
        }
        require(nodes.size()==actualFluid,"Fluid cache count mismatch");
        std::vector<Sample> samples; std::ifstream quad(run+"/contracts/flux_quadrature.tsv"); int label; T x,y,z,a,n0,n1,n2;
        while(quad>>label>>x>>y>>z>>a>>n0>>n1>>n2) {
            Sample s;s.port=label-1;s.area=a;s.normal=Vec(n0,n1,n2); int b[3]={int(std::floor(x)),int(std::floor(y)),int(std::floor(z))};T frac[3]={x-b[0],y-b[1],z-b[2]};int j=0;
            for(int ox=0;ox<2;++ox) for(int oy=0;oy<2;++oy) for(int oz=0;oz<2;++oz) {
                int ix=b[0]+ox,iy=b[1]+oy,iz=b[2]+oz;
                require(ix>=0&&ix<nx&&iy>=0&&iy<ny&&iz>=0&&iz<nz,"Flux sample outside grid");
                require(native[(iz*ny+iy)*nx+ix]>=1&&native[(iz*ny+iy)*nx+ix]<=4,"Flux sample unallocated");
                s.cells[j]=&lattice->get(ix,iy,iz);
                s.weights[j]=(ox?frac[0]:1-frac[0])*(oy?frac[1]:1-frac[1])*(oz?frac[2]:1-frac[2]);++j;
            }
            samples.push_back(s);
        }
        require(quad.eof()&&!samples.empty(),"Incomplete flux quadrature");native.clear();native.shrink_to_fit();
        T initialMass=inspect(nodes).rsum*massUnit;
        std::ofstream history(run+"/diagnostics/flow_history.csv"); history<<std::setprecision(17);
        history<<"iteration,time_s,rho_min,rho_max,rho_mean,u_max_m_s,u_max_lu,mach_max,q_in_m3_s,q_outlet_01_m3_s,q_outlet_02_m3_s,q_outlet_03_m3_s,flow_closure,inlet_target_error,total_mass,relative_mass_drift,inlet_outward_flux_m3_s,net_outward_flux_m3_s,startup_ramp,ghost_rho_min,ghost_rho_max,ghost_mach_max\n";
        std::ofstream stages(run+"/diagnostics/stages.csv"); stages<<"iteration,safety_status\n";
        for(int step=0;step<=maxSteps;++step) {
            if(step) {state.ramp=std::min(T(1),T(step)/rampSteps);lattice->collideAndStream();completed=step;}
            Stats s=inspect(nodes),g=inspect(ghosts);
            bool safe=!s.nan&&!s.inf&&!g.nan&&!g.inf&&s.rmin>=rhoLow&&s.rmax<=rhoHigh&&g.rmin>0&&s.umax/std::sqrt(DESCRIPTOR<T>::cs2)<=maLimit&&g.umax/std::sqrt(DESCRIPTOR<T>::cs2)<=maLimit;
            if(step==0||step==1||step%10==0||!safe) {
                T flux[4]={0,0,0,0};
                for(auto const& sample:samples) {
                    Vec u(0.,0.,0.);
                    for(int j=0;j<8;++j) {Vec v;sample.cells[j]->computeVelocity(v);u+=sample.weights[j]*v;}
                    flux[sample.port]+=sample.area*dot(u,sample.normal)*qUnit;
                }
                T qi=-flux[0], net=flux[0]+flux[1]+flux[2]+flux[3], closure=std::abs(net)/std::max(std::abs(qi),qTarget*1e-12),error=std::abs(qi-qTarget)/qTarget,mass=s.rsum*massUnit;
                history<<step<<','<<step*dt<<','<<s.rmin<<','<<s.rmax<<','<<s.rsum/nodes.size()<<','<<s.umax*velUnit<<','<<s.umax<<','<<s.umax/std::sqrt(DESCRIPTOR<T>::cs2)<<','<<qi<<','<<flux[1]<<','<<flux[2]<<','<<flux[3]<<','<<closure<<','<<error<<','<<mass<<','<<(mass-initialMass)/initialMass<<','<<flux[0]<<','<<net<<','<<state.ramp<<','<<g.rmin<<','<<g.rmax<<','<<g.umax/std::sqrt(DESCRIPTOR<T>::cs2)<<'\n';history.flush();
                if(step==0||step==1||step==10||step==100||step==1000||!safe) {
                    pcout<<"STAGE "<<step<<" safe="<<safe<<" rho=["<<s.rmin<<','<<s.rmax<<"] Ma="<<s.umax/std::sqrt(DESCRIPTOR<T>::cs2)<<" Qin="<<qi<<" inlet_error="<<error<<" closure="<<closure<<std::endl;
                    stages<<step<<','<<(safe?"PASS":"FAIL")<<'\n';stages.flush();
                }
            }
            if(!safe) { pcout<<"SAFETY_STOP NaN="<<s.nan+g.nan<<" Inf="<<s.inf+g.inf<<std::endl; throw std::runtime_error("Immediate runtime safety gate failure"); }
            if(step==0||step==10||step==100||step==1000) snapshot(run,step,nodes);
            if(step==0) {
                for(int k=0;k<5;++k) require(state.calls[k]>0,"Stage A: a surface profile was never used by native completion");
            }
        }
        std::ofstream status(run+"/diagnostics/solver_status.json"); status<<"{\"runtime_safety\":\"PASS\",\"timesteps\":"<<completed<<",\"cell_count\":0,\"mpi_ranks\":1,\"native_profile_calls\":[";
        for(int k=0;k<5;++k) status<<(k?",":"")<<state.calls[k];
        status<<"],\"note\":\"Flow acceptance is evaluated independently from measured CSV, not implied by runtime success\"}\n";
        pcout<<"PURE_FLUID_COMPLETE timesteps="<<completed<<" CELL_COUNT=0"<<std::endl;return 0;
    } catch(std::exception const& e) {
        pcout<<"STEP3_ERROR completed_timesteps="<<completed<<" "<<e.what()<<std::endl;
        if(!run.empty()) {std::ofstream f(run+"/diagnostics/solver_failure.txt");f<<"completed_timesteps="<<completed<<"\n"<<e.what()<<'\n';}
        return 2;
    }
}
