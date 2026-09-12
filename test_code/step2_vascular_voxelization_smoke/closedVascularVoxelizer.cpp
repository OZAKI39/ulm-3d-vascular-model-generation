// Test-local geometry-only path, based on helper/voxelizeDomain.cpp's
// TriangleSet -> DEFscaledMesh -> TriangleBoundary -> VoxelizedDomain stages.
// Existing HemoCell/Palabos source remains unchanged. No fluid lattice exists.
#include "palabos3D.h"
#include "palabos3D.hh"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <cstdint>
using namespace plb;

template<class T> void jsonVec(std::ostream& s, Array<T,3> const& v) {
    s << '[' << v[0] << ',' << v[1] << ',' << v[2] << ']';
}
int main(int argc, char** argv) {
    plbInit(&argc, &argv);
    try {
        if (argc != 5 || global::mpi().getSize() != 1)
            throw std::runtime_error("Usage: closed_voxelizer STL refDir refDirN RUN_DIR; MPI 1 only");
        const std::string stl=argv[1], out=argv[4];
        const plint refDir=std::stoi(argv[2]), refN=std::stoi(argv[3]);
        const plint margin=1, extraLayer=0, borderWidth=1, envelope=1, blockSize=16;
        global::directories().setOutputDir(out+"/");
        TriangleSet<double> triangles(stl, DBL);
        DEFscaledMesh<double> scaled(triangles, refN, refDir, margin, extraLayer);
        const auto origin=scaled.getPhysicalLocation();
        const double dx=scaled.getDx();
        Array<double,2> xr,yr,zr;
        scaled.getMesh().computeBoundingBox(xr,yr,zr);
        std::ofstream vertices(out+"/diagnostics/palabos_vertices_preinflate.f64",std::ios::binary);
        for (auto const& v:scaled.getVertexList())
            for (int j=0;j<3;++j) { double x=v[j]; vertices.write(reinterpret_cast<char*>(&x),sizeof x); }
        vertices.close();
        TriangleBoundary3D<double> boundary(scaled,false);
        boundary.getMesh().inflate(); // official default: 0.001 LU, recorded explicitly
        Array<double,2> ix,iy,iz;
        boundary.getMesh().computeBoundingBox(ix,iy,iz);
        pcout << std::setprecision(17) << "GEOMETRY_ONLY refDir=" << refDir << " refDirN=" << refN
              << " dx=" << dx << " margin=1 extraLayer=0 inflate_lu=0.001" << std::endl;
        VoxelizedDomain3D<double> domain(boundary,voxelFlag::inside,extraLayer,borderWidth,envelope,blockSize);
        auto& voxels=domain.getVoxelMatrix();
        MultiScalarField3D<int> flags(static_cast<MultiBlock3D&>(voxels));
        setToConstant(flags,flags.getBoundingBox(),0);
        setToConstant(flags,voxels,voxelFlag::inside,flags.getBoundingBox(),1);
        setToConstant(flags,voxels,voxelFlag::innerBorder,flags.getBoundingBox(),1);
        // STOP here. The official helper's CopyFromNeighbor X-end opening is omitted.
        const Box3D b=flags.getBoundingBox();
        const size_t nx=b.getNx(),ny=b.getNy(),nz=b.getNz(),n=nx*ny*nz;
        std::vector<uint8_t> dense(n,0),native(n,255);
        uint64_t counts[6]={0},fluid=0;
        auto const& bulks=flags.getMultiBlockManagement().getSparseBlockStructure().getBulks();
        for (auto const& entry:bulks) {
            auto const& f=flags.getComponent(entry.first);
            auto const& v=voxels.getComponent(entry.first);
            Dot3D loc=f.getLocation(); auto bb=entry.second;
            for(plint z=bb.z0;z<=bb.z1;++z) for(plint y=bb.y0;y<=bb.y1;++y) for(plint x=bb.x0;x<=bb.x1;++x) {
                size_t i=((z-b.z0)*ny+(y-b.y0))*nx+(x-b.x0);
                int flag=f.get(x-loc.x,y-loc.y,z-loc.z), code=v.get(x-loc.x,y-loc.y,z-loc.z);
                if(code<0 || code>5 || (flag!=0 && flag!=1)) throw std::runtime_error("Invalid native flag");
                dense[i]=flag; native[i]=code; counts[code]++; fluid+=flag;
            }
        }
        for(auto const& item:std::vector<std::pair<std::string,std::vector<uint8_t>*>>{
                {"closed_flag_matrix.u8",&dense},{"palabos_native_flags.u8",&native}}) {
            std::ofstream raw(out+"/diagnostics/"+item.first,std::ios::binary);
            raw.write(reinterpret_cast<char*>(item.second->data()),n);
            if(!raw) throw std::runtime_error("Raw output failure");
        }
        std::ofstream meta(out+"/diagnostics/palabos_geometry.json"); meta<<std::setprecision(17);
        meta<<"{\n\"effective_dx_m\":"<<dx<<",\n\"physical_origin_m\":"; jsonVec(meta,origin);
        meta<<",\n\"ref_dir\":"<<refDir<<",\"ref_dir_n\":"<<refN
            <<",\"margin\":1,\"extra_layer\":0,\"border_width\":1,\"envelope_width\":1,\"block_size\":16,\"inflate_lu\":0.001,\n"
            <<"\"lattice_shape\":["<<nx<<','<<ny<<','<<nz<<"],\"lattice_bbox\":["<<b.x0<<','<<b.x1<<','<<b.y0<<','<<b.y1<<','<<b.z0<<','<<b.z1<<"],\n"
            <<"\"preinflate_bbox_lu\":[["<<xr[0]<<','<<yr[0]<<','<<zr[0]<<"],["<<xr[1]<<','<<yr[1]<<','<<zr[1]<<"]],\n"
            <<"\"inflated_bbox_lu\":[["<<ix[0]<<','<<iy[0]<<','<<iz[0]<<"],["<<ix[1]<<','<<iy[1]<<','<<iz[1]<<"]],\n"
            <<"\"closed_fluid_voxels\":"<<fluid<<",\"closed_solid_voxels\":"<<n-fluid<<",\"sparse_blocks\":"<<bulks.size()<<",\n\"allocated_native_counts\":[";
        for(int i=0;i<6;++i) meta<<(i?",":"")<<counts[i];
        meta<<"],\"unallocated_native_code\":255,\"fluid_timesteps_run\":0,\"mpi_ranks\":1}\n";
        if(!meta || counts[0] || counts[5]) throw std::runtime_error("Incomplete voxelization/export");
        pcout<<"CLOSED_VOXELIZATION_COMPLETE shape="<<nx<<','<<ny<<','<<nz<<" fluid="<<fluid
             <<" inner_border="<<counts[4]<<" FLUID_TIMESTEPS_RUN=0"<<std::endl;
        return 0;
    } catch(std::exception const& e) { pcout<<"STEP2_ERROR "<<e.what()<<std::endl; return 2; }
}
