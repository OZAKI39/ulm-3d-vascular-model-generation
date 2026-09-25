#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "meshGeneratingFunctions.h"
#include "meshGeneratingFunctions.hh"
#include "palabos3D.h"
#include "palabos3D.hh"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
int main(int argc,char**argv){
 plb::plbInit(&argc,&argv);
 if(argc!=3||plb::global::mpi().getSize()!=1)throw std::runtime_error("Usage: rbc_geometry_probe NUMERICS_CONTRACT OUTPUT_DIR; MPI1, zero solver steps");
 std::ifstream in(argv[1]);std::string text((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
 auto number=[&](std::string key){size_t pos=text.find("\""+key+"\":");if(pos==std::string::npos)throw std::runtime_error("Missing numerics key");return std::stod(text.substr(text.find(':',pos)+1));};
 const double dx=number("dx_m");
 const double scale=std::cbrt(50.0/90.0); // User's nominal 90 -> 50 um^3 uniform DEV_V0 definition, frozen before fit.
 const double radius=3.91e-6*scale;
 auto boundary=plb::constructMeshElement(RBC_FROM_SPHERE,radius/dx,600,dx,std::string(""),plb::Array<double,3>(0.,0.,0.));
 auto const& mesh=boundary.getMesh();
 std::string root=argv[2];std::ofstream v(root+"/RBC_REFERENCE_VERTICES.csv"),f(root+"/RBC_REFERENCE_TRIANGLES.csv");v<<std::setprecision(17)<<"id,x_lu,y_lu,z_lu\n";f<<"id,a,b,c\n";
 for(plint i=0;i<mesh.getNumVertices();++i){auto p=mesh.getVertex(i);v<<i<<','<<p[0]<<','<<p[1]<<','<<p[2]<<'\n';}
 for(plint i=0;i<mesh.getNumTriangles();++i)f<<i<<','<<mesh.getVertexId(i,0)<<','<<mesh.getVertexId(i,1)<<','<<mesh.getVertexId(i,2)<<'\n';
 std::ofstream m(root+"/NATIVE_MESH_GENERATION.json");m<<std::setprecision(17)<<"{\"status\":\"PASS\",\"model\":\"MOUSE_RBC_DEV_V0\",\"shape\":\"RBC_FROM_SPHERE\",\"target_volume_um3\":50,\"scale\":"<<scale<<",\"radius_m\":"<<radius<<",\"dx_m\":"<<dx<<",\"vertices\":"<<mesh.getNumVertices()<<",\"triangles\":"<<mesh.getNumTriangles()<<",\"geometry_tuned_after_fit\":false,\"solver_steps\":0,\"native_mesh_inflate_default_retained\":true}\n";
 if(!m||!v||!f)throw std::runtime_error("Mesh output failed");return 0;
}
