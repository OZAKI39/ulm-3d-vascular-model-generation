#include "workflow_geometry.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("USAGE root output_csv");std::string root=argv[1];workflow::Geometry geom(root);frozen::FlowField field(root+"/fields/FROZEN_FLOW_FIELD_V0.h5");
 std::ofstream out(argv[2]);out<<std::setprecision(17)<<"mode,sample,x_m,y_m,z_m,radius_m,positive_inward_velocity_m_s,wall_distance_m,disclaimer\n";
 for(bool weighted:{false,true}){
  workflow::PositionSampler position(root,field,geom,field.dx,2e-9,4096,weighted);std::mt19937_64 rng(123);workflow::Candidate c{1,9.683592065545495e-7,0,.5};
  for(int i=0;i<3000;i++){rigid::Vec x{};if(!position.place(c,x,rng,{},{}))throw std::runtime_error("POSITION_TEST_EXHAUSTED");auto q=frozen::FlowFieldSampler(field).query(x);double un=-rigid::dot(q.u,geom.ports[0].normal);if(q.status!=frozen::Status::VALID||un<=0||geom.wall.distance(x)<c.radius+field.dx)throw std::runtime_error("INVALID_POSITION_ADMITTED");
   out<<(weighted?"VALID_SUPPORT_FLUX_WEIGHTED":"VALID_AREA_UNIFORM")<<','<<i;for(double v:x)out<<','<<v;out<<','<<c.radius<<','<<un<<','<<geom.wall.distance(x)<<",NOT EXPERIMENTAL CONCENTRATION\n";
  }
 }
 std::cout<<"{\"status\":\"PASS\",\"samples_per_mode\":3000,\"scope\":\"NATIVE_POSITION_DISTRIBUTION_DIAGNOSTIC_NOT_DYNAMIC_ADMISSIONS\",\"disclaimer\":\"NOT EXPERIMENTAL CONCENTRATION\"}\n";
 return 0;
 }catch(std::exception const&e){std::cerr<<e.what()<<std::endl;return 2;}}
