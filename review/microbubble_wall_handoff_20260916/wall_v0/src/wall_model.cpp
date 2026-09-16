#include "wall_model.hpp"
#include <hdf5.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace wallv0 {
using namespace rigid;
Frame frame(Vec n,const Frame*previous){double d=norm(n);if(!(d>0))throw std::runtime_error("INVALID_NORMAL");n=mul(n,1/d);int axis=0;for(int j=1;j<3;j++)if(std::abs(n[j])<std::abs(n[axis]))axis=j;Vec t{};t[axis]=1;t=sub(t,mul(n,dot(n,t)));t=mul(t,1/norm(t));
 if(previous){Vec old={(*previous)[0],(*previous)[3],(*previous)[6]};if(dot(t,old)<0)t=mul(t,-1);}Vec t2=cross(n,t);Frame Q;for(int i=0;i<3;i++){Q[3*i]=t[i];Q[3*i+1]=t2[i];Q[3*i+2]=n[i];}return Q;
}
Matrix6 rotate(const Matrix6&A,const Frame&Q){Matrix6 out{};for(int i=0;i<6;i++)for(int j=0;j<6;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++)out[6*i+j]+=Q[3*(i%3)+k]*A[6*(i/3*3+k)+j/3*3+l]*Q[3*(j%3)+l];return out;}
WallResistanceLookup::WallResistanceLookup(const std::string&path){
 hid_t file=H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);if(file<0)throw std::runtime_error("TABLE_OPEN_FAILED");
 hid_t ed=H5Dopen2(file,"epsilon",H5P_DEFAULT);hid_t es=H5Dget_space(ed);hsize_t n;H5Sget_simple_extent_dims(es,&n,nullptr);epsilon.resize(n);logs.resize(n);
 if(H5Dread(ed,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,epsilon.data())<0)throw std::runtime_error("TABLE_GRID_READ_FAILED");H5Sclose(es);H5Dclose(ed);
 hid_t rd=H5Dopen2(file,"R_wall_excess_scaled",H5P_DEFAULT);hid_t rs=H5Dget_space(rd);hsize_t shape[3];int nd=H5Sget_simple_extent_dims(rs,shape,nullptr);if(nd!=3||shape[0]!=n||shape[1]!=6||shape[2]!=6)throw std::runtime_error("TABLE_MATRIX_SHAPE_FAILED");excess.resize(n);
 if(H5Dread(rd,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,excess.data())<0)throw std::runtime_error("TABLE_MATRIX_READ_FAILED");H5Sclose(rs);H5Dclose(rd);H5Fclose(file);
 if(n<800||epsilon.front()!=.001||epsilon.back()!=20.)throw std::runtime_error("TABLE_RANGE_FAILED");
 for(size_t i=0;i<n;i++){if(i&&epsilon[i]<=epsilon[i-1])throw std::runtime_error("TABLE_ORDER_FAILED");logs[i]=std::log(epsilon[i]);for(double v:excess[i])if(!std::isfinite(v))throw std::runtime_error("NONFINITE_TABLE");}
}
Matrix6 WallResistanceLookup::scaled_excess(double e)const{if(!std::isfinite(e))throw std::runtime_error("NONFINITE_GAP");if(e>20)return{};if(e<=epsilon.front())return excess.front();auto it=std::lower_bound(epsilon.begin(),epsilon.end(),e);size_t j=it-epsilon.begin();if(*it==e)return excess[j];size_t i=j-1;double w=(std::log(e)-logs[i])/(logs[j]-logs[i]);Matrix6 r;for(int k=0;k<36;k++)r[k]=(1-w)*excess[i][k]+w*excess[j][k];return r;}
Matrix6 WallResistanceLookup::global_excess(double a,double mu,double h,const Frame&Q,bool enabled)const{if(!enabled)return{};auto r=scaled_excess(h/a);double b=6*std::acos(-1.)*mu*a,c=8*std::acos(-1.)*mu*a*a*a;double d[6]={std::sqrt(b),std::sqrt(b),std::sqrt(b),std::sqrt(c),std::sqrt(c),std::sqrt(c)};for(int i=0;i<6;i++)for(int j=0;j<6;j++)r[6*i+j]*=d[i]*d[j];return rotate(r,Q);}
passive::WallQuery STLNearestWallQuery::query(Vec x,double a)const{if(!flat)return mesh->query(x,a);passive::WallQuery q;q.closest={x[0],x[1],0};q.distance=x[2];q.gap=x[2]-a;q.normal=q.triangle_normal=q.triangle_normal_fluid=q.plane_normal={0,0,1};q.triangle_id=0;q.patch_triangles=1;q.valid=x[2]>0;q.planar=q.valid;return q;}
bool STLNearestWallQuery::swept_safe(Vec x,Vec y,double a)const{return flat?std::min(x[2],y[2])-a>=-1e-12:mesh->segment_safe(x,y,a-1e-12);}
std::vector<WallBlock> blocks(const std::vector<Entry>&stage,const std::vector<Entry>&base,bool hard){std::vector<WallBlock>b;for(size_t i=0;i<stage.size();i++)b.push_back({int(i),base[i].geometry.triangle_id,stage[i].excess,base[i].geometry.normal,base[i].geometry.gap,hard});return b;}
}
extern "C" void*wall_lookup_open(const char*p){try{return new wallv0::WallResistanceLookup(p);}catch(...){return nullptr;}}
extern "C" void wall_lookup_close(void*h){delete static_cast<wallv0::WallResistanceLookup*>(h);}
extern "C" int wall_lookup_evaluate(void*h,int n,const double*e,double*out){try{for(int k=0;k<n;k++){auto r=static_cast<wallv0::WallResistanceLookup*>(h)->scaled_excess(e[k]);std::copy(r.begin(),r.end(),out+36*k);}return 0;}catch(...){return 1;}}
extern "C" int wall_rotation_audit(int n,const double*normals,const double*matrices,const double*velocities,double*frames,double*f1,double*f2){try{for(int k=0;k<n;k++){auto Q=wallv0::frame({normals[3*k],normals[3*k+1],normals[3*k+2]});std::copy(Q.begin(),Q.end(),frames+9*k);wallv0::Matrix6 M;std::copy(matrices+36*k,matrices+36*(k+1),M.begin());auto G=wallv0::rotate(M,Q);double v[6]={},z[6]={};
 for(int i=0;i<6;i++)for(int j=0;j<3;j++)v[i]+=Q[3*j+i%3]*velocities[6*k+3*(i/3)+j];
 for(int i=0;i<6;i++)for(int j=0;j<6;j++){f1[6*k+i]+=G[6*i+j]*velocities[6*k+j];z[i]+=M[6*i+j]*v[j];}
 for(int i=0;i<6;i++)for(int j=0;j<3;j++)f2[6*k+i]+=Q[3*(i%3)+j]*z[3*(i/3)+j];}return 0;}catch(...){return 1;}}
extern "C" void*wall_geometry_open(const char*p){try{return new wallv0::STLNearestWallQuery(p);}catch(...){return nullptr;}}
extern "C" void wall_geometry_close(void*h){delete static_cast<wallv0::STLNearestWallQuery*>(h);}
extern "C" int wall_geometry_query(void*h,int n,const double*x,const double*a,double*out){try{for(int i=0;i<n;i++){auto q=static_cast<wallv0::STLNearestWallQuery*>(h)->query({x[3*i],x[3*i+1],x[3*i+2]},a[i]);double*p=out+27*i;for(int k=0;k<3;k++){p[k]=q.closest[k];p[3+k]=q.normal[k];p[6+k]=q.triangle_normal[k];p[9+k]=q.plane_normal[k];p[24+k]=q.triangle_normal_fluid[k];}p[12]=q.distance;p[13]=q.gap;p[14]=q.triangle_id;p[15]=q.rms_over_a;p[16]=q.normal_spread_deg;p[17]=q.nearest_normal_angle_deg;p[18]=q.patch_triangles;p[19]=q.valid;p[20]=q.planar;p[21]=q.ambiguous;p[22]=0;p[23]=0;}return 0;}catch(...){return 1;}}
extern "C" int wall_geometry_swept(void*h,int n,const double*x,const double*y,const double*a,int*out){try{for(int i=0;i<n;i++)out[i]=static_cast<wallv0::STLNearestWallQuery*>(h)->swept_safe({x[3*i],x[3*i+1],x[3*i+2]},{y[3*i],y[3*i+1],y[3*i+2]},a[i]);return 0;}catch(...){return 1;}}
extern "C" int wall_single_response(void*h,int n,const double*state,double*out){try{for(int k=0;k<n;k++){auto s=state+11*k;double a=s[0],gap=s[1];std::vector<rigid::Particle>p={{1,0,{0,0,a+gap},a}};std::vector<rigid::Background>b(1);for(int j=0;j<3;j++){b[0].du[j]=s[2+j]/(6*std::acos(-1.)*.001*a);b[0].dw[j]=s[5+j]/(8*std::acos(-1.)*.001*a*a*a);}auto Q=wallv0::frame({s[8],s[9],s[10]});auto E=static_cast<wallv0::WallResistanceLookup*>(h)->global_excess(a,.001,gap,Q);std::vector<rigid::WallBlock>w={{0,0,E,{s[8],s[9],s[10]},gap,false}};auto q=rigid::solve(p,b,{},nullptr,0,&w);std::copy(q.q.begin(),q.q.end(),out+6*k);}return 0;}catch(...){return 1;}}
