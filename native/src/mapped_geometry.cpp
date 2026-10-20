#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
static std::vector<int64_t> crossing_faces,crossing_bonds;
static bool triangle(const double* q,const double* n,const double* a,const double* b,const double* c,double lower,double upper){
 double e1[3],e2[3],h[3],s[3],cross[3];
 for(int k=0;k<3;k++){e1[k]=b[k]-a[k];e2[k]=c[k]-a[k];s[k]=q[k]-a[k];}
 for(int k=0;k<3;k++)h[k]=n[(k+1)%3]*e2[(k+2)%3]-n[(k+2)%3]*e2[(k+1)%3];
 double det=0;for(int k=0;k<3;k++)det+=e1[k]*h[k];if(std::abs(det)<1e-28)return false;
 double u=0;for(int k=0;k<3;k++)u+=s[k]*h[k];u/=det;if(u<-1e-9||u>1+1e-9)return false;
 for(int k=0;k<3;k++)cross[k]=s[(k+1)%3]*e1[(k+2)%3]-s[(k+2)%3]*e1[(k+1)%3];
 double v=0,t=0;for(int k=0;k<3;k++){v+=n[k]*cross[k];t+=e2[k]*cross[k];}v/=det;t/=det;
 return v>=-1e-9&&u+v<=1+1e-9&&t>=lower&&t<=upper;
}
extern "C" int64_t crossings_build(int nb,const double* X,const int64_t* bonds,const int64_t* ptr,const int64_t* candidates,const double* quads){
 crossing_faces.clear();crossing_bonds.clear();
 for(int b=0;b<nb;b++){
  const double* p=X+3*bonds[2*b];const double* q=X+3*bonds[2*b+1];double n[3];for(int k=0;k<3;k++)n[k]=q[k]-p[k];
  for(int64_t c=ptr[b];c<ptr[b+1];c++)for(int f=0;f<6;f++){
   int64_t face=6*candidates[c]+f;const double* v=quads+12*face;
   if(triangle(p,n,v,v+3,v+6,-1e-9,1+1e-9)||triangle(p,n,v,v+6,v+9,-1e-9,1+1e-9)){crossing_faces.push_back(face);crossing_bonds.push_back(b);}
  }
 }
 return crossing_faces.size();
}
extern "C" void crossings_copy(int64_t* faces,int64_t* bonds){
 std::copy(crossing_faces.begin(),crossing_faces.end(),faces);std::copy(crossing_bonds.begin(),crossing_bonds.end(),bonds);
 std::vector<int64_t>().swap(crossing_faces);std::vector<int64_t>().swap(crossing_bonds);
}
// The ray is clipped against each of six mapped reference tetrahedra. Invalid
// local fits remain circumsphere obstacles, matching the conservative V3 rule.
extern "C" void ray_mapped_cells(int nq,int nc,const double* q,const double* normal,const int64_t* owner,const double* x,const double* origin,const double* inverse,const double* radius,const unsigned char* valid,double probe,unsigned char* closed,unsigned char* blocked){
 for(int a=0;a<nq;a++){
  closed[a]=blocked[a]=0;
  for(int j=0;j<nc;j++){
   if(j==owner[a])continue;double dist=0,dot=0;
   for(int k=0;k<3;k++){double p=q[3*a+k]-x[3*j+k];dist+=p*p;dot+=p*normal[3*a+k];}
   double rad=radius[j],perp=dist-dot*dot;if(perp>rad*rad||dot>rad)continue;
   if(!valid[j]){double r=std::sqrt(std::max(0.,rad*rad-perp)),lo=-dot-r,hi=-dot+r;if(hi>=std::max(lo,probe)){blocked[a]=1;if(lo<=probe){closed[a]=1;break;}}continue;}
   for(int t=0;t<6;t++){
    int idx=j*6+t;double p[4]={0,0,0,1},v[4]={0,0,0,0};
    for(int k=0;k<3;k++){for(int l=0;l<3;l++){p[k]+=inverse[9*idx+3*k+l]*(q[3*a+l]-origin[3*idx+l]);v[k]+=inverse[9*idx+3*k+l]*normal[3*a+l];}p[3]-=p[k];v[3]-=v[k];}
    double lo=-1e300,hi=1e300;
    for(int k=0;k<4;k++){if(std::abs(v[k])<1e-8){if(p[k]<-1e-9){hi=-1e300;break;}}else {double hit=(-1e-9-p[k])/v[k];if(v[k]>0)lo=std::max(lo,hit);else hi=std::min(hi,hit);}}
    if(hi>=std::max(lo,probe)){blocked[a]=1;if(lo<=probe){closed[a]=1;break;}}
   }
   if(closed[a])break;
  }
 }
}
extern "C" void ray_triangles(int nq,int nt,const double* q,const double* normal,const double* tri,double probe,unsigned char* blocked){
 for(int a=0;a<nq;a++)for(int t=0;t<nt&&!blocked[a];t++)if(triangle(q+3*a,normal+3*a,tri+9*t,tri+9*t+3,tri+9*t+6,probe,1e300))blocked[a]=1;
}
