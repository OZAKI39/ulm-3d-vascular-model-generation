#include <cmath>
#include <algorithm>
// Conservative exterior straight-ray certificate for approximate material cells.
// Valid cells are affine images of the original dx^3 dual cubes. Rank-deficient
// cells are unresolved obstacles represented by their reference circumsphere.
extern "C" void ray_cells(int nq,int nc,const double* q,const double* normal,
 const long long* owner,const double* x,const double* inv,const double* radius,
 const unsigned char* valid,double half,double probe,unsigned char* closed,unsigned char* blocked){
 for(int a=0;a<nq;a++){
  closed[a]=0;blocked[a]=0;
  for(int j=0;j<nc;j++){
   if(j==owner[a])continue;
   double p[3],n[3],dist=0,dot=0;
   for(int k=0;k<3;k++){p[k]=q[3*a+k]-x[3*j+k];n[k]=normal[3*a+k];dist+=p[k]*p[k];dot+=p[k]*n[k];}
   double rad=radius[j],perp=dist-dot*dot;
   if(perp>rad*rad||dot>rad)continue;
   double lo=-1e300,hi=1e300;
   if(valid[j]){
    for(int k=0;k<3;k++){
     double pp=0,nn=0;
     for(int l=0;l<3;l++){pp+=inv[9*j+3*k+l]*p[l];nn+=inv[9*j+3*k+l]*n[l];}
     if(std::abs(nn)<1e-12){if(std::abs(pp)>half){hi=-1e300;break;}}
     else {double t0=(-half-pp)/nn,t1=(half-pp)/nn;lo=std::max(lo,std::min(t0,t1));hi=std::min(hi,std::max(t0,t1));}
    }
   }else {double r=std::sqrt(std::max(0.,rad*rad-perp));lo=-dot-r;hi=-dot+r;}
   if(hi>=std::max(lo,probe)){
    blocked[a]=1;
    if(lo<=probe&&hi>=probe){closed[a]=1;break;}
   }
  }
 }
}
extern "C" void ray_triangles(int nq,int nt,const double* q,const double* normal,const double* tri,double probe,unsigned char* blocked){
 for(int i=0;i<nq;i++){
  for(int j=0;j<nt&&!blocked[i];j++){
   const double *p=tri+9*j,*n=normal+3*i;
   double e1[3],e2[3],h[3],s[3],cross[3];
   for(int k=0;k<3;k++){e1[k]=p[3+k]-p[k];e2[k]=p[6+k]-p[k];s[k]=q[3*i+k]-p[k];}
   for(int k=0;k<3;k++)h[k]=n[(k+1)%3]*e2[(k+2)%3]-n[(k+2)%3]*e2[(k+1)%3];
   double a=0;for(int k=0;k<3;k++)a+=e1[k]*h[k];if(std::abs(a)<1e-24)continue;
   double u=0;for(int k=0;k<3;k++)u+=s[k]*h[k];u/=a;if(u<0||u>1)continue;
   for(int k=0;k<3;k++)cross[k]=s[(k+1)%3]*e1[(k+2)%3]-s[(k+2)%3]*e1[(k+1)%3];
   double v=0,t=0;for(int k=0;k<3;k++){v+=n[k]*cross[k];t+=e2[k]*cross[k];}v/=a;t/=a;
   if(v>=0&&u+v<=1&&t>=probe)blocked[i]=1;
  }
 }
}
