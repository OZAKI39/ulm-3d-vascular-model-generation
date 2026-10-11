// Read-only energy/force evaluator for the locked equal-volume, uniform-K/G LPS.
// Never linked into LAMMPS; no trajectory changes. m is frozen initial weighted volume.
#include <cmath>
#include <vector>
static double core(int n,int nb,const double*x,const int*b,const unsigned char*alive,
 const double*r0,const double*h,const double*m,double V,double K,double G,double*force,double L) {
 std::vector<double> theta(n,0.0), e(nb,0.0);
 for(int k=0;k<nb;k++)if(alive[k]){int i=b[2*k],j=b[2*k+1];double r=0;
  for(int d=0;d<3;d++){double z=x[3*j+d]-x[3*i+d];if(d<2 && L>0)z-=L*std::round(z/L);r+=z*z;}e[k]=std::sqrt(r)-r0[k];
  if(std::abs(e[k])<2.2204e-16)e[k]=0;
  theta[i]+=e[k]*V*h[k];theta[j]+=e[k]*V*h[k];}
 double U=0;
 for(int i=0;i<n;i++){theta[i]*=3/m[i];U+=0.5*(K-5*G/3)*V*theta[i]*theta[i];}
 if(force)for(int q=0;q<3*n;q++)force[q]=0;
 for(int k=0;k<nb;k++)if(alive[k]){int i=b[2*k],j=b[2*k+1];
  U+=7.5*G*V*V*h[k]/r0[k]*e[k]*e[k]*(1/m[i]+1/m[j]);
  if(force){double r=r0[k]+e[k];double f=V*V*h[k]*((3*K-5*G)*(theta[i]/m[i]+theta[j]/m[j])+15*G/r0[k]*(1/m[i]+1/m[j])*e[k]);
   for(int d=0;d<3;d++){double delta=x[3*j+d]-x[3*i+d];if(d<2 && L>0)delta-=L*std::round(delta/L);double z=f*delta/r;force[3*i+d]+=z;force[3*j+d]-=z;}}}
 return U;
}

extern "C" double evaluate_periodic(int n,int nb,const double*x,const int*b,const unsigned char*alive,const double*r0,const double*h,const double*m,double V,double K,double G,double*force,double L){return core(n,nb,x,b,alive,r0,h,m,V,K,G,force,L);}
