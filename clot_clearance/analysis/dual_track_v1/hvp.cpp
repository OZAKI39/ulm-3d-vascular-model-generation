// Exact directional derivative of isolated equal-volume LPS force, nonperiodic.
// Recovery diagnostic only. Production energy and LAMMPS remain unchanged.
#include <cmath>
#include <vector>
extern "C" void hessian_vector(int n,int nb,const double*x,const int*b,const unsigned char*a,const double*r0,const double*h,const double*m,double V,double K,double G,const double*v,double*out){
 std::vector<double> th(n,0),dt(n,0),r(nb),e(nb),de(nb);
 for(int k=0;k<nb;k++)if(a[k]){int i=b[2*k],j=b[2*k+1];double r2=0,dot=0;for(int d=0;d<3;d++){double q=x[3*j+d]-x[3*i+d];r2+=q*q;dot+=q*(v[3*j+d]-v[3*i+d]);}r[k]=std::sqrt(r2);e[k]=r[k]-r0[k];if(std::abs(e[k])<2.2204e-16)e[k]=0;de[k]=dot/r[k];th[i]+=h[k]*e[k];th[j]+=h[k]*e[k];dt[i]+=h[k]*de[k];dt[j]+=h[k]*de[k];}
 for(int i=0;i<n;i++){th[i]*=3*V/m[i];dt[i]*=3*V/m[i];for(int d=0;d<3;d++)out[3*i+d]=0;}
 for(int k=0;k<nb;k++)if(a[k]){int i=b[2*k],j=b[2*k+1];double c=V*V*h[k]*((3*K-5*G)*(th[i]/m[i]+th[j]/m[j])+15*G*e[k]/r0[k]*(1/m[i]+1/m[j]));double dc=V*V*h[k]*((3*K-5*G)*(dt[i]/m[i]+dt[j]/m[j])+15*G*de[k]/r0[k]*(1/m[i]+1/m[j]));for(int d=0;d<3;d++){double q=(x[3*j+d]-x[3*i+d])/r[k];double w=dc*q+c/r[k]*(v[3*j+d]-v[3*i+d]-q*de[k]);out[3*i+d]-=w;out[3*j+d]+=w;}}
}
