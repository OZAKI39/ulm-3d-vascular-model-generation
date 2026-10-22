#include <cstdint>
#include <limits>
// Original face + clamped-edge arithmetic, same strict argmin tie order.
// -ffp-contract=off: no fused multiply-add or fast-math reassociation.
static double dot(const double*a,const double*b){return (a[0]*b[0]+a[1]*b[1])+a[2]*b[2];}
// NumPy's active stride-3 einsum loop reduces (0+2)+1, whereas np.sum
// uses (0+1)+2. Keep both orders explicit; verified on the production host.
static double eindot(const double*a,const double*b){return (a[0]*b[0]+a[2]*b[2])+a[1]*b[1];}
extern "C" void closest(const double*p,const double*t,int64_t count,double*out,double*bary){
  for(int64_t k=0;k<count;++k){
    const double*a=t+9*k;double e[3],f[3],v[3];
    for(int j=0;j<3;++j){e[j]=a[3+j]-a[j];f[j]=a[6+j]-a[j];v[j]=p[j]-a[j];}
    double ee=eindot(e,e),ef=eindot(e,f),ff=eindot(f,f),ev=eindot(e,v),fv=eindot(f,v),den=ee*ff-ef*ef;
    double u=den>0?(ff*ev-ef*fv)/den:0.,w=den>0?(ee*fv-ef*ev)/den:0.;
    double best[3],bw[3]={1-u-w,u,w},d[3];
    for(int j=0;j<3;++j){best[j]=(a[j]+u*e[j])+w*f[j];d[j]=p[j]-best[j];}
    double dist=(den>0&&bw[0]>=0&&u>=0&&w>=0)?dot(d,d):std::numeric_limits<double>::infinity();
    const int edges[3][2]={{0,1},{1,2},{2,0}};
    for(auto &edge:edges){
      int i=edge[0],j=edge[1];double e2[3],v2[3],cp[3];
      for(int h=0;h<3;++h){e2[h]=a[3*j+h]-a[3*i+h];v2[h]=p[h]-a[3*i+h];}
      double len=dot(e2,e2),s=len>0?dot(v2,e2)/len:0.;s=s<0?0:s>1?1:s;
      for(int h=0;h<3;++h){cp[h]=a[3*i+h]+s*e2[h];d[h]=p[h]-cp[h];}
      double ds=dot(d,d);
      if(ds<dist){dist=ds;for(int h=0;h<3;++h){best[h]=cp[h];bw[h]=0.;}bw[i]=1-s;bw[j]=s;}
    }
    for(int j=0;j<3;++j){out[3*k+j]=best[j];bary[3*k+j]=bw[j];}
  }
}
