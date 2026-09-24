// P8.2 diagnostic only: double-precision P1 advection. No MB/wall physics.
#include <cmath>
#include <algorithm>
#include <cstdint>
extern "C" {
struct Mesh { const double *o,*inv,*tol,*vel; const int64_t *tet,*nb; const int *role; int64_t n; };
static Mesh m;
void p82_init(const double* o,const double* inv,const double* tol,const double* v,
              const int64_t* tet,const int64_t* nb,const int* role,int64_t n) {m={o,inv,tol,v,tet,nb,role,n};}
static void weights(const double*x,int64_t c,double*w) {
    w[0]=1.; for(int i=0;i<3;i++){w[i+1]=0.;for(int j=0;j<3;j++)w[i+1]+=m.inv[9*c+3*i+j]*(x[j]-m.o[3*c+j]);w[0]-=w[i+1];}
}
static bool field(const double*x,int64_t &c,double*v,int64_t &outside,int &face) {
    for(int iteration=0;iteration<256;iteration++) {
        if(c<0||c>=m.n)return false;
        double w[4];weights(x,c,w);int f=0;for(int j=1;j<4;j++)if(w[j]<w[f])f=j;
        if(w[f]>=-m.tol[c]){for(int j=0;j<3;j++){v[j]=0.;for(int k=0;k<4;k++)v[j]+=w[k]*m.vel[3*m.tet[4*c+k]+j];}return true;}
        int64_t next=m.nb[4*c+f];if(next<0){outside=c;face=f;return false;}c=next;
    } return false;
}
int p82_sample(const double*x,int64_t cell,double*v){int64_t outside=-1;int face=-1;return field(x,cell,v,outside,face)?int(cell):-1;}
// status: 1 outlet, 2 time horizon, 3 step budget, 4 dt floor, 5 locator.
int p82_trace(const double*initial,int64_t cell,double abs_tol,double spatial_step,
              double max_time,double max_path,int max_steps,double*out,int*count,int*outlet) {
    double x[3]={initial[0],initial[1],initial[2]},t=0.,dt=1e-4,path_length=0.;*count=1;*outlet=0;
    out[0]=0.;for(int j=0;j<3;j++)out[j+1]=x[j];
    for(int attempt=0;attempt<max_steps*20;attempt++) {
        if(*count>=max_steps)return 3;if(t>=max_time)return 2;if(path_length>=max_path)return 6;
        double k1[3],k2[3],k3[3],k4[3],z[3],y[3];int64_t outside=-1;int face=-1;
        if(!field(x,cell,k1,outside,face))return 5;
        double speed=std::sqrt(k1[0]*k1[0]+k1[1]*k1[1]+k1[2]*k1[2]);
        dt=std::min(dt,std::min(.01,std::min(max_time-t,spatial_step/std::max(speed,1e-30))));
        if(dt<1e-13||t+dt==t)return 4;
        int64_t c2=cell,c3=cell,c4=cell;
        for(int j=0;j<3;j++)z[j]=x[j]+.5*dt*k1[j];
        bool ok=field(z,c2,k2,outside,face);
        if(ok){for(int j=0;j<3;j++)z[j]=x[j]+.75*dt*k2[j];ok=field(z,c3,k3,outside,face);}
        if(ok){for(int j=0;j<3;j++)y[j]=x[j]+dt*(2./9*k1[j]+1./3*k2[j]+4./9*k3[j]);ok=field(y,c4,k4,outside,face);}
        if(!ok) {
            // Resolve only an actual open cap crossed by a local Euler ODE step.
            // Wall encounters reject and halve dt; never project onto the wall.
            if(outside>=0 && face>=0 && m.role[4*outside+face]>0 && m.role[4*outside+face]<=3) {
                double w[4],we[4],end[3];weights(x,outside,w);
                for(int j=0;j<3;j++)end[j]=x[j]+dt*k1[j];weights(end,outside,we);
                double slope=we[face]-w[face];double alpha=slope<0?-w[face]/slope:-1;
                if(alpha>=0&&alpha<=1 && dt*alpha*speed<=spatial_step*.02) {
                    double wh[4];for(int j=0;j<3;j++)end[j]=x[j]+alpha*dt*k1[j];weights(end,outside,wh);
                    bool on=true;for(int k=0;k<4;k++)if(wh[k]<-m.tol[outside]*4)on=false;
                    if(on){t+=alpha*dt;int n=(*count)++;out[4*n]=t;for(int j=0;j<3;j++)out[4*n+j+1]=end[j];*outlet=m.role[4*outside+face];return 1;}
                }
            }
            dt*=.5;continue;
        }
        double error=0.;for(int j=0;j<3;j++){double low=x[j]+dt*(7./24*k1[j]+.25*k2[j]+1./3*k3[j]+.125*k4[j]);error=std::max(error,std::abs(y[j]-low));}
        if(error>abs_tol){dt*=std::max(.1,.8*std::cbrt(abs_tol/error));continue;}
        double ds=0.;for(int j=0;j<3;j++)ds+=(y[j]-x[j])*(y[j]-x[j]);path_length+=std::sqrt(ds);
        t+=dt;cell=c4;int n=(*count)++;out[4*n]=t;for(int j=0;j<3;j++){x[j]=y[j];out[4*n+j+1]=x[j];}
        dt*=error>0?std::min(2.,.9*std::cbrt(abs_tol/error)):2.;
    }return 3;
}
}
