#include <cmath>
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>
using std::log;using std::pow;
constexpr double MY_PI=3.141592653589793238462643383279502884;
std::array<double,2> target_poly(double radi,double radj,double gap,double mu,int flaglog){
 double a_sq=0,a_sh=0,a_pu=0,h_sep=gap/radi;
 double beta0=radj/radi,beta1=1+beta0;
 double beta[2][5]={{0},{0}},pre[2];beta[0][1]=beta0;beta[1][1]=beta1;pre[1]=8.0*(pre[0]=MY_PI*mu*radi)*radi*radi;pre[0]*=6.0;
        if (flaglog) {
          a_sq = beta0*beta0/beta1/beta1/h_sep +
            (1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3.0)*log(1.0/h_sep);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0 *
                   pow(beta0,3.0)+pow(beta0,4.0))/21.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
          a_sq *= 6.0*MY_PI*mu*radi;
          a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/pow(beta1,3.0) *
            log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3.0) +
                       16.0*pow(beta0,4.0))/375.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
          a_sh *= 6.0*MY_PI*mu*radi;
          // old invalid eq for pumping term
          // changed 29Jul16 from eq 9.25 -> 9.27 in Kim and Karilla
//          a_pu = beta0*(4.0+beta0)/10.0/beta1/beta1*log(1.0/h_sep);
//          a_pu += (32.0-33.0*beta0+83.0*beta0*beta0+43.0 *
//                   pow(beta0,3.0))/250.0/pow(beta1,3.0)*h_sep*log(1.0/h_sep);
//          a_pu *= 8.0*MY_PI*mu*pow(radi,3.0);
          a_pu = 2.0*beta0/5.0/beta1*log(1.0/h_sep);
          a_pu += 2.0*(8.0+6.0*beta0+33.0*beta0*beta0)/125.0/beta1/beta1*
                   h_sep*log(1.0/h_sep);
          a_pu *= 8.0*MY_PI*mu*pow(radi,3.0);
        } else a_sq = 6.0*MY_PI*mu*radi*(beta0*beta0/beta1/beta1/h_sep);
 return {a_sq,a_sh};}
std::array<double,2> target_upoly(double radi,double radj,double gap,double mu,int flaglog){
 double a_sq=0,a_sh=0,a_pu=0,h_sep=gap/radi;
 double beta0=radj/radi,beta1=1+beta0;
 double beta[2][5]={{0},{0}},pre[2];beta[0][1]=beta0;beta[1][1]=beta1;pre[1]=8.0*(pre[0]=MY_PI*mu*radi)*radi*radi;pre[0]*=6.0;
        if (flaglog) {
          beta[0][2] = beta[0][1]*beta[0][1];
          beta[0][3] = beta[0][2]*beta[0][1];
          beta[0][4] = beta[0][3]*beta[0][1];
          beta[1][2] = beta[1][1]*beta[1][1];
          beta[1][3] = beta[1][2]*beta[1][1];
          double log_h_sep_beta13 = log(1.0/h_sep)/beta[1][3];
          double h_sep_beta11 = h_sep/beta[1][1];

          a_sq = pre[0]*(beta[0][2]/beta[1][2]/h_sep
                +((0.2+1.4*beta[0][1]+0.2*beta[0][2])
                  +(1.0+18.0*(beta[0][1]+beta[0][3])-29.0*beta[0][2]
                    +beta[0][4])*h_sep_beta11/21.0)*log_h_sep_beta13);

          a_sh = pre[0]*((8.0*(beta[0][1]+beta[0][3])+4.0*beta[0][2])/15.0
                +(64.0-180.0*(beta[0][1]+beta[0][3])+232.0*beta[0][2]
                  +64.0*beta[0][4])*h_sep_beta11/375.0)*log_h_sep_beta13;

          /*a_sq = beta0*beta0/beta1/beta1/h_sep
                  +(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3)*log(1.0/h_sep);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0*pow(beta0,3)
                  +pow(beta0,4))/21.0/pow(beta1,4)*h_sep*log(1.0/h_sep);
          a_sq *= 6.0*MY_PI*mu*radi;

          a_sh = 4.0*beta0*(2.0+beta0
                  +2.0*beta0*beta0)/15.0/pow(beta1,3)*log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3)
                  +16.0*pow(beta0,4))/375.0/pow(beta1,4)*h_sep*log(1.0/h_sep);
          a_sh *= 6.0*MY_PI*mu*radi;*/
        } else {
          //a_sq = 6.0*MY_PI*mu*radi*(beta0*beta0/beta1/beta1/h_sep);
          a_sq = pre[0]*(beta[0][1]*beta[0][1]/(beta[1][1]*beta[1][1]*h_sep));
        }
 return {a_sq,a_sh};}
std::array<double,2> old_poly(double radi,double radj,double gap,double mu,int flaglog){
 double a_sq=0,a_sh=0,a_pu=0,h_sep=gap/radi;
 double beta0=radj/radi,beta1=1+beta0;
 double beta[2][5]={{0},{0}},pre[2];beta[0][1]=beta0;beta[1][1]=beta1;pre[1]=8.0*(pre[0]=MY_PI*mu*radi)*radi*radi;pre[0]*=6.0;
        if (flaglog) {
          // Jeffrey & Onishi (1984) near-field resistance functions are
          // expansions in the *symmetric* dimensionless gap
          // xi = 2*gap/(radi+radj) = 2*h_sep/beta1, not in the per-particle
          // h_sep = gap/radi.  Using h_sep in the log-order terms (and dropping
          // the beta0 prefactor of the squeeze log term) makes the resistance
          // depend on which particle is "i", breaking pairwise force symmetry
          // (Newton's 3rd law) for unequal radii.  See GitHub issue #1933.
          double xi = 2.0*h_sep/beta1;
          a_sq = beta0*beta0/beta1/beta1/h_sep +
            beta0*(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3.0)*log(1.0/xi);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0 *
                   pow(beta0,3.0)+pow(beta0,4.0))/21.0/pow(beta1,4.0) *
            h_sep*log(1.0/xi);
          a_sq *= 6.0*MY_PI*mu*radi;
          a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/pow(beta1,3.0) *
            log(1.0/xi);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3.0) +
                       16.0*pow(beta0,4.0))/375.0/pow(beta1,4.0) *
            h_sep*log(1.0/xi);
          a_sh *= 6.0*MY_PI*mu*radi;
          // old invalid eq for pumping term
          // changed 29Jul16 from eq 9.25 -> 9.27 in Kim and Karilla
//          a_pu = beta0*(4.0+beta0)/10.0/beta1/beta1*log(1.0/h_sep);
//          a_pu += (32.0-33.0*beta0+83.0*beta0*beta0+43.0 *
//                   pow(beta0,3.0))/250.0/pow(beta1,3.0)*h_sep*log(1.0/h_sep);
//          a_pu *= 8.0*MY_PI*mu*pow(radi,3.0);
          a_pu = 2.0*beta0/5.0/beta1*log(1.0/xi);
          a_pu += 2.0*(8.0+6.0*beta0+33.0*beta0*beta0)/125.0/beta1/beta1*
                   h_sep*log(1.0/xi);
          a_pu *= 8.0*MY_PI*mu*pow(radi,3.0);
        } else a_sq = 6.0*MY_PI*mu*radi*(beta0*beta0/beta1/beta1/h_sep);
 return {a_sq,a_sh};}
int main(int argc,char**argv){std::ifstream in(argv[1]);std::ofstream out(argv[2]);out<<std::setprecision(17);out<<"id,flaglog,target_poly_sq,target_poly_sq_swap,target_poly_sh,target_poly_sh_swap,target_upoly_sq,target_upoly_sq_swap,target_upoly_sh,target_upoly_sh_swap,old_poly_sq,old_poly_sq_swap,old_poly_sh,old_poly_sh_swap\n";long id;double a,b,h,mu;
while(in>>id>>a>>b>>h>>mu)for(int f=0;f<=1;++f){auto x=target_poly(a,b,h,mu,f),y=target_poly(b,a,h,mu,f),u=target_upoly(a,b,h,mu,f),v=target_upoly(b,a,h,mu,f),o=old_poly(a,b,h,mu,f),p=old_poly(b,a,h,mu,f);out<<id<<","<<f<<","<<x[0]<<","<<y[0]<<","<<x[1]<<","<<y[1]<<","<<u[0]<<","<<v[0]<<","<<u[1]<<","<<v[1]<<","<<o[0]<<","<<p[0]<<","<<o[1]<<","<<p[1]<<"\n";}}
