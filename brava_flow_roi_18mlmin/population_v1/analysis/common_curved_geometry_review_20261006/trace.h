// Observation only: emitted in separate trace builds, never a production hook.
#pragma once
#include "fs.h"
#include <fstream>
#include <iomanip>
namespace validation_trace {
inline std::ofstream points, parameters;
inline void sample(int opt,int g,const std::array<fsType,2>& fs,double jac,
                   const Array<double>& metric,const Array<double>& Dw,
                   const Array<double>& Dp,const Array<double>& Hw) {
  if(!points.is_open())return;
  points<<opt<<' '<<g<<' '<<fs[0].eNoN<<' '<<fs[1].eNoN<<' '
        <<fs[opt-1].w(g)<<' '<<jac<<' ';
  for(int j=0;j<3;j++)points<<fs[0].xi(j,g)<<' ';
  for(int a=0;a<fs[0].eNoN;a++)points<<fs[0].N(a,g)<<' ';
  for(int a=0;a<fs[1].eNoN;a++)points<<fs[1].N(a,g)<<' ';
  for(int a=0;a<Dw.ncols();a++)for(int j=0;j<3;j++)points<<Dw(j,a)<<' ';
  for(int a=0;a<Dp.ncols();a++)for(int j=0;j<3;j++)points<<Dp(j,a)<<' ';
  for(int a=0;a<Hw.ncols();a++)for(int j=0;j<6;j++)points<<Hw(j,a)<<' ';
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)points<<metric(i,j)<<' ';
  points<<'\n';
}
inline void tau(int opt,int vms,double m,double c,double b) {
  if(parameters.is_open())parameters<<opt<<' '<<vms<<' '<<m<<' '<<c<<' '<<b<<'\n';
}
}
