#include "frozen_flow.hpp"
#include <iomanip>
#include <iostream>
int main(int argc,char**argv){
  if(argc!=2)return 2;
  try{
    frozen::FlowField field(argv[1]);frozen::FlowFieldSampler sampler(field);
    frozen::Vec p;std::cout<<std::setprecision(17);
    while(std::cin>>p[0]>>p[1]>>p[2]){
      auto q=sampler.query(p);
      std::cout<<static_cast<int>(q.status)<<' '<<q.u[0]<<' '<<q.u[1]<<' '<<q.u[2]<<'\n';
    }
    return std::cin.eof()?0:3;
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
