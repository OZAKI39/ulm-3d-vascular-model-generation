#include "frozen_flow.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <memory>
using namespace frozen;
int main(int argc,char**argv){try{if(argc!=5)throw std::runtime_error("usage: flow_audit_cli query|force field.h5 input.csv output.csv");std::string mode=argv[1];std::unique_ptr<FlowField> field;if(mode=="query")field=std::make_unique<FlowField>(argv[2]);std::ifstream input(argv[3]);std::ofstream out(argv[4]);if(!input||!out)throw std::runtime_error("Cannot open audit CSV");out<<std::setprecision(17);std::string line;std::getline(input,line);out<<(mode=="query"?"id,status,ux,uy,uz\n":"id,fx,fy,fz\n");while(std::getline(input,line)){if(line.empty())continue;std::stringstream ss(line);std::string token;std::vector<double> a;while(std::getline(ss,token,','))a.push_back(std::stod(token));out<<int64_t(a.at(0));if(mode=="query"){auto q=FlowFieldSampler(*field).query({a.at(1),a.at(2),a.at(3)});out<<','<<int(q.status);for(double v:q.u)out<<','<<v;}else{auto f=TechnicalStokesDrag::force(.001,a.at(1),{a.at(2),a.at(3),a.at(4)},{a.at(5),a.at(6),a.at(7)});for(double v:f)out<<','<<v;}out<<'\n';}return 0;}catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 1;}}
