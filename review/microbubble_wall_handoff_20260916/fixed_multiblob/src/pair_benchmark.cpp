#include "pair_baseline_copy/rigid_math.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
int main(){double a=9.683592065545495e-7;std::vector<rigid::Particle> p={{0,0,{0,0,0},a},{1,0,{0,0,2.1*a},a}};std::vector<rigid::Background> bg(2);std::vector<std::pair<int,int>> pairs={{0,1}};volatile double checksum=0;auto warm=rigid::assemble(p,bg,pairs);std::cout<<std::setprecision(17)<<"{\"queries_per_repeat\":10000,\"seconds_per_query\":[";for(int r=0;r<5;r++){auto t=std::chrono::steady_clock::now();for(int k=0;k<10000;k++){auto s=rigid::assemble(p,bg,pairs);for(int i=0;i<6;i++)for(int j=0;j<6;j++)checksum+=s.R[i*12+j];}double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();std::cout<<(r?",":"")<<secs/10000;}std::cout<<"],\"checksum\":"<<checksum<<",\"R12_row_major\":[";for(int k=0;k<144;k++)std::cout<<(k?",":"")<<warm.R[k];std::cout<<"]}\n";}
