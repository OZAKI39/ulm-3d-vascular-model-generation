#include "wall_distance.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
int main(int argc,char**argv){if(argc!=4)return 2;passive::WallDistance wall(argv[1]);std::ifstream in(argv[2]);std::ofstream out(argv[3]);std::string row;std::getline(in,row);out<<std::setprecision(17)<<"id,distance_m\n";while(std::getline(in,row)){std::replace(row.begin(),row.end(),',',' ');std::istringstream s(row);double id,ref;frozen::Vec p;s>>id>>p[0]>>p[1]>>p[2]>>ref;out<<id<<','<<wall.distance(p)<<'\n';}}
