#pragma once
#include <fstream>
#include <string>
namespace dispatch_trace {
inline std::ofstream out;
inline std::string stage="setup";
inline void hit(const std::string& f,const std::string& face="*") {if(out.is_open())out<<stage<<','<<face<<','<<f<<'\n';}
}
