#include <iostream>
#include <tuple>
static_assert(__cplusplus==201703L);
static_assert(__GNUC__==12);
int main(){auto [a,b]=std::tuple<int,int>{19,23};std::cout<<"CXX17_OK gcc="<<__GNUC__<<"."<<__GNUC_MINOR__<<" result="<<a+b<<"\n";return a+b!=42;}
