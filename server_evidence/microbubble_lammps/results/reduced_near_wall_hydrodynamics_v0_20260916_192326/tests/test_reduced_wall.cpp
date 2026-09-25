#include "reduced_wall_kernel.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
using namespace reducedwall;
namespace {
struct Gate { double limit,worst=0;size_t count=0,failures=0; };
std::map<std::string,Gate> gates;
void check(const std::string&name,double value) {
    auto&g=gates.at(name);++g.count;
    if(!std::isfinite(value)||value>g.limit)++g.failures;
    g.worst=std::max(g.worst,std::isfinite(value)?value:1e300);
}
double norm(const Vec6&v){double s=0;for(double x:v)s+=x*x;return std::sqrt(s);}
double mnorm(const Matrix6&a){double s=0;for(double x:a)s+=x*x;return std::sqrt(s);}
Vec6 scaled(const Vec6&v,const Vec6&bulk) {
    Vec6 z{};for(int i=0;i<6;++i)z[i]=v[i]/std::sqrt(bulk[i]);return z;
}
Matrix6 scaled(const Matrix6&m,const Vec6&bulk) {
    Matrix6 z{};for(int i=0;i<6;++i)for(int j=0;j<6;++j)z[6*i+j]=m[6*i+j]/std::sqrt(bulk[i])/std::sqrt(bulk[j]);return z;
}
double matrix_diff(const Matrix6&a,const Matrix6&b,const Vec6&bulk) {
    Matrix6 diff{};for(int i=0;i<36;++i)diff[i]=a[i]-b[i];
    return mnorm(scaled(diff,bulk))/std::max(1.,mnorm(scaled(b,bulk)));
}
double rhs_diff(const Vec6&a,const Vec6&b,const Vec6&bulk,double scale) {
    Vec6 v{};for(int i=0;i<6;++i)v[i]=a[i]-b[i];return norm(scaled(v,bulk))/scale;
}
Vec6 solve(const Matrix6&r,const Vec6&rhs,const Vec6&bulk) {
    // Test-only Gaussian elimination with partial pivoting. Work in the
    // conjugate scaled coordinates to avoid SI translation/rotation imbalance.
    auto m=scaled(r,bulk);auto b=scaled(rhs,bulk);
    for(int k=0;k<6;++k){
        int p=k;for(int i=k+1;i<6;++i)if(std::abs(m[6*i+k])>std::abs(m[6*p+k]))p=i;
        if(std::abs(m[6*p+k])<1e-15)throw std::runtime_error("SINGULAR_SOLVE");
        for(int j=0;j<6;++j)std::swap(m[6*k+j],m[6*p+j]);
        std::swap(b[k],b[p]);
        for(int i=k+1;i<6;++i){const double f=m[6*i+k]/m[6*k+k];
            for(int j=k;j<6;++j)m[6*i+j]-=f*m[6*k+j];
            b[i]-=f*b[k];}
    }
    Vec6 z{},q{};
    for(int i=5;i>=0;--i){double v=b[i];for(int j=i+1;j<6;++j)v-=m[6*i+j]*z[j];z[i]=v/m[6*i+i];}
    for(int i=0;i<6;++i)q[i]=z[i]/std::sqrt(bulk[i]);
    return q;
}
bool spd(const Matrix6&a) {
    Matrix6 l{};
    for(int i=0;i<6;++i)for(int j=0;j<=i;++j){
        double v=a[6*i+j];for(int k=0;k<j;++k)v-=l[6*i+k]*l[6*j+k];
        if(i==j){if(!std::isfinite(v)||v<=0)return false;l[6*i+j]=std::sqrt(v);}
        else l[6*i+j]=v/l[6*j+j];
    }
    // Explicit positive power for deterministic mixed translation/rotation probes.
    for(int k=1;k<=7;++k){Vec6 z{};for(int i=0;i<6;++i)z[i]=std::sin((k+1.)*(i+1.));
        double power=0;for(int i=0;i<6;++i)for(int j=0;j<6;++j)power+=z[i]*a[6*i+j]*z[j];
        if(!std::isfinite(power)||power<=0)return false;}
    return true;
}
ReducedWallShearContribution evaluate(const ReducedWallKernel&k,const ReducedWallInput&in) {
    auto out=k.evaluate_shear_rhs(in);if(out.status!=Status::OK)throw std::runtime_error(status_name(out.status));return out;
}
Vec6 run_tests(const ReducedWallKernel&k,const ReducedWallInput&in,const ReducedWallShearContribution&o,bool dense) {
    const auto&b=o.resistance.bulk_diagonal;const auto&r=o.resistance.total_SI;
    auto m=scaled(r,b);Matrix6 transpose{};
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)transpose[6*i+j]=r[6*j+i];
    check("A",matrix_diff(r,transpose,b));check("B",spd(m)?0:1);
    auto other=in;other.radius_m*=2;other.gap_m*=2;other.viscosity_Pa_s*=3;
    const auto rr=k.evaluate_resistance(other.radius_m,other.viscosity_Pa_s,other.gap_m,other.frame);
    if(rr.status!=Status::OK)throw std::runtime_error("SCALING_QUERY_FAILED");
    Matrix6 expected{};
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)
        expected[6*i+j]=o.resistance.wall_excess_SI[6*i+j]*(i<3?(j<3?6:12):(j<3?12:24));
    check("C",matrix_diff(rr.wall_excess_SI,expected,rr.bulk_diagonal));
    Vec6 characteristic{};
    for(int i=0;i<6;++i)characteristic[i]=std::sqrt(b[i])*o.q_free[i];
    double gnorm=0;for(double v:in.g_wall)gnorm+=v*v;
    double char_rhs=std::max(norm(characteristic),std::sqrt(b[0])*in.radius_m*std::sqrt(gnorm));
    if(char_rhs==0)char_rhs=std::sqrt(b[0])*in.radius_m; // unit shear for zero-shear probe
    for(int mode=0;mode<2;++mode){
        other=in;const double angle=mode?std::acos(-1.):.371;
        for(int i=0;i<3;++i){
            other.frame.t1[i]=std::cos(angle)*in.frame.t1[i]+std::sin(angle)*in.frame.t2[i];
            other.frame.t2[i]=-std::sin(angle)*in.frame.t1[i]+std::cos(angle)*in.frame.t2[i];
            if(mode){other.frame.t1[i]=-in.frame.t1[i];other.frame.t2[i]=-in.frame.t2[i];}
        }
        const auto v=evaluate(k,other);
        check(mode?"E":"D",std::max(matrix_diff(v.resistance.wall_excess_SI,o.resistance.wall_excess_SI,b),
                                    rhs_diff(v.b_wall_shear,o.b_wall_shear,b,char_rhs)));
    }
    other=in;other.g_wall={0,0,0};const auto zero=evaluate(k,other);
    check("F",norm(scaled(zero.b_wall_shear,b))/char_rhs);
    for(double factor:{-1.,2.}){
        other=in;for(int i=0;i<3;++i)other.g_wall[i]*=factor;
        const auto v=evaluate(k,other);Vec6 target{};
        for(int i=0;i<6;++i)target[i]=factor*o.b_wall_shear[i];
        check(factor<0?"G":"H",rhs_diff(v.b_wall_shear,target,b,char_rhs));
    }
    double fn=0;for(int i=0;i<3;++i)fn+=in.frame.n[i]*o.b_wall_shear[i];
    check("I",std::abs(fn)/std::sqrt(b[0])/char_rhs);
    Vec6 rhs{};for(int i=0;i<6;++i)rhs[i]=b[i]*o.q_free[i]+o.b_wall_shear[i];
    const Vec6 q=solve(r,rhs,b);Vec6 error{},target{};
    for(int i=0;i<6;++i){error[i]=std::sqrt(b[i])*(q[i]-o.q_target[i]);target[i]=std::sqrt(b[i])*o.q_target[i];}
    check("J",norm(error)/std::max(norm(target),char_rhs*1e-30));
    if(dense){
        const double us=q[0]/o.q_free[0],os=q[4]/o.q_free[4];
        check("K",std::max(std::abs(us/o.FU-1),std::abs(os/o.FOMEGA-1)));
    }
    return q;
}
std::vector<std::string> split(const std::string&s) {
    std::vector<std::string>v;std::stringstream stream(s);std::string x;
    while(std::getline(stream,x,','))v.push_back(x);
    return v;
}
template<class A> void values(std::ostream&f,const A&a){for(double x:a)f<<','<<x;}
}
int main(int argc,char**argv) try {
    if(argc!=5)throw std::runtime_error("usage: test_reduced_wall TABLE INPUT.csv OUTPUT.csv VALIDATION.json");
    for(const auto&x:std::vector<std::pair<std::string,double>>{{"A",1e-10},{"B",0},{"C",1e-10},{"D",1e-9},{"E",1e-9},
        {"F",1e-12},{"G",1e-10},{"H",1e-10},{"I",1e-12},{"J",1e-8},{"K",1e-8},{"INPUT_REJECTION",0}})gates.emplace(x.first,Gate{x.second});
    ReducedWallKernel kernel(argv[1]);std::ifstream input(argv[2]);std::ofstream output(argv[3]);
    if(!input||!output)throw std::runtime_error("FILE_OPEN_FAILED");
    output<<std::setprecision(17)<<"id,status,epsilon,FU,FOMEGA";
    for(const auto&x:std::vector<std::pair<std::string,int>>{{"bulk",6},{"E",36},{"R",36},{"qfree",6},{"qtarget",6},{"b",6},{"qsolved",6}})
        for(int i=0;i<x.second;++i)output<<','<<x.first<<i;
    output<<'\n';std::string line;std::getline(input,line);size_t valid=0,invalid=0;
    while(std::getline(input,line)){
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        if(line.empty())continue;
        auto row=split(line);if(row.size()!=19)throw std::runtime_error("BAD_INPUT_ROW");
        ReducedWallInput in{};in.radius_m=std::stod(row[1]);in.viscosity_Pa_s=std::stod(row[2]);in.gap_m=std::stod(row[3]);
        for(int i=0;i<3;++i){in.frame.t1[i]=std::stod(row[4+i]);in.frame.t2[i]=std::stod(row[7+i]);
            in.frame.n[i]=std::stod(row[10+i]);in.g_wall[i]=std::stod(row[13+i]);}
        in.source=static_cast<WallShearSource>(std::stoi(row[16]));auto out=kernel.evaluate_shear_rhs(in);
        check("INPUT_REJECTION",status_name(out.status)==row[17]?0:1);
        Vec6 solved{};
        if(out.status==Status::OK){++valid;solved=run_tests(kernel,in,out,row[18]=="dense");}else ++invalid;
        output<<row[0]<<','<<status_name(out.status)<<','<<out.resistance.epsilon<<','<<out.FU<<','<<out.FOMEGA;
        values(output,out.resistance.bulk_diagonal);values(output,out.resistance.wall_excess_SI);values(output,out.resistance.total_SI);
        values(output,out.q_free);values(output,out.q_target);values(output,out.b_wall_shear);values(output,solved);output<<'\n';
    }
    bool pass=true;for(const auto&kv:gates)pass=pass&&kv.second.count>0&&kv.second.failures==0;
    std::ofstream report(argv[4]);report<<std::setprecision(17)<<"{\n  \"ALGEBRA_PASS\": "<<(pass?"true":"false")
        <<",\n  \"valid_cases\": "<<valid<<",\n  \"invalid_cases\": "<<invalid<<",\n  \"tests\": {\n";
    bool first=true;
    for(const auto&kv:gates){if(!first)report<<",\n";first=false;const auto&g=kv.second;
        report<<"    \""<<kv.first<<"\": {\"count\": "<<g.count<<", \"failures\": "<<g.failures
              <<", \"worst\": "<<g.worst<<", \"limit\": "<<g.limit<<'}';}
    report<<"\n  }\n}\n";
    std::cout<<"CPP_A_K="<<(pass?"PASS":"FAIL")<<" valid="<<valid<<" invalid="<<invalid<<'\n';return pass?0:1;
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
