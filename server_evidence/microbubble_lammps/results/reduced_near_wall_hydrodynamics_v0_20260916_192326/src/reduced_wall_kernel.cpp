#include "reduced_wall_kernel.hpp"
#include "cf2003_coefficients_generated.hpp"
#include <hdf5.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace reducedwall {
namespace {
struct H5Handle {
    hid_t id; herr_t (*close)(hid_t);
    H5Handle(hid_t value,herr_t (*closer)(hid_t)):id(value),close(closer) {
        if(id<0) throw std::runtime_error("RMBW_HDF5_OPEN_ERROR");
    }
    ~H5Handle(){close(id);}
    H5Handle(const H5Handle&)=delete;
    H5Handle& operator=(const H5Handle&)=delete;
};
template<class A> bool finite(const A& a) {
    return std::all_of(a.begin(),a.end(),[](double v){return std::isfinite(v);});
}
double dot(const Vec3&a,const Vec3&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vec3 cross(const Vec3&a,const Vec3&b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
bool valid_frame(const WallFrame&f) {
    const Vec3 v[]={f.t1,f.t2,f.n};
    for(int i=0;i<3;++i){
        if(!finite(v[i]))return false;
        for(int j=0;j<3;++j)
            if(std::abs(dot(v[i],v[j])-(i==j?1.:0.))>1e-12)return false;
    }
    return std::abs(dot(cross(f.t1,f.t2),f.n)-1.)<=1e-12;
}
Matrix6 rotate(const Matrix6&a,const WallFrame&f) {
    Matrix6 r{};const Vec3 col[]={f.t1,f.t2,f.n};
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)
        for(int k=0;k<3;++k)for(int l=0;l<3;++l)
            r[6*i+j]+=col[k][i%3]*a[6*(3*(i/3)+k)+3*(j/3)+l]*col[l][j%3];
    return r;
}
double response(const std::array<double,20>&c,double x) {
    double p=c.back();
    for(int j=18;j>=0;--j)p=p*x+c[j];
    return 1./p;
}
}
const char* status_name(Status s) {
    switch(s){
        case Status::OK:return "OK";
        case Status::INVALID_GAP:return "INVALID_GAP";
        case Status::INVALID_INPUT:return "INVALID_INPUT";
        case Status::INVALID_FRAME:return "INVALID_FRAME";
        case Status::OUTSIDE_RESISTANCE_DOMAIN:return "OUTSIDE_RESISTANCE_DOMAIN";
        case Status::OUTSIDE_SHEAR_DOMAIN:return "OUTSIDE_SHEAR_DOMAIN";
        case Status::UNSUPPORTED_SOURCE:return "UNSUPPORTED_SOURCE";
        case Status::NUMERICAL_FAILURE:return "NUMERICAL_FAILURE";
    }
    return "INVALID_STATUS";
}
ReducedWallKernel::ReducedWallKernel(const std::string&path) {
    H5Handle file(H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
    H5Handle ds(H5Dopen2(file.id,"epsilon",H5P_DEFAULT),H5Dclose);
    H5Handle space(H5Dget_space(ds.id),H5Sclose);
    hsize_t n=0;
    if(H5Sget_simple_extent_ndims(space.id)!=1)throw std::runtime_error("INVALID_RMBW_GRID_SHAPE");
    H5Sget_simple_extent_dims(space.id,&n,nullptr);
    if(n!=8211)throw std::runtime_error("UNEXPECTED_RMBW_GRID_SIZE");
    epsilon_.resize(n);log_epsilon_.resize(n);excess_.resize(n);
    if(H5Dread(ds.id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,epsilon_.data())<0)
        throw std::runtime_error("RMBW_GRID_READ_FAILED");
    H5Handle md(H5Dopen2(file.id,"R_wall_excess_scaled",H5P_DEFAULT),H5Dclose);
    H5Handle ms(H5Dget_space(md.id),H5Sclose);hsize_t shape[3]{};
    if(H5Sget_simple_extent_ndims(ms.id)!=3)throw std::runtime_error("INVALID_RMBW_MATRIX_RANK");
    H5Sget_simple_extent_dims(ms.id,shape,nullptr);
    if(shape[0]!=n||shape[1]!=6||shape[2]!=6)throw std::runtime_error("INVALID_RMBW_MATRIX_SHAPE");
    if(H5Dread(md.id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,excess_.data())<0)
        throw std::runtime_error("RMBW_MATRIX_READ_FAILED");
    if(epsilon_.front()!=.001||epsilon_.back()!=20.)throw std::runtime_error("INVALID_RMBW_DOMAIN");
    for(size_t i=0;i<n;++i){
        if(!std::isfinite(epsilon_[i])||epsilon_[i]<=0||(i&&epsilon_[i]<=epsilon_[i-1])||!finite(excess_[i]))
            throw std::runtime_error("INVALID_RMBW_DATA");
        log_epsilon_[i]=std::log(epsilon_[i]);
    }
}
Matrix6 ReducedWallKernel::interpolate_scaled_excess(double e) const {
    auto it=std::lower_bound(epsilon_.begin(),epsilon_.end(),e);
    const size_t j=static_cast<size_t>(it-epsilon_.begin());
    if(*it==e)return excess_[j]; // caller already rejects all exterior values
    const size_t i=j-1;
    const double w=(std::log(e)-log_epsilon_[i])/(log_epsilon_[j]-log_epsilon_[i]);
    Matrix6 a{};
    for(int k=0;k<36;++k)a[k]=(1.-w)*excess_[i][k]+w*excess_[j][k];
    return a;
}
ReducedWallResistance ReducedWallKernel::evaluate_resistance(double a,double mu,double h,const WallFrame&f) const {
    ReducedWallResistance r;
    if(!std::isfinite(a)||!std::isfinite(mu)||!std::isfinite(h)||a<=0||mu<=0)return r;
    if(h<=0){r.status=Status::INVALID_GAP;return r;}
    if(!valid_frame(f)){r.status=Status::INVALID_FRAME;return r;}
    r.epsilon=h/a;
    if(!std::isfinite(r.epsilon)||r.epsilon<epsilon_.front()||r.epsilon>epsilon_.back()){
        r.status=Status::OUTSIDE_RESISTANCE_DOMAIN;return r;
    }
    const double pi=std::acos(-1.);
    const double tt=6*pi*mu*a, rr=8*pi*mu*a*a*a;
    if(!std::isfinite(tt)||!std::isfinite(rr)||tt<=0||rr<=0){r.status=Status::NUMERICAL_FAILURE;return r;}
    for(int i=0;i<6;++i)r.bulk_diagonal[i]=(i<3?tt:rr);
    auto local=interpolate_scaled_excess(r.epsilon);
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)
        local[6*i+j]*=std::sqrt(r.bulk_diagonal[i])*std::sqrt(r.bulk_diagonal[j]);
    r.wall_excess_SI=rotate(local,f);r.total_SI=r.wall_excess_SI;
    for(int i=0;i<6;++i)r.total_SI[6*i+i]+=r.bulk_diagonal[i];
    r.status=finite(r.total_SI)?Status::OK:Status::NUMERICAL_FAILURE;return r;
}
ReducedWallShearContribution ReducedWallKernel::evaluate_shear_rhs(const ReducedWallInput&in) const {
    ReducedWallShearContribution out;
    out.resistance=evaluate_resistance(in.radius_m,in.viscosity_Pa_s,in.gap_m,in.frame);
    out.status=out.resistance.status;if(out.status!=Status::OK)return out;
    if(in.source!=WallShearSource::ANALYTIC){out.status=Status::UNSUPPORTED_SOURCE;return out;}
    if(!finite(in.g_wall)){out.status=Status::INVALID_INPUT;return out;}
    const double e=out.resistance.epsilon;
    if(e<reference::epsilon_min||e>reference::epsilon_max){out.status=Status::OUTSIDE_SHEAR_DOMAIN;return out;}
    out.FU=response(reference::u,std::log(e));out.FOMEGA=response(reference::omega,std::log(e));
    const double d=in.radius_m+in.gap_m,gn=dot(in.g_wall,in.frame.n);
    for(int i=0;i<3;++i)out.g_t[i]=in.g_wall[i]-gn*in.frame.n[i];
    const auto c=cross(in.frame.n,out.g_t);
    for(int i=0;i<3;++i){
        out.q_free[i]=d*out.g_t[i];out.q_free[i+3]=.5*c[i];
        out.q_target[i]=out.FU*out.q_free[i];out.q_target[i+3]=out.FOMEGA*out.q_free[i+3];
    }
    for(int i=0;i<6;++i){
        out.b_wall_shear[i]=-out.resistance.bulk_diagonal[i]*out.q_free[i];
        for(int j=0;j<6;++j)out.b_wall_shear[i]+=out.resistance.total_SI[6*i+j]*out.q_target[j];
    }
    if(!std::isfinite(out.FU)||!std::isfinite(out.FOMEGA)||out.FU<=0||out.FU>=1||out.FOMEGA<=0||out.FOMEGA>=1
       ||!finite(out.b_wall_shear)||!finite(out.q_free)||!finite(out.q_target))out.status=Status::NUMERICAL_FAILURE;
    return out;
}
} // namespace reducedwall
