#include "flat_wall_engine.hpp"
#include <cmath>
#include <stdexcept>

namespace phase2 {
double center_distance(const Plane& plane, rigid::Vec x) {
    long double d=0;
    for(int k=0;k<3;++k)
        d+=(static_cast<long double>(x[k])-plane.origin[k])*plane.frame.n[k];
    return static_cast<double>(d);
}
double gap(const Plane& plane, rigid::Vec x, double a) {
    return center_distance(plane,x)-a;
}
bool segment_safe(const Plane& plane, rigid::Vec start, rigid::Vec end, double a) {
    if(!std::isfinite(a)||a<=0) return false;
    const double h0=gap(plane,start,a),h1=gap(plane,end,a);
    return std::isfinite(h0)&&std::isfinite(h1)&&h0>0&&h1>0;
}
const char* engine_status(reducedwall::Status status) {
    if(status==reducedwall::Status::OUTSIDE_RESISTANCE_DOMAIN||
       status==reducedwall::Status::OUTSIDE_SHEAR_DOMAIN)
        return "OUTSIDE_CERTIFIED_WALL_DOMAIN";
    return reducedwall::status_name(status);
}
FlatWallEngine::FlatWallEngine(const std::string& table,Plane plane)
    :kernel_(table),plane_(plane) {
    for(double v:plane.origin) if(!std::isfinite(v)) throw std::runtime_error("INVALID_PLANE");
    if(!std::isfinite(plane.gamma_dot)) throw std::runtime_error("INVALID_ANALYTIC_SHEAR");
    // Use the frozen frame validator. No local-plane validity or geometry inference.
    auto check=kernel_.evaluate_resistance(1.,.001,.01,plane.frame);
    if(check.status!=reducedwall::Status::OK) throw std::runtime_error(engine_status(check.status));
}
rigid::Background FlatWallEngine::background(const rigid::Particle& p) const {
    rigid::Background b;
    const auto g=rigid::mul(plane_.frame.t1,plane_.gamma_dot);
    b.u=rigid::mul(g,center_distance(plane_,p.x));
    b.omega=rigid::mul(rigid::cross(plane_.frame.n,g),.5);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
        b.grad[3*i+j]=g[i]*plane_.frame.n[j];
        b.strain[3*i+j]=.5*(g[i]*plane_.frame.n[j]+g[j]*plane_.frame.n[i]);
    }
    return b;
}
reducedwall::ReducedWallShearContribution FlatWallEngine::evaluate(const rigid::Particle& p) const {
    return kernel_.evaluate_shear_rhs({p.a,.001,gap(plane_,p.x,p.a),plane_.frame,
        rigid::mul(plane_.frame.t1,plane_.gamma_dot),reducedwall::WallShearSource::ANALYTIC});
}
}
