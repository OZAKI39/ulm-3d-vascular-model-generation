#pragma once
#include <array>

namespace reducedwall {
using Vec3=std::array<double,3>;
using Vec6=std::array<double,6>;
using Matrix6=std::array<double,36>; // row-major; [Vx,Vy,Vz,Ox,Oy,Oz]
struct WallFrame { Vec3 t1, t2, n; }; // right handed; n points wall -> center
enum class WallShearSource { ANALYTIC=0, GRADIENT_PROXY=1, PALABOS_VALIDATED=2 };
enum class Status { OK, INVALID_GAP, INVALID_INPUT, INVALID_FRAME,
                    OUTSIDE_RESISTANCE_DOMAIN, OUTSIDE_SHEAR_DOMAIN,
                    UNSUPPORTED_SOURCE, NUMERICAL_FAILURE };
const char* status_name(Status s);
struct ReducedWallInput {
    double radius_m, viscosity_Pa_s, gap_m;
    WallFrame frame;
    Vec3 g_wall; // inverse seconds; explicitly projected onto tangent plane
    WallShearSource source=WallShearSource::ANALYTIC;
};
struct ReducedWallResistance {
    Status status=Status::INVALID_INPUT;
    double epsilon=0;
    Vec6 bulk_diagonal{};
    Matrix6 wall_excess_SI{}, total_SI{};
};
struct ReducedWallShearContribution {
    Status status=Status::INVALID_INPUT;
    ReducedWallResistance resistance;
    double FU=0, FOMEGA=0;
    Vec3 g_t{};
    Vec6 q_free{}, q_target{}, b_wall_shear{};
};
} // namespace reducedwall
