#pragma once
#include "reduced_wall_types.hpp"
#include <string>
#include <vector>

namespace reducedwall {
class ReducedWallKernel {
public:
    // Read-only frozen RMBW table. No production engine dependency.
    explicit ReducedWallKernel(const std::string& table_path);
    ReducedWallResistance evaluate_resistance(double radius_m, double viscosity_Pa_s,
                                             double gap_m, const WallFrame& frame) const;
    ReducedWallShearContribution evaluate_shear_rhs(const ReducedWallInput& input) const;
private:
    std::vector<double> epsilon_, log_epsilon_;
    std::vector<Matrix6> excess_;
    Matrix6 interpolate_scaled_excess(double epsilon) const;
};
} // namespace reducedwall
