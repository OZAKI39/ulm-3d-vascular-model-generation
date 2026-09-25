#pragma once
#include "rigid_math.hpp"
#include "reduced_wall_kernel.hpp"

namespace phase2 {
struct Plane {
    rigid::Vec origin;
    reducedwall::WallFrame frame;
    double gamma_dot;
};
// A plane has an affine signed distance: endpoint checks certify its entire segment.
// Safety neither queries the hydrodynamic kernel nor changes a proposed position.
double center_distance(const Plane&, rigid::Vec);
double gap(const Plane&, rigid::Vec, double radius);
bool segment_safe(const Plane&, rigid::Vec start, rigid::Vec end, double radius);
const char* engine_status(reducedwall::Status);
class FlatWallEngine {
public:
    FlatWallEngine(const std::string& table, Plane plane);
    rigid::Background background(const rigid::Particle&) const;
    reducedwall::ReducedWallShearContribution evaluate(const rigid::Particle&) const;
    const Plane& plane() const { return plane_; }
private:
    reducedwall::ReducedWallKernel kernel_;
    Plane plane_;
};
}
