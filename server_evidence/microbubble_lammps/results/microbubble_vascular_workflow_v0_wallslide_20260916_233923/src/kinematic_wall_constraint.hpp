#pragma once
#include "rigid_math.hpp"
#include "wall_distance.hpp"
#include <functional>
namespace workflow {
using rigid::Vec;
struct WallConstraintInput {
 Vec evaluation,base,raw_velocity;
 double radius,margin,stage_dt;
 bool allow_offset=false;
};
struct WallConstraintResult {
 bool active=false,raw_segment_safe=false,safe=false;
 double gap=0,raw_normal_velocity=0,removed_normal_velocity=0;
 Vec normal{},raw_velocity{},constrained_velocity{},endpoint{};
 bool position_projection_used=false;
 double position_projection_distance=0,gap_before=0,gap_after=0;
 Vec correction_normal{},uncorrected_endpoint{};
 std::string status="NO_CONTACT_ACTION";
};
using WallQuery=std::function<passive::ClosestWallPoint(Vec)>;
using WallSegment=std::function<bool(Vec,Vec,double)>;
WallConstraintResult constrain_wall(const WallConstraintInput&,const WallQuery&,const WallSegment&);
WallConstraintResult constrain_wall(const WallConstraintInput&,const passive::WallDistance&);
}
