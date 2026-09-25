#include "kinematic_wall_constraint.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace workflow {
using namespace rigid;
WallConstraintResult constrain_wall(const WallConstraintInput&i,const WallQuery&query,const WallSegment&segment){
 WallConstraintResult r;r.raw_velocity=i.raw_velocity;r.constrained_velocity=i.raw_velocity;
 auto q=query(i.evaluation);if(!(q.distance>0)||!std::isfinite(q.distance))throw std::runtime_error("INVALID_WALL_QUERY");
 r.normal=mul(sub(i.evaluation,q.point),1/q.distance);r.gap=q.distance-i.radius;r.raw_normal_velocity=dot(i.raw_velocity,r.normal);
 r.endpoint=add(i.base,mul(i.raw_velocity,i.stage_dt));r.raw_segment_safe=segment(i.base,r.endpoint,i.radius+i.margin);
 if(!r.raw_segment_safe&&r.raw_normal_velocity<0){
  if(q.ambiguous)throw std::runtime_error("STOP_WALL_NORMAL_AMBIGUITY");
  r.active=true;r.removed_normal_velocity=-r.raw_normal_velocity;r.constrained_velocity=sub(i.raw_velocity,mul(r.normal,r.raw_normal_velocity));r.status="INWARD_COMPONENT_REMOVED";
 }else if(r.raw_normal_velocity>0)r.status="OUTWARD_MOTION_UNCHANGED";
 else if(r.raw_normal_velocity==0)r.status="TANGENTIAL_MOTION_UNCHANGED";
 r.endpoint=add(i.base,mul(r.constrained_velocity,i.stage_dt));r.uncorrected_endpoint=r.endpoint;
 r.safe=segment(i.base,r.endpoint,i.radius+i.margin);
 if(!r.safe){
  r.status="CURVED_GEOMETRY_RETRY";
  auto e=query(r.endpoint);r.gap_before=e.distance-i.radius;r.gap_after=r.gap_before;
  double deficit=i.margin-r.gap_before,cap=i.margin/4;
  double cushion=64*std::numeric_limits<double>::epsilon()*std::max(norm(r.endpoint),i.radius+i.margin);
  if(r.active&&deficit>0&&i.allow_offset){
   if(deficit+cushion>cap){r.status="DEEP_PENETRATION_REJECTED";return r;}
   if(e.ambiguous)throw std::runtime_error("STOP_WALL_NORMAL_AMBIGUITY");
   if(!(e.distance>0))throw std::runtime_error("INVALID_WALL_QUERY");
   r.correction_normal=mul(sub(r.endpoint,e.point),1/e.distance);
   r.position_projection_distance=deficit+cushion;r.position_projection_used=true;
   r.endpoint=add(r.endpoint,mul(r.correction_normal,r.position_projection_distance));
   r.gap_after=query(r.endpoint).distance-i.radius;
   r.safe=segment(i.base,r.endpoint,i.radius+i.margin);
   r.status=r.safe?"BOUNDED_POSITION_PROJECTION":"CURVED_GEOMETRY_RETRY";
  }
 }
 return r;
}
WallConstraintResult constrain_wall(const WallConstraintInput&i,const passive::WallDistance&w){return constrain_wall(i,[&](Vec x){return w.closest(x);},[&](Vec a,Vec b,double d){return w.segment_safe_exact(a,b,d);});}
}
