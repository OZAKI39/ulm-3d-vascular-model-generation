#include "kinematic_wall_constraint.hpp"
#include <cmath>
#include <limits>
namespace workflow {
using namespace rigid;
WallConstraintResult constrain_surface_wall(const WallConstraintInput&i,const passive::WallDistance&w,const SurfacePatchTopology&topology,double factor,double theta){
 WallConstraintResult r;r.raw_velocity=i.raw_velocity;r.constrained_velocity=i.raw_velocity;r.endpoint=add(i.base,mul(i.raw_velocity,i.stage_dt));r.contact_point=i.evaluation;r.gap=w.distance(i.evaluation)-i.radius;
 r.raw_segment_safe=w.segment_safe_exact(i.base,r.endpoint,i.radius+i.margin);
 if(!r.raw_segment_safe){
  double lo=0,hi=1;Vec rawend=r.endpoint;
  // Prefix safety is monotone even when the raw chord leaves and re-enters.
  if(w.segment_safe_exact(i.base,i.base,i.radius+i.margin))for(int k=0;k<40;k++){double mid=(lo+hi)/2;Vec x=add(i.base,mul(sub(rawend,i.base),mid));if(w.segment_safe_exact(i.base,x,i.radius+i.margin))lo=mid;else hi=mid;}
  else hi=0;
  r.first_contact_fraction=hi;r.contact_point=add(i.base,mul(sub(rawend,i.base),hi));
 }
 r.surface=topology.query(r.contact_point,factor,theta);
 if(r.surface.status=="FAIL_NONMANIFOLD"||r.surface.status.find("STOP_")==0){r.mode=r.surface.status;r.status=r.surface.status;return r;}
 r.normal=r.surface.clusters.at(0).normal;r.raw_normal_velocity=dot(i.raw_velocity,r.normal);
 if(!r.raw_segment_safe){std::vector<Vec>normals;for(auto const&c:r.surface.clusters)normals.push_back(c.normal);r.projection=project_surface_velocity(i.raw_velocity,normals);r.constrained_velocity=r.projection.used;r.active=norm(sub(r.constrained_velocity,i.raw_velocity))>0;r.removed_normal_velocity=norm(sub(r.constrained_velocity,i.raw_velocity));r.mode=r.surface.status;r.status="SURFACE_KINEMATIC_PROJECTION";}
 r.endpoint=add(i.base,mul(r.constrained_velocity,i.stage_dt));r.uncorrected_endpoint=r.endpoint;r.safe=w.segment_safe_exact(i.base,r.endpoint,i.radius+i.margin);
 if(!r.safe){
  r.status="CURVED_GEOMETRY_RETRY";r.mode="TIME_REFINEMENT";r.endpoint_surface=topology.query(r.endpoint,factor,theta);r.gap_before=r.endpoint_surface.minimum_distance-i.radius;r.gap_after=r.gap_before;
  if(r.endpoint_surface.status=="FAIL_NONMANIFOLD"||r.endpoint_surface.status.find("STOP_")==0){r.mode=r.endpoint_surface.status;r.status=r.endpoint_surface.status;return r;}
  double swept_gap=w.segment_distance_exact(i.base,r.endpoint)-i.radius;
  double deficit=i.margin-std::min(r.gap_before,swept_gap),cap=i.margin/4,cushion=64*std::numeric_limits<double>::epsilon()*std::max(norm(r.endpoint),i.radius+i.margin);
  // Multi-surface contact always refines time; no position correction there.
  if(r.active&&i.allow_offset&&r.surface.clusters.size()==1&&r.endpoint_surface.clusters.size()==1&&deficit>0){
   if(deficit+cushion>cap){r.status="DEEP_PENETRATION_REJECTED";return r;}
   Vec d=sub(r.endpoint,r.endpoint_surface.clusters[0].point);r.correction_normal=mul(d,1/norm(d));r.position_projection_distance=deficit+cushion;r.position_projection_used=true;r.endpoint=add(r.endpoint,mul(r.correction_normal,r.position_projection_distance));r.gap_after=w.distance(r.endpoint)-i.radius;
   r.safe=w.segment_safe_exact(i.base,r.endpoint,i.radius+i.margin);r.status=r.safe?"BOUNDED_POSITION_PROJECTION":"CURVED_GEOMETRY_RETRY";r.mode=r.safe?"BOUNDED_POSITION_CORRECTION":"TIME_REFINEMENT";
  }
 }
 return r;
}
}
