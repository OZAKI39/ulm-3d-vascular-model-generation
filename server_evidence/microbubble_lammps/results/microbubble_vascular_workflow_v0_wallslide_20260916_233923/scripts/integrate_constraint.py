from pathlib import Path
S=Path(__file__).resolve().parents[1];p=S/'src/workflow_main.cpp';s=p.read_text()
s=s.replace('#include "lifecycle.hpp"','#include "lifecycle.hpp"\n#include "kinematic_wall_constraint.hpp"')
s=s.replace('retrylog,pairlog;','retrylog,pairlog,walllog,projectionlog;')
s=s.replace('   traj.open(','''   walllog.open("wall_constraint_timeseries.csv");projectionlog.open("GEOMETRIC_PROJECTION_EVENTS.csv");
   walllog<<"time,step,stage,particle_id,gap,gap_over_radius,nx,ny,nz,vx_raw,vy_raw,vz_raw,vx_used,vy_used,vz_used,vn_raw,vn_used,tangential_speed_raw,tangential_speed_used,removed_normal_speed,constraint_status,dt,wall_constraint_active,raw_segment_safe,omega_raw_x,omega_raw_y,omega_raw_z,omega_used_x,omega_used_y,omega_used_z,x_eval,y_eval,z_eval,disclaimer\\n";
   projectionlog<<"particle_id,time,step,stage,attempt,accepted,gap_before,gap_after,correction_distance,nx,ny,nz,before_x,before_y,before_z,after_x,after_y,after_z,reason,disclaimer\\n";
   traj.open(''')
s=s.replace('&retrylog,&pairlog}', '&retrylog,&pairlog,&walllog,&projectionlog}')
pos=s.index('  auto background=')
s=s[:pos]+r'''  auto write_wall=[&](int step,int stage,double t,double dt,const std::vector<Particle>&p,const Solution&raw,const Solution&used,const std::vector<workflow::WallConstraintResult>&wr){if(rank)return;for(size_t i=0;i<p.size();i++){auto const&r=wr[i];walllog<<t<<','<<step<<','<<stage<<','<<p[i].id<<','<<r.gap<<','<<r.gap/p[i].a;for(auto v:{r.normal,r.raw_velocity,r.constrained_velocity})for(double x:v)walllog<<','<<x;double vn=dot(r.constrained_velocity,r.normal);walllog<<','<<r.raw_normal_velocity<<','<<vn<<','<<norm(sub(r.raw_velocity,mul(r.normal,r.raw_normal_velocity)))<<','<<norm(sub(r.constrained_velocity,mul(r.normal,vn)))<<','<<r.removed_normal_velocity<<','<<r.status<<','<<dt<<','<<r.active<<','<<r.raw_segment_safe;for(int j=3;j<6;j++)walllog<<','<<raw.q[6*i+j];for(int j=3;j<6;j++)walllog<<','<<used.q[6*i+j];for(double x:p[i].x)walllog<<','<<x;walllog<<','<<workflow::disclaimer<<'
';}};
  auto write_projection=[&](int step,int stage,int attempt,bool accepted,double t,const std::vector<Particle>&p,const std::vector<workflow::WallConstraintResult>&wr){if(rank)return;for(size_t i=0;i<wr.size();i++){auto const&r=wr[i];if(!r.position_projection_used)continue;projectionlog<<p[i].id<<','<<t<<','<<step<<','<<stage<<','<<attempt<<','<<accepted<<','<<r.gap_before<<','<<r.gap_after<<','<<r.position_projection_distance;for(auto v:{r.correction_normal,r.uncorrected_endpoint,r.endpoint})for(double x:v)projectionlog<<','<<x;projectionlog<<",BOUNDED_GEOMETRIC_OFFSET_PROJECTION,"<<workflow::disclaimer<<'
';}};
  auto constrain=[&](const std::vector<Particle>&eval,const std::vector<Particle>&base,const Solution&raw,Solution&used,std::vector<Particle>&end,std::vector<workflow::WallConstraintResult>&wr,double h,bool offset){used=raw;end=base;wr.clear();bool safe=true;for(size_t i=0;i<eval.size();i++){Vec v{raw.q[6*i],raw.q[6*i+1],raw.q[6*i+2]};auto r=workflow::constrain_wall({eval[i].x,base[i].x,v,eval[i].a,margin,h,offset},geometry.wall);wr.push_back(r);for(int j=0;j<3;j++)used.q[6*i+j]=r.constrained_velocity[j];end[i].x=r.endpoint;safe=safe&&r.safe;}return safe;};
''' +s[pos:]
s=s.replace('Solution first,second;','Solution first,second,first_used,second_used;std::vector<workflow::WallConstraintResult>wallfirst,wallsecond;')
s=s.replace('bool ok=false;int retries=0;','bool ok=false;int retries=0;bool wall_retry=false;')
s=s.replace('    std::vector<workflow::Admission> admissions;', '    wallfirst.clear();wallsecond.clear();\n    std::vector<workflow::Admission> admissions;')
a=s.index('      first=solve(');b=s.index('       else{\n        crossings.clear()',a)
s=s[:a]+'''      first=solve(p,bg,pairs,&p,.5*dt);
      bool safe=constrain(p,p,first,first_used,mid,wallfirst,.5*dt,wall_retry);
      if(!safe)rejection="CURVED_GEOMETRY_MIDPOINT_RETRY";
      else if(swept(p,mid,pairs)<-1e-12)rejection="PAIR_SWEPT_MIDPOINT_REJECT";
      else{
       second=solve(mid,background(mid),pairs,&p,dt);
       safe=constrain(mid,p,second,second_used,next,wallsecond,dt,wall_retry);
       sg=swept(p,next,pairs);if(!safe)rejection="CURVED_GEOMETRY_FINAL_RETRY";else if(sg< -1e-12)rejection="PAIR_SWEPT_REJECT";
'''+s[b:]
s=s.replace('    if(ok)break;','''    write_projection(accepted,0,retries,ok,time,p,wallfirst);write_projection(accepted,1,retries,ok,time+.5*dt,mid,wallsecond);
    if(ok)break;
    if(rejection.find("CURVED_GEOMETRY")!=std::string::npos)wall_retry=true;''')
s=s.replace('"WALL_CONSTRAINT_STALL");dt=nextdt;', '"CURVED_GEOMETRY_TANGENTIAL_STALL");dt=nextdt;')
s=s.replace('   // Admission starts', '''   write_wall(accepted,0,time,dt,p,first,first_used,wallfirst);write_wall(accepted,1,time+.5*dt,dt,mid,second,second_used,wallsecond);
   // Admission starts''')
s=s.replace('record(accepted,time,p,first,{})','record(accepted,time,p,first_used,{})').replace('a[j+3]=second.q[6*i+j]','a[j+3]=second_used.q[6*i+j]').replace('record(accepted,time,actual,second,crossings)','record(accepted,time,actual,second_used,crossings)')
p.write_text(s)
p=S/'CMakeLists.txt';s=p.read_text().replace('project(BCFluxWorkflowV0','project(WallSlidingWorkflowV0').replace('src/wall_distance.cpp)','src/wall_distance.cpp src/kinematic_wall_constraint.cpp)');s+='\nadd_executable(test_wall_constraint src/test_wall_constraint.cpp)\ntarget_compile_options(test_wall_constraint PRIVATE -UNDEBUG)\ntarget_link_libraries(test_wall_constraint PRIVATE rigid_math)\n';p.write_text(s)
