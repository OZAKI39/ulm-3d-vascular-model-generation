#pragma once
#include "fix.h"
#include "frozen_flow.hpp"
#include "wall_distance.hpp"
#include <memory>
#include <unordered_map>
#include <fstream>
namespace LAMMPS_NS {
class FixPassiveFlow:public Fix {
 public:
 FixPassiveFlow(LAMMPS*,int,char**);
 int setmask()override;
 void setup(int)override;
 void initial_integrate(int)override;
 void end_of_step()override;
 void post_run()override;
 private:
 std::unique_ptr<frozen::FlowField> field;
 std::unique_ptr<passive::WallDistance> wall;
 struct Proposal {frozen::Vec old,mid,next,velocity;double radius;};
 std::unordered_map<tagint,Proposal> pending;
 int stride,rank,check_pairs;double margin;
 bool stopped=false;int reason=0;
 long long queries=0,invalid=0,nonfinite=0,accepted=0;
 bigint last_written=-1;
 std::ofstream trace;
 std::string prefix;
 void prepare_next();void record(bool force=false);void stop(int);
 bool sample(frozen::Vec,frozen::Vec&);
};
}
