#pragma once
#include "fix.h"
#include "frozen_flow.hpp"
#include <fstream>
#include <memory>
#include <unordered_map>
namespace LAMMPS_NS {
class FixFrozenFlowDrag:public Fix{
public:
 FixFrozenFlowDrag(LAMMPS*,int,char**);
 int setmask()override;
 void setup(int)override;
 void post_force(int)override;
 void end_of_step()override;
 void post_run()override;
private:
 std::unique_ptr<frozen::FlowField> field;
 double mu,wall_at_spawn,margin;frozen::Vec spawn;
 int stride,rank;std::string prefix;std::ofstream trace;
 bigint last_written=-1;long long calls=0,queries=0,invalid=0,nonfinite=0;
 double minimum_clearance=1e300;
 struct Applied{frozen::Vec force,velocity;};std::unordered_map<tagint,Applied> applied;
 void record();
};
}
