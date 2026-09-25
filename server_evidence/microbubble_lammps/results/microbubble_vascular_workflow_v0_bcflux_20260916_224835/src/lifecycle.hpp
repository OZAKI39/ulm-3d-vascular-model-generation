#pragma once
#include "lammps.h"
#include "input.h"
#include <vector>
#include <string>
namespace workflow {
// Shared actual production removal path. compress no preserves all remaining
// globally assigned tags and prevents reuse after a deletion or trial rollback.
inline void remove_atoms(LAMMPS_NS::LAMMPS*l,const std::vector<long long>&ids){
 if(ids.empty())return;std::string cmd="group workflow_remove id";for(auto id:ids)cmd+=" "+std::to_string(id);
 l->input->one(cmd);l->input->one("delete_atoms group workflow_remove compress no");l->input->one("group workflow_remove delete");l->input->one("run 0 post no");
}
}
