"""Apply a small additive interface change to this stage's stable baseline copies."""
from pathlib import Path
r=Path(__file__).resolve().parents[1]
h=(r/'provenance/baseline_src/rigid_math.hpp').read_text().replace('#include "frozen_flow.hpp"','#include "frozen_flow.hpp"\n#include "reduced_wall_types.hpp"')
h=h.replace('>> &','>>&')
h=h.replace('const std::vector<std::pair<int,int>>&pairs);','const std::vector<std::pair<int,int>>&pairs,const std::vector<reducedwall::ReducedWallShearContribution>*walls=nullptr);')
h=h.replace('double dt=0);','double dt=0,const std::vector<reducedwall::ReducedWallShearContribution>*walls=nullptr);')
cpp=(r/'provenance/baseline_src/rigid_math.cpp').read_text()
cpp=cpp.replace('const std::vector<std::pair<int,int>>&pairs){\n int d=', 'const std::vector<std::pair<int,int>>&pairs,const std::vector<reducedwall::ReducedWallShearContribution>*walls){\n int d=',1)
addition=''' // Phase 2: preserve the existing bulk and pair RHS; add only wall excess once.
 if(walls){
  if(walls->size()!=p.size())throw std::runtime_error("WALL_LOAD_SIZE_MISMATCH");
  for(size_t i=0;i<p.size();++i){
   const auto&w=(*walls)[i];
   if(w.status!=reducedwall::Status::OK)throw std::runtime_error("INVALID_WALL_LOAD");
   for(int k=0;k<6;++k){
    s.rhs[6*i+k]+=w.b_wall_shear[k];
    for(int j=0;j<6;++j)s.R[(6*i+k)*d+6*i+j]+=w.resistance.wall_excess_SI[6*k+j];
   }
  }
 }
'''
assert cpp.count(' return s;\n}')==1
cpp=cpp.replace(' return s;\n}',addition+' return s;\n}',1)
cpp=cpp.replace('const std::vector<Particle>*base,double dt){\n auto sys=assemble(p,bg,pairs);','const std::vector<Particle>*base,double dt,const std::vector<reducedwall::ReducedWallShearContribution>*walls){\n auto sys=assemble(p,bg,pairs,walls);',1)
(r/'src/rigid_math.hpp').write_text(h);(r/'src/rigid_math.cpp').write_text(cpp)
print('Stable assembly modified additively; PCG, scaling, constraints and pair formulas preserved.')
