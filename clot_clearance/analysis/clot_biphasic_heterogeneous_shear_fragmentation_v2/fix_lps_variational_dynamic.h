// Internal isolated research coordinator; GPL-compatible LAMMPS extension.
#ifdef FIX_CLASS
FixStyle(LPS/VARIATIONAL/DYNAMIC,FixLPSVariationalDynamic);
#else
#ifndef LMP_FIX_LPS_VARIATIONAL_DYNAMIC_H
#define LMP_FIX_LPS_VARIATIONAL_DYNAMIC_H
#include "fix.h"
namespace LAMMPS_NS {
class BondBPMPeriVariational;
class FixLPSVariationalDynamic : public Fix {
 public:
  FixLPSVariationalDynamic(class LAMMPS *, int, char **);
  ~FixLPSVariationalDynamic() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_integrate() override;
  void post_force(int) override;
 private:
  BondBPMPeriVariational *material;
  double rho;
};
}
#endif
#endif
