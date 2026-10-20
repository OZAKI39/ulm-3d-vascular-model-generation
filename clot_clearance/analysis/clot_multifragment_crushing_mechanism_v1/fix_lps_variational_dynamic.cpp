// Physical SI force conversion and one damage commit per native Verlet drift.
#include "fix_lps_variational_dynamic.h"
#include "bond_bpm_peri_variational.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "update.h"
#include "utils.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
using namespace FixConst;

FixLPSVariationalDynamic::FixLPSVariationalDynamic(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), material(nullptr), rho(0)
{
  if(narg!=4 || igroup!=0)
    error->all(FLERR,"LPS/VARIATIONAL/DYNAMIC requires group all and physical density");
  rho=utils::numeric(FLERR,arg[3],false,lmp);
  if(!(rho>0))error->all(FLERR,"Variational dynamic density must be positive");
  force_reneighbor=1;
  next_reneighbor=-1;
}
FixLPSVariationalDynamic::~FixLPSVariationalDynamic()
{
  if(material && force && force->bond==material)material->attach_dynamic(false);
}
int FixLPSVariationalDynamic::setmask()
{
  return POST_INTEGRATE | POST_FORCE;
}
void FixLPSVariationalDynamic::init()
{
  if(std::strcmp(update->integrate_style,"verlet")!=0)
    error->all(FLERR,"Variational dynamic coordinator supports Verlet only");
  auto *bond=dynamic_cast<BondBPMPeriVariational *>(force->bond);
  if(!bond)error->all(FLERR,"Variational dynamic coordinator requires its isolated bond style");
  if(material && material!=bond)
    error->all(FLERR,"Cannot replace variational material during dynamic history");
  if(!material){material=bond;material->attach_dynamic(true);}
  int nve_count=0,external_count=0,indent_count=0;bool converted=false;
  for(auto *fix : modify->get_fix_list()) {
    if(fix==this){converted=true;continue;}
    if(std::strcmp(fix->style,"nve")==0){nve_count++;continue;}
    if(std::strcmp(fix->style,"move")==0 || std::strcmp(fix->style,"property/atom")==0)continue;
    if(std::strcmp(fix->style,"external")==0 && std::strcmp(fix->id,"v3forces")==0 && converted){external_count++;continue;}
    if(std::strcmp(fix->style,"indent")==0 && std::strcmp(fix->id,"indenter")==0 && converted){indent_count++;continue;}
    error->all(FLERR,"Unadmitted fix or external force before physical SI conversion");
  }
  if(indent_count>1 || (indent_count && external_count))error->all(FLERR,"One audited indenter, without external traction, is required");
  if(external_count>1)error->all(FLERR,"V3 permits one audited external callback");
  if(nve_count!=1)error->all(FLERR,"Variational dynamics requires one free-node fix nve");
  const double mass=rho*material->physical_volume();
  for(int i=0;i<atom->nlocal;i++) {
    double actual=atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]];
    if(std::abs(actual-mass)>1e-12*mass)
      error->all(FLERR,"Variational dynamics requires physical rho times volume masses");
  }
}
void FixLPSVariationalDynamic::setup(int vflag)
{
  // Setup computed a fresh density force. It is not an accepted physical step.
  post_force(vflag);
}
void FixLPSVariationalDynamic::post_integrate()
{
  int failures=material->dynamic_commit(update->ntimestep);
  if(failures<0)error->all(FLERR,"Invalid or duplicate variational physical-step commit");
  if(failures)next_reneighbor=update->ntimestep;
}
void FixLPSVariationalDynamic::post_force(int)
{
  // After reverse communication; both pair and bond use the same volume convention.
  const double volume=material->physical_volume();
  for(int i=0;i<atom->nlocal;i++)
    for(int j=0;j<3;j++)atom->f[i][j]*=volume;
}
