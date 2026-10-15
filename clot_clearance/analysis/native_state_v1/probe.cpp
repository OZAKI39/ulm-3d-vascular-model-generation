// Read-only internal-state observer against the exact production ABI.
// No LAMMPS sources or forces are modified. grow() is an explicit test stimulus.
#include "lammps.h"
#include "force.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#define protected public
#include "BPM/bond_bpm_peri.h"
#undef protected
extern "C" double *observe(void *p,int field) {
 auto *b=dynamic_cast<LAMMPS_NS::BondBPMPeri*>(static_cast<LAMMPS_NS::LAMMPS*>(p)->force->bond);
 if (!b) return nullptr;
 return field==0 ? b->wvolume : (field==1 ? b->theta : b->winv);
}
extern "C" void grow_atom_storage(void *p) {
 auto *a=static_cast<LAMMPS_NS::LAMMPS*>(p)->atom;a->avec->grow(a->nmax+100);
}
