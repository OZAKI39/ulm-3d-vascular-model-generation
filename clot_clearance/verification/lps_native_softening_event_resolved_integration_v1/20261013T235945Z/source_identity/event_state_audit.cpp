// Isolated read-only rollback audit. No model, integration or topology mutation.
#include "lammps.h"
#include "atom.h"
extern "C" int lps_event_special(void *ptr, double *out)
{
  auto *lmp = static_cast<LAMMPS_NS::LAMMPS *>(ptr);
  auto *a = lmp->atom;
  int count = 0;
  for (int i = 0; i < a->nlocal; ++i) {
    if (out) {
      out[count] = static_cast<double>(a->tag[i]);
      for (int j = 0; j < 3; ++j) out[count+1+j] = a->nspecial[i][j];
    }
    count += 4;
    for (int j = 0; j < a->nspecial[i][2]; ++j) {
      if (out) out[count] = static_cast<double>(a->special[i][j]);
      ++count;
    }
  }
  return count;
}
