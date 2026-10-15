// Read-only ABI observer compiled against frozen final headers; no physics writes.
#include <vector>
#include <string>
#include <mutex>
#include <unordered_map>
#include "lammps.h"
#include "force.h"
#include "bond.h"
#include "fix.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"
#define protected public
#include "atom.h"
#include "BPM/bond_bpm_peri.h"
#undef protected
using namespace LAMMPS_NS;
static void (*callback)(long long,int)=nullptr;
class FixAuditReadOnly : public Fix {
 public:
 FixAuditReadOnly(LAMMPS *p,int n,char **a):Fix(p,n,a){nevery=1;}
 int setmask() override {return FixConst::END_OF_STEP;}
 void setup(int) override {if(callback)callback(update->ntimestep,0);}
 void end_of_step() override {if(callback)callback(update->ntimestep,1);}
};
static Fix* creator(LAMMPS*p,int n,char**a){return new FixAuditReadOnly(p,n,a);}
extern "C" void register_observer(void (*fn)(long long,int)) {
 callback=fn;Modify::fix_styles().set_plugin("audit/readonly/finalv3",creator);
}
extern "C" double* observe(void*p,int field){
 auto*b=dynamic_cast<BondBPMPeri*>(static_cast<LAMMPS*>(p)->force->bond);
 return !b?nullptr:field==0?b->wvolume:field==1?b->theta:b->winv;
}
extern "C" long long counter(void*p,int field){
 auto*l=static_cast<LAMMPS*>(p);auto*a=l->atom;
 switch(field){case 0:return a->nextsort;case 1:return l->neighbor->ago;case 2:return a->nmax;case 3:return a->nbins;case 4:return a->sortfreq;default:return -1;}
}
extern "C" int* permutation(void*p){return static_cast<LAMMPS*>(p)->atom->permute;}
// Return all owned endpoint records. Never dereference ghost bond/history arrays.
extern "C" int endpoint_records(void*p,int*out,int capacity){
 auto*a=static_cast<LAMMPS*>(p)->atom;int k=0;
 for(int i=0;i<a->nlocal;i++){
  if(a->num_bond[i]<0 || a->num_bond[i]>a->bond_per_atom)return -1;
  for(int j=0;j<a->num_bond[i];j++){
   if(k>=capacity)return -2;
   out[3*k]=a->tag[i];out[3*k+1]=a->bond_atom[i][j];out[3*k+2]=a->bond_type[i][j];k++;
  }
 }return k;
}
#include "pair.h"
#include "neigh_list.h"
// Pure Pair::single readout of the actual neighbor list and special masks.
// This reconstructs contact forces in an audit buffer, never atom->f.
extern "C" int contact_readout(void*p,double*out,int*ends,int capacity){
 auto*l=static_cast<LAMMPS*>(p);auto*a=l->atom;auto*pair=l->force->pair;auto*list=pair->list;
 if(!list)return -1;
 int nt=a->nlocal+a->nghost;for(int i=0;i<3*nt;i++)out[i]=0.0;int count=0;
 for(int ii=0;ii<list->inum;ii++){
  int i=list->ilist[ii];
  for(int jj=0;jj<list->numneigh[i];jj++){
   int encoded=list->firstneigh[i][jj];int sb=(encoded >> SBBITS)&3;
   if(l->force->special_lj[sb]==0.0)continue;
   int j=encoded&NEIGHMASK;double dx[3],rsq=0.;for(int k=0;k<3;k++){dx[k]=a->x[i][k]-a->x[j][k];rsq+=dx[k]*dx[k];}
   if(rsq>=pair->cutsq[a->type[i]][a->type[j]])continue;
   double fp=0.;pair->single(i,j,a->type[i],a->type[j],rsq,1.,1.,fp);if(fp==0.)continue;
   if(count>=capacity)return -2;ends[2*count]=a->tag[i];ends[2*count+1]=a->tag[j];count++;
   for(int k=0;k<3;k++){out[3*i+k]+=dx[k]*fp;if(l->force->newton_pair || j<a->nlocal)out[3*j+k]-=dx[k]*fp;}
  }
 }return count;
}
