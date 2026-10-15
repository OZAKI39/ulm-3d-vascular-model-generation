// Isolated research extension; LAMMPS GPL-compatible.
#ifdef BOND_CLASS
BondStyle(bpm/peri/variational,BondBPMPeriVariational);
#else
#ifndef LMP_BOND_BPM_PERI_VARIATIONAL_H
#define LMP_BOND_BPM_PERI_VARIATIONAL_H
#include "bond.h"
#include <vector>
#include <map>
#include <string>
namespace LAMMPS_NS {
class BondBPMPeriVariational : public Bond {
 public:
  explicit BondBPMPeriVariational(class LAMMPS *);
  ~BondBPMPeriVariational() override;
  void settings(int,char **) override;
  void coeff(int,char **) override;
  void init_style() override;
  void compute(int,int) override;
  double equilibrium_distance(int) override;
  double single(int,double,int,int,double &) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  int action(int,long long,const double *,const double *);
  int output(int,double *);
 private:
  int n=0,nb=0;
  bool ready=false,imported=false;
  double V=0,K=0,G=0,delta=0,s0=0,sc=0,A=0,U=0,Diss=0,lastdrop=0,seconds=0;
  long long token=-1,commits=0,events=0,evaluations=0;
  std::vector<long long> bid;
  std::vector<tagint> tags;
  std::map<tagint,int> ordinal;
  std::map<std::pair<tagint,tagint>,int> edgekey;
  std::vector<int> ei,ej,types,failed;
  std::vector<double> m,r,h,kap,d,x,accepted_x,e,dir,L,D,k,v,tdir,theta,fn;
  void allocate();
  void validate();
  void gather();
  void evaluate();
  void sync_topology();
  std::string serialize() const;
  void deserialize(const std::string &);
};
}
#endif
#endif
