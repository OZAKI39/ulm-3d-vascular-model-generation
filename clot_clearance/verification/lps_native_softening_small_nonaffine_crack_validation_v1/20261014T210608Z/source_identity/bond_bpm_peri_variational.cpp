// Exact port of frozen 20261013T181824Z endpoint-star compliance, not shared v_b.
#include "bond_bpm_peri_variational.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "lammps.h"
#include "memory.h"
#include "update.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <chrono>
#include <cstring>
using namespace LAMMPS_NS;
BondBPMPeriVariational::BondBPMPeriVariational(LAMMPS *p):Bond(p){partial_flag=1;writedata=0;reinitflag=0;}
BondBPMPeriVariational::~BondBPMPeriVariational(){if(allocated)memory->destroy(setflag);}
void BondBPMPeriVariational::allocate(){if(!allocated){memory->create(setflag,atom->nbondtypes+1,"variational:setflag");allocated=1;std::fill(setflag,setflag+atom->nbondtypes+1,0);}}
void BondBPMPeriVariational::settings(int argc,char **argv){
 if(argc!=1 || ready)error->all(FLERR,"variational requires one immutable reference file; reinitialization forbidden");
 std::ifstream f(argv[0]);std::string magic;f>>magic>>n>>nb>>V;
 if(magic!="LPSVAR_REFERENCE_V1" || n<1 || n>4096 || nb<1 || nb>1000000 || !(V>0))error->all(FLERR,"Invalid variational reference header");
 tags.resize(n);m.resize(n);bid.resize(nb);ei.resize(nb);ej.resize(nb);types.resize(nb);r.resize(nb);h.resize(nb);
 for(int a=0;a<n;a++){f>>tags[a]>>m[a];if(!ordinal.emplace(tags[a],a).second)error->all(FLERR,"Duplicate reference atom ID");}
 for(int b=0;b<nb;b++){tagint i,j;int im[3];f>>bid[b]>>i>>j>>types[b]>>r[b]>>h[b]>>im[0]>>im[1]>>im[2];
  if(!ordinal.count(i)||!ordinal.count(j)||i==j||im[0]||im[1]||im[2])error->all(FLERR,"Invalid or unsupported periodic reference identity");ei[b]=ordinal.at(i);ej[b]=ordinal.at(j);
 }
 if(!f)error->all(FLERR,"Incomplete variational reference file");kap.assign(nb,0);d.assign(nb,0);failed.assign(nb,0);ready=true;
}
void BondBPMPeriVariational::coeff(int argc,char **argv){
 if(argc!=6 || std::string(argv[0])!="1" || atom->nbondtypes!=1)error->all(FLERR,"variational supports one frozen type: 1 K G horizon s0 sc");
 K=utils::numeric(FLERR,argv[1],false,lmp);G=utils::numeric(FLERR,argv[2],false,lmp);delta=utils::numeric(FLERR,argv[3],false,lmp);s0=utils::numeric(FLERR,argv[4],false,lmp);sc=utils::numeric(FLERR,argv[5],false,lmp);
 allocate();setflag[1]=1;validate();
}
void BondBPMPeriVariational::validate(){
 if(!ready || !(G>0 && K>5*G/3 && sc>s0 && s0>0 && delta>0))error->all(FLERR,"Invalid variational positive-energy parameters or absent history");
 std::set<long long> ids;edgekey.clear();ordinal.clear();for(int a=0;a<n;a++)if(tags[a]<=0||!(m[a]>0)||!ordinal.emplace(tags[a],a).second)error->all(FLERR,"Invalid reference node");
 std::vector<double> mr(n,0);
 for(int b=0;b<nb;b++){
  if(bid[b]<=0||!ids.insert(bid[b]).second||ei[b]<0||ej[b]<0||ei[b]>=n||ej[b]>=n||types[b]!=1||!(r[b]>0 && h[b]>0)||!std::isfinite(kap[b])||kap[b]<0||!std::isfinite(d[b])||d[b]<0||d[b]>1||failed[b]!=(d[b]==1))error->all(FLERR,"Invalid original bond/history record");
  auto key=std::minmax(tags[ei[b]],tags[ej[b]]);if(!edgekey.emplace(key,b).second)error->all(FLERR,"Duplicate endpoint pairs not supported; refusing identity collision");
  mr[ei[b]]+=h[b]*r[b]*V;mr[ej[b]]+=h[b]*r[b]*V;
 }
 for(int a=0;a<n;a++)if(std::abs(mr[a]-m[a])>2e-12*m[a])error->all(FLERR,"Reference weighted volume inconsistent with full original graph");
 A=V*(K-5*G/3);x.resize(3*n);e.resize(nb);dir.resize(3*nb);L.resize(2*nb);D.resize(2*nb);k.resize(2*nb);v.resize(2*nb);tdir.resize(2*nb);theta.resize(n);fn.resize(3*n);
 for(int q=0;q<2*nb;q++){int b=q%nb,a=q<nb?ei[b]:ej[b];L[q]=3*V*h[b]/m[a];D[q]=15*G*V*V*h[b]/(r[b]*m[a]);k[q]=9*K*V*V*h[b]/(r[b]*m[a]);}
}
void BondBPMPeriVariational::init_style(){
 if(domain->xperiodic||domain->yperiodic||domain->zperiodic||domain->triclinic)error->all(FLERR,"variational research scope requires nonperiodic orthogonal box");
 if(force->newton_bond)error->all(FLERR,"variational requires newton bond off");
 if(atom->natoms!=n || std::string(update->unit_style)!="si")error->all(FLERR,"variational requires exact reference atoms and SI units");
 if(force->special_lj[1]!=0 || force->special_lj[2]!=1 || force->special_lj[3]!=1)error->all(FLERR,"variational requires special_bonds lj 0 1 1");
 validate();
}
void BondBPMPeriVariational::gather(){
 std::fill(x.begin(),x.end(),0);std::vector<int> count(n,0);
 for(int a=0;a<atom->nlocal;a++){auto it=ordinal.find(atom->tag[a]);if(it==ordinal.end())error->all(FLERR,"Unknown owned atom tag");int i=it->second;count[i]++;for(int c=0;c<3;c++)x[3*i+c]=atom->x[a][c];}
 MPI_Allreduce(MPI_IN_PLACE,x.data(),3*n,MPI_DOUBLE,MPI_SUM,world);MPI_Allreduce(MPI_IN_PLACE,count.data(),n,MPI_INT,MPI_SUM,world);
 for(int i=0;i<n;i++)if(count[i]!=1)error->all(FLERR,"Owned atom identity not unique across ranks");
 for(int b=0;b<nb;b++){double rr=0;for(int c=0;c<3;c++){double y=x[3*ej[b]+c]-x[3*ei[b]+c];dir[3*b+c]=y;rr+=y*y;}rr=std::sqrt(rr);if(!(rr>0))error->all(FLERR,"Coincident original interaction endpoints");e[b]=rr-r[b];for(int c=0;c<3;c++)dir[3*b+c]/=rr;}
}
void BondBPMPeriVariational::evaluate(){
 auto start=std::chrono::steady_clock::now();gather();
 if(std::all_of(failed.begin(),failed.end(),[](int q){return q==1;})) {
  U=0;std::fill(fn.begin(),fn.end(),0);std::fill(theta.begin(),theta.end(),0);
  std::fill(v.begin(),v.end(),0);std::fill(tdir.begin(),tdir.end(),0);
  zero_graph_skips++;cache_x=x;cache_valid=true;
  seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return;
 }
 std::vector<double> numerator(n,0),denominator(n,1/A),w(2*nb),hh(2*nb),den(2*nb);
 for(int q=0;q<2*nb;q++){int b=q%nb,a=q<nb?ei[b]:ej[b];den[q]=D[q]*d[b]+k[q]*(1-d[b]);w[q]=k[q]*(1-d[b])/den[q];hh[q]=d[b]/den[q];numerator[a]+=L[q]*w[q]*e[b];denominator[a]+=L[q]*L[q]*hh[q];}
 for(int a=0;a<n;a++)numerator[a]/=denominator[a];
 U=0;std::fill(theta.begin(),theta.end(),0);std::fill(fn.begin(),fn.end(),0);
 for(int q=0;q<2*nb;q++){int b=q%nb,a=q<nb?ei[b]:ej[b];double coupling=L[q]*numerator[a],g=D[q]*e[b]+coupling;v[q]=w[q]*e[b]-hh[q]*coupling;tdir[q]=w[q]*g;theta[a]+=L[q]*v[q];double z=g/den[q];U+=.5*D[q]*v[q]*v[q]+.5*k[q]*(1-d[b])*d[b]*z*z;}
 for(double t:theta)U+=.5*A*t*t;
 for(int b=0;b<nb;b++){double t=tdir[b]+tdir[b+nb];for(int c=0;c<3;c++){double f=t*dir[3*b+c];fn[3*ei[b]+c]+=f;fn[3*ej[b]+c]-=f;}}
 cache_x=x;cache_valid=true;evaluations++;seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
void BondBPMPeriVariational::compute(int ef,int vf){
 ev_init(ef,vf);
 if(cache_valid){gather();if(x==cache_x)cache_hits++;else evaluate();}
 else evaluate();
 // Native BPM uses force density. Exported fn,U retain physical N,J.
 for(int a=0;a<atom->nlocal;a++){int i=ordinal.at(atom->tag[a]);for(int c=0;c<3;c++)atom->f[a][c]+=fn[3*i+c]/V;}
 if(eflag_global)energy=comm->me==0?U/V:0;
 if(eflag_atom)error->all(FLERR,"Per-atom energy not exposed; use complete variational global energy");
 if(vflag_global && comm->me==0)for(int b=0;b<nb;b++){double t=(tdir[b]+tdir[b+nb])*(r[b]+e[b]);double *nvec=dir.data()+3*b;virial[0]-=t*nvec[0]*nvec[0];virial[1]-=t*nvec[1]*nvec[1];virial[2]-=t*nvec[2]*nvec[2];virial[3]-=t*nvec[0]*nvec[1];virial[4]-=t*nvec[0]*nvec[2];virial[5]-=t*nvec[1]*nvec[2];}
}
void BondBPMPeriVariational::sync_topology(){
 // Only OWNED arrays are accessed. No history is indexed by these mutable slots.
 for(int a=0;a<atom->nlocal;a++){
  for(int b=0;b<atom->num_bond[a];){auto key=std::minmax(atom->tag[a],atom->bond_atom[a][b]);auto it=edgekey.find(key);if(it==edgekey.end())error->all(FLERR,"Native topology contains unknown original bond");
   if(failed[it->second]){int last=--atom->num_bond[a];atom->bond_atom[a][b]=atom->bond_atom[a][last];atom->bond_type[a][b]=atom->bond_type[a][last];}else b++;
  }
  // All higher-order weights are 1, so only direct exclusions need removal.
  for(int s=0;s<atom->nspecial[a][0];){auto key=std::minmax(atom->tag[a],atom->special[a][s]);auto it=edgekey.find(key);
   if(it!=edgekey.end()&&failed[it->second]){for(int q=s;q<atom->nspecial[a][2]-1;q++)atom->special[a][q]=atom->special[a][q+1];for(int c=0;c<3;c++)atom->nspecial[a][c]--;}else s++;
  }
 }
 atom->nbonds=nb-std::count(failed.begin(),failed.end(),1);
}
int BondBPMPeriVariational::action(int op,long long tok,const double *kk,const double *dd){
 if(update->whichflag!=0)return -1; // Only between accepted run segments, never inside force evaluation.
 if(op==1){
  if(dynamic_attached)return -7;
  if(imported||commits)return -2;
  if(comm->me==0){std::copy(kk,kk+nb,kap.begin());std::copy(dd,dd+nb,d.begin());}
  MPI_Bcast(kap.data(),nb,MPI_DOUBLE,0,world);MPI_Bcast(d.data(),nb,MPI_DOUBLE,0,world);
  for(int b=0;b<nb;b++)failed[b]=(d[b]==1);validate();imported=true;cache_valid=false;sync_topology();return 0;
 }
 if(op!=2)return -3;
 if(dynamic_attached)return -7;
 return commit_state(tok);
}
void BondBPMPeriVariational::attach_dynamic(bool attached){
 if(attached && dynamic_attached)error->all(FLERR,"Only one variational dynamic coordinator is allowed");
 dynamic_attached=attached;
}
int BondBPMPeriVariational::dynamic_commit(long long step){
 if(!dynamic_attached || update->whichflag!=1 || update->setupflag)return -8;
 if(commits && step!=token+1)return -9;
 int result=commit_state(step);
 return result==0?int(last_new_failures):-10;
}
int BondBPMPeriVariational::commit_state(long long tok){
 long long low=tok,high=tok;MPI_Allreduce(MPI_IN_PLACE,&low,1,MPI_LONG_LONG,MPI_MIN,world);MPI_Allreduce(MPI_IN_PLACE,&high,1,MPI_LONG_LONG,MPI_MAX,world);if(low!=high || tok<token)return -4;
 evaluate();if(tok==token)return x==accepted_x?1:-5;
 double oldU=U;std::vector<double> nextk(nb),nextd(nb);
 for(int b=0;b<nb;b++){nextk[b]=std::max(kap[b],std::max(e[b]/r[b],0.));double z=nextk[b];nextd[b]=z<=s0?0:(z>=sc?1:sc*(z-s0)/(z*(sc-s0)));if(nextd[b]<d[b]-2e-13)return -6;}
 bool changed=nextd!=d;
 kap.swap(nextk);d.swap(nextd);
 if(changed)evaluate();
 lastdrop=oldU-U;Diss+=lastdrop;last_new_failures=0;
 for(int b=0;b<nb;b++)if(d[b]==1&&!failed[b]){failed[b]=1;events++;last_new_failures++;}
 if(last_new_failures)sync_topology();token=tok;commits++;accepted_x=x;return 0;
}
int BondBPMPeriVariational::output(int field,double *out){
 const std::vector<double>* p=nullptr;
 switch(field){case 0:p=&fn;break;case 1:p=&x;break;case 2:p=&kap;break;case 3:p=&d;break;case 4:p=&m;break;case 5:p=&v;break;case 6:p=&tdir;break;case 7:p=&theta;break;
 case 8:{double s[]={U,Diss,lastdrop,double(events),double(commits),double(token),double(evaluations),seconds,double(n),double(nb)};std::copy(s,s+10,out);return 10;}
 case 11:{double s[]={double(cache_hits),double(zero_graph_skips),double(last_new_failures),V,double(dynamic_attached)};std::copy(s,s+5,out);return 5;}
 case 9:for(int b=0;b<nb;b++)out[b]=double(bid[b]);return nb;
 case 10:for(int b=0;b<nb;b++)out[b]=failed[b];return nb;
 default:return -1;}
 std::copy(p->begin(),p->end(),out);return p->size();
}
double BondBPMPeriVariational::equilibrium_distance(int){return delta;}
double BondBPMPeriVariational::single(int,double,int,int,double &ff){ff=0;error->all(FLERR,"variational is state-coupled; single-bond energy query unsupported");return 0;}
std::string BondBPMPeriVariational::serialize()const{
 std::ostringstream s;s<<std::setprecision(17)<<"LPSVAR_RESTART_V1\n"<<n<<' '<<nb<<' '<<V<<' '<<K<<' '<<G<<' '<<delta<<' '<<s0<<' '<<sc<<'\n'<<token<<' '<<commits<<' '<<events<<' '<<Diss<<' '<<lastdrop<<' '<<imported<<'\n';
 for(int a=0;a<n;a++)s<<tags[a]<<' '<<m[a]<<'\n';
 for(int b=0;b<nb;b++)s<<bid[b]<<' '<<ei[b]<<' '<<ej[b]<<' '<<types[b]<<' '<<r[b]<<' '<<h[b]<<' '<<kap[b]<<' '<<d[b]<<' '<<failed[b]<<'\n';
 s<<accepted_x.size()<<'\n';for(double z:accepted_x)s<<z<<' ';s<<'\n';return s.str();
}
void BondBPMPeriVariational::deserialize(const std::string &text){
 std::istringstream s(text);std::string magic;s>>magic>>n>>nb>>V>>K>>G>>delta>>s0>>sc;
 if(magic!="LPSVAR_RESTART_V1"||n<1||n>4096||nb<1||nb>1000000)error->all(FLERR,"Missing or incompatible variational restart history");
 s>>token>>commits>>events>>Diss>>lastdrop>>imported;tags.resize(n);m.resize(n);bid.resize(nb);ei.resize(nb);ej.resize(nb);types.resize(nb);r.resize(nb);h.resize(nb);kap.resize(nb);d.resize(nb);failed.resize(nb);
 for(int a=0;a<n;a++)s>>tags[a]>>m[a];for(int b=0;b<nb;b++)s>>bid[b]>>ei[b]>>ej[b]>>types[b]>>r[b]>>h[b]>>kap[b]>>d[b]>>failed[b];
 size_t count;s>>count;if(count!=0&&count!=size_t(3*n))error->all(FLERR,"Incomplete variational accepted-state coordinates");accepted_x.resize(count);for(double &z:accepted_x)s>>z;
 if(!s)error->all(FLERR,"Truncated variational restart");ready=true;cache_valid=false;allocate();setflag[1]=1;validate();
}
void BondBPMPeriVariational::write_restart(FILE *fp){std::string s=serialize();int size=s.size();fwrite(&size,sizeof(int),1,fp);fwrite(s.data(),1,size,fp);}
void BondBPMPeriVariational::read_restart(FILE *fp){int size=0;if(comm->me==0)utils::sfread(FLERR,&size,sizeof(int),1,fp,nullptr,error);MPI_Bcast(&size,1,MPI_INT,0,world);if(size<32||size>200000000)error->all(FLERR,"Invalid variational restart payload");std::string s(size,'\0');if(comm->me==0)utils::sfread(FLERR,&s[0],1,size,fp,nullptr,error);MPI_Bcast(&s[0],size,MPI_CHAR,0,world);deserialize(s);}
extern "C" {
__attribute__((visibility("default"))) int lpsvar_action(void *p,int op,long long tok,const double *k,const double *d){auto *b=dynamic_cast<BondBPMPeriVariational*>(static_cast<LAMMPS*>(p)->force->bond);return b?b->action(op,tok,k,d):-99;}
__attribute__((visibility("default"))) int lpsvar_output(void *p,int field,double *out){auto *b=dynamic_cast<BondBPMPeriVariational*>(static_cast<LAMMPS*>(p)->force->bond);return b?b->output(field,out):-99;}
}
