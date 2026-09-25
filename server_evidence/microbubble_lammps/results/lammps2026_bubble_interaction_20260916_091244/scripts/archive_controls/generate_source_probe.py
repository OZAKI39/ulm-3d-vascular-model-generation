from pathlib import Path
import json,hashlib,re,subprocess,sys
R=Path('/home/lzy/projects/compre_output/lammps2026_bubble_interaction/20260916_091244');P=R/'provenance/upstream_lubrication';O=R/'source_correction_audit';O.mkdir(exist_ok=True)
source=P/'pair_lubricate_poly.cpp';assert source.read_bytes()==(P/'TARGET_RAW_pair_lubricate_poly.cpp').read_bytes()
extracts={}
for name,filename,upoly in [('target_poly','pair_lubricate_poly.cpp',False),('target_upoly','pair_lubricateU_poly.cpp',True),('old_poly','OLD_RAW_pair_lubricate_poly.cpp',False)]:
 text=(P/filename).read_text();marker='        if (flaglog) {';start=text.index(marker)
 # The upstream block runs to the velocity/force section immediately after its else branch.
 end=text.index('        // Relative',start) if upoly else text.index('        // relative',start)
 if upoly:end=min(end,text.index('        // Find force',start)) if '        // Find force' in text[start:] else end
 block=text[start:end]
 # Exact source snippet; cut at the end of the complete flaglog/else coefficient expression.
 last=block.index('a_sq = pre[0]*',block.index('} else')) if upoly else block.index('} else a_sq')
 stop=block.index(';',last)+1
 if upoly:stop=block.index('}',stop)+1
 block=block[:stop]
 (O/(name+'_literal_source.txt')).write_text(block+'\n');extracts[name]={'source_file':filename,'source_sha256':hashlib.sha256((P/filename).read_bytes()).hexdigest(),'first_line':text[:start].count('\n')+1,'literal_sha256':hashlib.sha256(block.encode()).hexdigest(),'bytes':block,'upoly':upoly}
head='''#include <cmath>
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>
using std::log;using std::pow;
constexpr double MY_PI=3.141592653589793238462643383279502884;
'''
for name,e in extracts.items():
 head+='std::array<double,2> '+name+'(double radi,double radj,double gap,double mu,int flaglog){\n double a_sq=0,a_sh=0,a_pu=0,h_sep=gap/radi;\n double beta0=radj/radi,beta1=1+beta0;\n double beta[2][5]={{0},{0}},pre[2];beta[0][1]=beta0;beta[1][1]=beta1;pre[1]=8.0*(pre[0]=MY_PI*mu*radi)*radi*radi;pre[0]*=6.0;\n'+e['bytes']+'\n return {a_sq,a_sh};}\n'
head+='''int main(int argc,char**argv){std::ifstream in(argv[1]);std::ofstream out(argv[2]);out<<std::setprecision(17);out<<"id,flaglog,target_poly_sq,target_poly_sq_swap,target_poly_sh,target_poly_sh_swap,target_upoly_sq,target_upoly_sq_swap,target_upoly_sh,target_upoly_sh_swap,old_poly_sq,old_poly_sq_swap,old_poly_sh,old_poly_sh_swap\\n";long id;double a,b,h,mu;
while(in>>id>>a>>b>>h>>mu)for(int f=0;f<=1;++f){auto x=target_poly(a,b,h,mu,f),y=target_poly(b,a,h,mu,f),u=target_upoly(a,b,h,mu,f),v=target_upoly(b,a,h,mu,f),o=old_poly(a,b,h,mu,f),p=old_poly(b,a,h,mu,f);out<<id<<","<<f<<","<<x[0]<<","<<y[0]<<","<<x[1]<<","<<y[1]<<","<<u[0]<<","<<v[0]<<","<<u[1]<<","<<v[1]<<","<<o[0]<<","<<p[0]<<","<<o[1]<<","<<p[1]<<"\\n";}}
'''
(O/'source_coefficient_probe.cpp').write_text(head)
(O/'EXTRACTION_PROVENANCE.json').write_text(json.dumps({k:{x:y for x,y in v.items() if x!='bytes'} for k,v in extracts.items()},indent=2)+'\n')
(O/'AUDIT_CONTRACT.json').write_text(json.dumps({'role':'PHASE_A_UPSTREAM_SOURCE_CORRECTION_AUDIT_ONLY','pairs':10000,'seed':20260916,'relative_swap_symmetry_gate':1e-12,'radii_m':[.375e-6,2.625e-6],'gap_over_radius_sum':[.001,.1],'mu_Pa_s':.001,'flags':[0,1],'phase_C_custom_implementation':'NOT_STARTED','math_snippets':'verbatim from official source, not rewritten formulas'},indent=2)+'\n')
