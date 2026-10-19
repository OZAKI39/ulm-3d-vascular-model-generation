"""Build V4 in isolation, with nonuniform endpoint-star quadrature and SI forces."""
from common import *
import subprocess,shutil,difflib,time
def replace(text,old,new,count=1):
 assert text.count(old)==count,(old,text.count(old));return text.replace(old,new)
def main():
 start=time.perf_counter();build=O/'build';build.mkdir(exist_ok=True);assert not (build/'liblammps.so').exists()
 frozen=O/'provenance/v3_source';cpp=(frozen/'bond_bpm_peri_variational.cpp').read_text();header=(frozen/'bond_bpm_peri_variational.h').read_text();fix=(frozen/'fix_lps_variational_dynamic.cpp').read_text()
 cpp=replace(cpp,'magic!="LPSVAR_REFERENCE_V1"','(magic!="LPSVAR_REFERENCE_V1" && magic!="LPSVAR_REFERENCE_V2")')
 cpp=replace(cpp,'tags.resize(n);m.resize(n);bid.resize(nb);','tags.resize(n);m.resize(n);volumes.assign(n,V);Ai.resize(n);bid.resize(nb);',2)
 cpp=replace(cpp,'f>>tags[a]>>m[a];','f>>tags[a]>>m[a];if(magic=="LPSVAR_REFERENCE_V2")f>>volumes[a];')
 cpp=replace(cpp,'!(m[a]>0)','!(m[a]>0 && volumes[a]>0 && std::isfinite(volumes[a]))')
 cpp=replace(cpp,'mr[ei[b]]+=h[b]*r[b]*V;mr[ej[b]]+=h[b]*r[b]*V;','mr[ei[b]]+=h[b]*r[b]*volumes[ej[b]];mr[ej[b]]+=h[b]*r[b]*volumes[ei[b]];')
 cpp=replace(cpp,'A=V*(K-5*G/3);','A=V*(K-5*G/3);for(int a=0;a<n;a++)Ai[a]=volumes[a]*(K-5*G/3);')
 cpp=replace(cpp,'L[q]=3*V*h[b]/m[a];D[q]=15*G*V*V*h[b]/(r[b]*m[a]);k[q]=9*K*V*V*h[b]/(r[b]*m[a]);','int other=q<nb?ej[b]:ei[b];double vv=volumes[a]*volumes[other];L[q]=3*volumes[other]*h[b]/m[a];D[q]=15*G*vv*h[b]/(r[b]*m[a]);k[q]=9*K*vv*h[b]/(r[b]*m[a]);')
 cpp=replace(cpp,'std::vector<double> numerator(n,0),denominator(n,1/A),w(2*nb),hh(2*nb),den(2*nb);','std::vector<double> numerator(n,0),denominator(n),w(2*nb),hh(2*nb),den(2*nb);for(int a=0;a<n;a++)denominator[a]=1/Ai[a];')
 cpp=replace(cpp,'for(double t:theta)U+=.5*A*t*t;','for(int a=0;a<n;a++)U+=.5*Ai[a]*theta[a]*theta[a];')
 cpp=replace(cpp,'fn[3*i+c]/V','fn[3*i+c]')
 cpp=replace(cpp,'?U/V:0','?U:0')
 cpp=cpp.replace('// Native BPM uses force density. Exported fn,U retain physical N,J.','// V4 uses physical N,J throughout; no global-volume conversion.')
 cpp=replace(cpp,'case 14:std::copy(complete_by_bond.begin(),complete_by_bond.end(),out);return nb;','case 14:std::copy(complete_by_bond.begin(),complete_by_bond.end(),out);return nb;\n case 15:std::copy(volumes.begin(),volumes.end(),out);return n;')
 cpp=replace(cpp,'<<"LPSVAR_ID_SOFTENING_RESTART_V4\\n"','<<"LPSVAR_NONUNIFORM_RESTART_V5\\n"')
 cpp=replace(cpp,"s<<tags[a]<<' '<<m[a]<<'\\n';","s<<tags[a]<<' '<<m[a]<<' '<<volumes[a]<<'\\n';")
 cpp=replace(cpp,'magic!="LPSVAR_ID_SOFTENING_RESTART_V4")','magic!="LPSVAR_ID_SOFTENING_RESTART_V4" && magic!="LPSVAR_NONUNIFORM_RESTART_V5")')
 cpp=replace(cpp,'magic=="LPSVAR_ID_SOFTENING_RESTART_V4")','magic=="LPSVAR_ID_SOFTENING_RESTART_V4" || magic=="LPSVAR_NONUNIFORM_RESTART_V5")',2)
 cpp=replace(cpp,'for(int a=0;a<n;a++)s>>tags[a]>>m[a];','for(int a=0;a<n;a++){s>>tags[a]>>m[a];if(magic=="LPSVAR_NONUNIFORM_RESTART_V5")s>>volumes[a];}')
 header=replace(header,'double physical_volume() const { return V; }','double physical_volume() const { return V; }\n  double nodal_volume(tagint id) const { return volumes[ordinal.at(id)]; }')
 header=replace(header,'std::vector<double> m,r,h,kap,d,x,accepted_x,e,dir,L,D,k,v,tdir,theta,fn;','std::vector<double> volumes,Ai;\n  std::vector<double> m,r,h,kap,d,x,accepted_x,e,dir,L,D,k,v,tdir,theta,fn;')
 fix=replace(fix,'  const double mass=rho*material->physical_volume();','')
 fix=replace(fix,'    double actual=atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]];','    const double mass=rho*material->nodal_volume(atom->tag[i]);\n    double actual=atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]];')
 fix=replace(fix,'  const double volume=material->physical_volume();\n  for(int i=0;i<atom->nlocal;i++)\n    for(int j=0;j<3;j++)atom->f[i][j]*=volume;','  // Both the isolated V4 material and contact kernels already return N.\n  // External forces follow this coordinator unchanged; never multiply by V twice.')
 source=R/'verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z/source/src/BPM'
 pair=(source/'pair_bpm_peri.cpp').read_text()
 pair=replace(pair,'kspring[itype][jtype] * vfrac_eff * (dr / cut[itype][jtype])','kspring[itype][jtype] * vfrac_eff * vfrac_eff * (dr / cut[itype][jtype])',2)
 pair=replace(pair,'fpair * vfrac_eff, delx, dely, delz','fpair, delx, dely, delz')
 pair=pair.replace('// The mean nodal volume keeps the contact Newton-third-law balanced; for','// V4 physical symmetric mean-volume-squared contact quadrature; for')
 generated={'bond_bpm_peri_variational.cpp':cpp,'bond_bpm_peri_variational.h':header,'fix_lps_variational_dynamic.cpp':fix,'fix_lps_variational_dynamic.h':(frozen/'fix_lps_variational_dynamic.h').read_text(),'pair_bpm_peri.cpp':pair,'pair_bpm_peri.h':(source/'pair_bpm_peri.h').read_text()}
 for name,text in generated.items():(A/name).write_text(text)
 audit=json.loads((V2_O/'BUILD_AUDIT.json').read_text());objects=list(audit['object_sha256'])
 for path,want in audit['object_sha256'].items():assert sha(path)==want,path
 commands=[]
 with (build/'build.log').open('w') as log:
  for old in audit['commands'][:-1]:
   if '-c' not in old:continue
   src=Path(old[old.index('-c')+1]);name=src.name
   if name not in ['bond_bpm_peri_variational.cpp','fix_lps_variational_dynamic.cpp','force.cpp','style_bond.cpp']:continue
   newsource=A/name if name in generated else src
   flags=old[:old.index('-c')];flags[1:1]=['-I'+str(A)];obj=build/(name+'.o');cmd=flags+['-c',str(newsource),'-o',str(obj)]
   subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True);commands.append(cmd)
   matches=[i for i,p in enumerate(objects) if Path(p).name==name+'.o'];assert len(matches)==1;objects[matches[0]]=str(obj)
  flags=audit['commands'][0][:audit['commands'][0].index('-c')];flags[1:1]=['-I'+str(A)];obj=build/'pair_bpm_peri.cpp.o';cmd=flags+['-c',str(A/'pair_bpm_peri.cpp'),'-o',str(obj)]
  subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True);commands.append(cmd)
  matches=[i for i,p in enumerate(objects) if Path(p).name=='pair_bpm_peri.cpp.o'];assert len(matches)==1;objects[matches[0]]=str(obj)
  rsp=build/'objects.rsp';rsp.write_text('\n'.join(objects)+'\n');link=audit['commands'][-1].copy();link[link.index('-o')+1]=str(build/'liblammps.so.0');link=[('@'+str(rsp)) if p.startswith('@') else p for p in link]
  subprocess.run(link,cwd=audit['link_cwd'],stdout=log,stderr=subprocess.STDOUT,check=True);commands.append(link)
 (build/'liblammps.so').symlink_to('liblammps.so.0')
 dump(O/'BUILD_AUDIT.json',dict(status='BUILT_NOT_YET_ADMITTED',commands=commands,link_cwd=audit['link_cwd'],object_sha256={p:sha(p) for p in objects},source_sha256={name:sha(A/name) for name in generated},library_sha256=sha(LIB),baseline_library_unchanged=sha(V2_O/'build/liblammps.so')==audit['library_sha256'],wall_s=time.perf_counter()-start,nonuniform_contact='Same native linear repulsion, DX onset and live-bond exclusion, symmetric arithmetic-mean nodal volume squared. Force and potential in N,J; regular-volume limit recovers V3.',restart='V5 includes every original-ID nodal volume and all prior history fields'))
 print('V4_NATIVE_BUILT',sha(LIB),flush=True)
if __name__=='__main__':main()
