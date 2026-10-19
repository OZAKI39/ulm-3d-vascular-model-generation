"""Freeze V3 inputs and build an isolated, timing-only native library."""
from pathlib import Path
import json,hashlib,shutil,subprocess,datetime,sys,platform,importlib.metadata
R=Path('/workspace/y_clot/project')
TASK='clot_irregular_lps_discretization_v4';RUN='20261019T221000Z'
O=R/'verification'/TASK/RUN;A=R/'analysis'/TASK/'code'
V3=R/'verification/clot_fragment_surface_coupling_v3/20261018T170400Z'
V2=R/'verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,z):p.write_text(json.dumps(z,indent=2)+'\n')
def main():
    freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
    old=json.loads((V3/'SOURCE_FREEZE.json').read_text())
    for path,want in old['implementation_sha256'].items():assert sha(R/path)==want,path
    freeze['v3_implementation_sha256']=old['implementation_sha256']
    freeze['v3_input_sha256']=old['input_sha256']
    freeze['server_git_head']=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip()
    freeze['server_git_branch']=subprocess.check_output(['git','-c',f'safe.directory={R}','branch','--show-current'],cwd=R,text=True).strip()
    freeze['environment']={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'native_sha256':sha(V2/'build/liblammps.so')}
    dump(O/'SOURCE_FREEZE.json',freeze)
    shutil.copy2('/workspace/y_clot/env.sh',O/'provenance/environment.sh')
    (O/'provenance/package_versions.txt').write_text('\n'.join(sorted(f'{d.metadata["Name"]}=={d.version}' for d in importlib.metadata.distributions()))+'\n')
    p0=O/'p0';(p0/'provenance').mkdir(parents=True,exist_ok=True)
    for name in sorted(set(old['input_sha256'])|{'SOURCE_FREEZE.json'}):
        src=V3/name;dst=p0/name
        if src.exists():
            if dst.exists():assert sha(dst)==sha(src)
            else:shutil.copy2(src,dst)
    build=p0/'build';build.mkdir(exist_ok=True)
    source=(O/'provenance/v3_source/bond_bpm_peri_variational.cpp').read_text()
    timer=r'''
#include <ctime>
namespace {
struct V4Profile { double wall[10]={},cpu[10]={},calls[10]={}; } v4_profile;
struct V4Timer {
 int slot;double child_wall=0,child_cpu=0;V4Timer *parent;
 std::chrono::steady_clock::time_point start;std::clock_t cp;
 static thread_local V4Timer *top;
 V4Timer(int s):slot(s),parent(top),start(std::chrono::steady_clock::now()),cp(std::clock()){top=this;}
 ~V4Timer(){double w=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();double c=double(std::clock()-cp)/CLOCKS_PER_SEC;
  v4_profile.wall[slot]+=w-child_wall;v4_profile.cpu[slot]+=c-child_cpu;v4_profile.calls[slot]++;
  if(parent){parent->child_wall+=w;parent->child_cpu+=c;}top=parent;}
};
thread_local V4Timer *V4Timer::top=nullptr;
}
extern "C" __attribute__((visibility("default"))) int lpsvar_profile(double *out,int reset){
 for(int k=0;k<10;k++){out[3*k]=v4_profile.wall[k];out[3*k+1]=v4_profile.cpu[k];out[3*k+2]=v4_profile.calls[k];}
 if(reset)v4_profile=V4Profile{};return 30;
}
'''
    source=source.replace('using namespace LAMMPS_NS;',timer+'\nusing namespace LAMMPS_NS;')
    signatures=[('void BondBPMPeriVariational::gather(){',0),('void BondBPMPeriVariational::evaluate(){',1),('void BondBPMPeriVariational::compute(int ef,int vf){',2),('int BondBPMPeriVariational::commit_state(long long tok){',4),('void BondBPMPeriVariational::sync_topology(){',5),('void BondBPMPeriVariational::settings(int argc,char **argv){',6),('int BondBPMPeriVariational::output(int field,double *out){',7),('std::string BondBPMPeriVariational::serialize()const{',8),('void BondBPMPeriVariVariational::deserialize(const std::string &text){',9)]
    signatures[-1]=('void BondBPMPeriVariational::deserialize(const std::string &text){',9)
    for sig,slot in signatures:
        assert source.count(sig)==1,sig;source=source.replace(sig,sig+f' V4Timer v4_scope({slot});')
    begin=' for(int b=0;b<nb;b++){nextk[b]=';end=' bool changed=nextd!=d;'
    assert source.count(begin)==source.count(end)==1
    source=source.replace(begin,' {V4Timer candidate_scope(3);\n'+begin).replace(end,' }\n'+end)
    cpp=build/'bond_bpm_peri_variational.cpp';cpp.write_text(source)
    oldbuild=json.loads((V2/'BUILD_AUDIT.json').read_text());objects=list(oldbuild['object_sha256'])
    for path,want in oldbuild['object_sha256'].items():assert sha(path)==want,path
    flags=oldbuild['commands'][0][:oldbuild['commands'][0].index('-c')]
    flags[1:1]=['-I'+str(O/'provenance/v3_source')]
    command=flags+['-c',str(cpp),'-o',str(build/'bond_bpm_peri_variational.cpp.o')]
    with (build/'build.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        matches=[i for i,p in enumerate(objects) if Path(p).name=='bond_bpm_peri_variational.cpp.o'];assert len(matches)==1
        objects[matches[0]]=str(build/'bond_bpm_peri_variational.cpp.o')
        rsp=build/'objects.rsp';rsp.write_text('\n'.join(objects)+'\n')
        link=oldbuild['commands'][-1].copy();link[link.index('-o')+1]=str(build/'liblammps.so.0');link=[('@'+str(rsp)) if x.startswith('@') else x for x in link]
        subprocess.run(link,cwd=oldbuild['link_cwd'],stdout=log,stderr=subprocess.STDOUT,check=True)
    lib=build/'liblammps.so'
    if not lib.exists():lib.symlink_to('liblammps.so.0')
    dump(O/'P0_TIMING_BUILD.json',{'status':'PASS','physics_changed':False,'timing_only':True,'source_sha256':sha(cpp),'library_sha256':sha(lib),'baseline_library_sha256':sha(V2/'build/liblammps.so'),'commands':[command,link],'native_profile_slots':['gather','evaluate','compute','damage_candidate','history_commit','sync_topology','initialize','state_output','serialize','deserialize']})
    print('P0_TIMING_BUILD_PASS',flush=True)
if __name__=='__main__':main()
