from common import *
import subprocess,re,difflib

def main():
 old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');rt=HERE/'runtime';rt.mkdir(exist_ok=True)
 previous=read(PREVIOUS/'COMMON_GEOMETRY_PATCH_IDENTITY.json')
 for p,h in previous['reused_objects'].items():assert sha(p)==h
 assert sha(PREVIOUS/'runtime/patched.o')==previous['objects']['runtime/patched.o']['object_sha256']
 cmd0=[a.replace(old['source_checkout'],str(RETAINED)) for a in old['compile_command']]+['-I'+str(HERE),'-I'+str(PREVIOUS)]
 meta=dict(compile_commands=[],link_commands=[],objects={},binaries={},production_modified=False,patch_identity_sha256=sha(PREVIOUS/'COMMON_GEOMETRY_PATCH_IDENTITY.json'))
 def compile(src,obj):
  cmd=list(cmd0);cmd[cmd.index('-c')+1]=str(src);cmd[cmd.index('-o')+1]=str(obj)
  meta['compile_commands'].append(cmd);put('DISPATCH_BUILD_IDENTITY.json',meta)
  subprocess.run(cmd,cwd=old['link_cwd'],check=True)
  meta['objects'][str(obj.relative_to(HERE))]=dict(source_sha256=sha(src),object_sha256=sha(obj))
 compile(HERE/'native_dispatch.cpp',rt/'native_dispatch.o')
 targets={'set_bc.cpp':{'set_bc_dir':'','set_bc_dir_l':'lFa.name','set_bc_dir_w':'','set_bc_dir_wl':'lFa.name','set_bc_neu':'','set_bc_neu_l':'lFa.name','set_bc_trac_l':'lFa.name'},'nn.cpp':{'get_nnx':'std::to_string(eNoN)'},'eq_assem.cpp':{'b_assem_neu_bc':'lFa.name'}}
 patch=''
 for filename,functions in targets.items():
  original=(SOURCE/filename).read_text();text=original;off=0
  for fn,face in functions.items():
   match=re.search(r'\bvoid '+fn+r'\([^;]*?\)\s*\{',text,re.S)
   assert match,(filename,fn)
   at=match.end();stmt='\n  dispatch_trace::hit("'+fn+'"'+(','+face if face else '')+');'
   text=text[:at]+stmt+text[at:]
  text='#include "dispatch_trace.h"\n'+text;(rt/filename).write_text(text)
  patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),fromfile='a/'+filename,tofile='b/'+filename))
  compile(rt/filename,rt/(filename+'.o'))
 (HERE/'DISPATCH_OBSERVATION.diff').write_text(patch)
 for traced in [False,True]:
  link=list(old['link_command'])
  for i,a in enumerate(link):
   if a.endswith('/main.cpp.o'):link[i]=str(rt/'native_dispatch.o')
   if a.endswith('/fluid.cpp.o'):link[i]=str(PREVIOUS/'runtime/patched.o')
   if traced:
    for f in targets:
     if a.endswith('/'+f+'.o'):link[i]=str(rt/(f+'.o'))
  binary=rt/('dispatch_trace' if traced else 'dispatch_original');link[link.index('-o')+1]=str(binary)
  meta['link_commands'].append(link);put('DISPATCH_BUILD_IDENTITY.json',meta)
  subprocess.run(link,cwd=old['link_cwd'],check=True);meta['binaries'][binary.name]={'path':str(binary),'sha256':sha(binary)}
 meta['status']='BUILT';put('DISPATCH_BUILD_IDENTITY.json',meta)
if __name__=='__main__':main()
