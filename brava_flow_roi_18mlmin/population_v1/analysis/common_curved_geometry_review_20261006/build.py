"""Patch only an independent copy of fluid.cpp, retaining the exact production objects."""
from common import *
import difflib
import subprocess
import shutil

def patched_source(original):
    inc='#include "fluid.h"'
    text=original.replace(inc,inc+'\n#ifdef SVMP_VALIDATION_COMMON_GEOMETRY\n#include "common_geometry.h"\n#endif',1)
    begin=text.index('void construct_fluid(');end=text.index('/// @brief Reproduces Fortran',begin)
    f=text[begin:end]
    for opt in [1,2]:
        loop=f'    for (int g = 0; g < fs[{opt-1}].nG; g++) {{'
        pos=f.index(loop)
        start=f.index('      if (g == 0 || !fs[',pos)
        stop=f.index(f'      double w = fs[{opt-1}].w(g) * Jac;',start)
        old=f[start:stop]
        common='''#ifdef SVMP_VALIDATION_COMMON_GEOMETRY
      if (nsd == 3 && lM.nFs == 2 && lM.eType == ElementType::TET10) {
        auto Dgeom = fs[0].Nx.rslice(g);
        const auto geometry = validation_common_geometry::evaluate(fs[0].N.rcol(g), Dgeom, xwl);
        validation_common_geometry::transform(Dgeom, geometry, Nwx);
        const auto Dpressure = fs[1].Nx.rslice(g);
        validation_common_geometry::transform(Dpressure, geometry, Nqx);
        Jac = geometry.determinant;
        ksix = geometry.metric;
        auto Hreference = fs[0].Nxx.rslice(g);
        // Same mathematical chain rule; reuse this geometry's inverse and
        // centered curvature to avoid cancellation in uncentered coordinates.
        validation_common_geometry::transform_hessian(Hreference, xwl, geometry, Nwx, Nwxx);
      } else
#endif
      {
'''+old+'      }\n'
        f=f[:start]+common+f[stop:]
    return text[:begin]+f+text[end:]

def instrument(text):
    text=text.replace('#include "fluid.h"','#include "fluid.h"\n#include "trace.h"',1)
    for opt in [1,2]:
        needle=f'      double w = fs[{opt-1}].w(g) * Jac;'
        assert text.count(needle)==1
        text=text.replace(needle,needle+f'\n      validation_trace::sample({opt},g,fs,Jac,ksix,Nwx,Nqx,Nwxx);',1)
    # Capture actual parameters after their unchanged native calculation.
    a=text.index('void fluid_3d_c(');b=text.index('void fluid_3d_m(',a)
    part=text[a:b];needle='  //  Local residual'
    assert needle in part
    part=part.replace(needle,'  validation_trace::tau(2,vmsFlag,tauM,0.0,0.0);\n'+needle,1)
    text=text[:a]+part+text[b:]
    a=text.index('void fluid_3d_m(');part=text[a:]
    needle='  rV[0] = tauB*'
    assert needle in part
    part=part.replace(needle,'  validation_trace::tau(1,vmsFlag,tauM,tauC,tauB);\n\n'+needle,1)
    return text[:a]+part

def main():
    assert read(HERE/'CURVED_MIXED_BLOCKER_BASELINE.json')['status']=='BLOCKER_REPRODUCED'
    old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json')
    rt=HERE/'runtime';rt.mkdir(exist_ok=True)
    source=(SOURCE/'fluid.cpp').read_text();patched=patched_source(source)
    (rt/'fluid_patched.cpp').write_text(patched)
    diff=''.join(difflib.unified_diff(source.splitlines(True),patched.splitlines(True),
        fromfile='a/Code/Source/solver/fluid.cpp',tofile='b/Code/Source/solver/fluid.cpp'))
    diff+=''.join(difflib.unified_diff([], (HERE/'common_geometry.h').read_text().splitlines(True),
        fromfile='/dev/null',tofile='b/Code/Source/solver/common_geometry.h'))
    (HERE/'COMMON_GEOMETRY_PATCH.diff').write_text(diff)
    (rt/'fluid_baseline_trace.cpp').write_text(instrument(source))
    (rt/'fluid_patched_trace.cpp').write_text(instrument(patched))
    cmd0=[a.replace(old['source_checkout'],str(RETAINED)) for a in old['compile_command']]
    cmd0+=['-I'+str(HERE)]
    meta=dict(status='BUILDING',production_source_sha256=sha(SOURCE/'fluid.cpp'),
        protocol_sha256=sha(HERE/'COMMON_GEOMETRY_HARNESS_PROTOCOL.json'),patch_sha256=sha(HERE/'COMMON_GEOMETRY_PATCH.diff'),
        production_modified=False,compile_commands=[],link_commands=[],binaries={},objects={})
    def compile_file(src,obj,patch=False):
        cmd=list(cmd0);cmd[cmd.index('-c')+1]=str(src);cmd[cmd.index('-o')+1]=str(obj)
        if patch:cmd+=['-DSVMP_VALIDATION_COMMON_GEOMETRY']
        meta['compile_commands'].append(cmd);put('COMMON_GEOMETRY_PATCH_IDENTITY.json',meta)
        subprocess.run(cmd,cwd=old['link_cwd'],check=True)
        meta['objects'][str(obj.relative_to(HERE))]=dict(source_sha256=sha(src),object_sha256=sha(obj))
    compile_file(HERE/'native_harness.cpp',rt/'native_harness.o')
    variants=[('original',None,False),('baseline_trace','fluid_baseline_trace.cpp',False),
              ('patched_trace','fluid_patched_trace.cpp',True),('patched','fluid_patched.cpp',True)]
    for name,src,patch in variants:
        if src:compile_file(rt/src,rt/(name+'.o'),patch)
        link=list(old['link_command'])
        for i,a in enumerate(link):
            if a.endswith('/main.cpp.o'):link[i]=str(rt/'native_harness.o')
            if src and a.endswith('/fluid.cpp.o'):link[i]=str(rt/(name+'.o'))
        link[link.index('-o')+1]=str(rt/('harness_'+name))
        meta['link_commands'].append(link);put('COMMON_GEOMETRY_PATCH_IDENTITY.json',meta)
        subprocess.run(link,cwd=old['link_cwd'],check=True)
        meta['binaries'][name]=dict(path=str(rt/('harness_'+name)),sha256=sha(rt/('harness_'+name)))
    # All reused objects explicitly identified; no modification/rebuild in production tree.
    meta['reused_objects']={str((Path(old['link_cwd'])/a).resolve()):sha(Path(old['link_cwd'])/a)
        for a in old['link_command'] if a.endswith(('.o','.a')) and not a.endswith('/main.cpp.o')}
    meta['status']='BUILT';put('COMMON_GEOMETRY_PATCH_IDENTITY.json',meta)
    print('Built original, trace-only baseline, gated patched trace, gated patched native binaries')

if __name__=='__main__':main()
