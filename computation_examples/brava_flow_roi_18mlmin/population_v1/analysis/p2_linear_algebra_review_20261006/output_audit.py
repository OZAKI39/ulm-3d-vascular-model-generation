"""Read-only native postprocessor gate: no assembly, linear solve or timestep."""
from common import *
import subprocess,shutil
old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');rt=HERE/'runtime/output_audit';rt.mkdir(parents=True,exist_ok=True)
# Expose a private method to the audit translation unit only. No class layout,
# method implementation, production header or production object is modified.
h=(SOURCE/'Integrator.h').read_text().replace('private:','public:');(rt/'Integrator.h').write_text(h)
s=(SOURCE/'main.cpp').read_text();needle='    // Create LinearAlgebra objects for each equation.'
s=s.replace(needle,'''    {
      simulation->com_mod.cEq=0;
      std::ofstream o("NATIVE_PRESSURE_CORRECTOR.json");
      try { simulation->get_integrator().corrector_taylor_hood();o<<"{\\\"status\\\":\\\"PASS\\\",\\\"exception\\\":null}\\n"; }
      catch(const std::exception& e){o<<"{\\\"status\\\":\\\"FAIL\\\",\\\"exception\\\":"<<std::quoted(e.what())<<"}\\n";}
      o.close();MPI_Finalize();return 0;
    }
'''+needle);(rt/'main.cpp').write_text(s)
cmd=[a.replace(old['source_checkout'],str(RETAINED)) for a in old['compile_command']]+['-I'+str(SOURCE)];cmd[cmd.index('-c')+1]=str(rt/'main.cpp');cmd[cmd.index('-o')+1]=str(rt/'main.o');
if '--collect-only' not in sys.argv:subprocess.run(cmd,cwd=old['link_cwd'],check=True)
link=list(old['link_command'])
for i,a in enumerate(link):
 if a.endswith('/main.cpp.o'):link[i]=str(rt/'main.o')
 if a.endswith('/fluid.cpp.o'):link[i]=str(PREVIOUS/'runtime/patched.o')
link[link.index('-o')+1]=str(rt/'pressure_corrector_audit');
if '--collect-only' not in sys.argv:subprocess.run(link,cwd=old['link_cwd'],check=True)
run=rt/'run';run.mkdir(exist_ok=True)
if not (rt/'SV_MESH').exists():(rt/'SV_MESH').symlink_to(HERE/'runtime/p2/SV_MESH',target_is_directory=True)
for n in ['solver.xml','inlet_exact_velocity.dat','outlet_exact_traction.vtp']:shutil.copy2(HERE/'runtime/p2_graph/run'/n,run/n)
identity=read(BASE/'production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json');launch=list(identity['command']);launch[launch.index('-n')+1]='1';launch[-2]=str(rt/'pressure_corrector_audit');
if '--collect-only' not in sys.argv:subprocess.run(launch,cwd=run,env=native_env(),check=True)
r=json.loads((run/'NATIVE_PRESSURE_CORRECTOR.json').read_text(),strict=False);r.update(CFD_timesteps=0,linear_solves=0,volume_assembly_calls=0,native_method='Integrator::corrector_taylor_hood',MPI_ranks=1,MPI_role='Serial postprocessor replay only; not solver-rank candidate sweep',pressure_values='initial zero; inversion path is geometry-only and independent of pressure values',source_sha256=sha(SOURCE/'Integrator.cpp'),header_visibility_only=True,binary_sha256=sha(rt/'pressure_corrector_audit'),compile_command=cmd,link_command=link,curved_mesh_sha256=sha(TARGET/'runtime/PIPE_CURVED_ELEVATION.npz'),scope='H high-order output gate only; no geometry/representation or BC validation repeated',production_modified=False)
put('HIGH_ORDER_OUTPUT_GATE.json',r);print('Native pressure output gate',r['status'],r.get('exception'),flush=True)
