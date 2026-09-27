import glob, json, runpy
files=glob.glob('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/mesh_generate/scripts/sv_*.py')
for name in files:
    with open(name) as stream: compile(stream.read(), name, "exec")
print("EMBEDDED_SYNTAX_PASS="+json.dumps(files))
runpy.run_path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/mesh_generate/scripts/sv_api_probe.py',run_name="__main__")
