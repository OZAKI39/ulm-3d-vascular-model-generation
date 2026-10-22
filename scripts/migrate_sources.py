"""Extract local dependency closures; archive original bytes and record every relocation."""
from pathlib import Path
import ast, hashlib, json, re, zipfile, difflib

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;CE=BASE/'computation_examples/_shared'
P=CE/'particle_workflow/particle_3d/src/particle_3d'
N=CE/'active_network_development/ulm_3D_vascular'
F=CE/'formal_3D_flow_solver/FEM_SimVascular'
V=BASE/'vascular_printing'
selected={}; mapping=[]

def add(src,dst):
    src=Path(src);dst=Path(dst)
    if dst in selected and selected[dst].resolve()!=src.resolve():
        if selected[dst].read_bytes()!=src.read_bytes():raise ValueError(('conflicting sources',dst,src))
    selected[dst]=src

def particle_closure():
    seen=set();stack=['formal_dynamics_p9a5','population_inlet_p9a4','formal_cohort_p9a5',
        'particle82_point_native','field','rbc_distribution','rbc_integrator','microbubble','validation_cases']
    while stack:
        name=stack.pop();p=P/(name+'.py')
        if name in seen or not p.is_file():continue
        seen.add(name);add(p,Path('particles')/p.name)
        for a in ast.walk(ast.parse(p.read_text())):
            if isinstance(a,ast.ImportFrom) and a.level==1 and a.module:stack.append(a.module.split('.')[0])
    add(P/'particle82_point_native.cpp','particles/particle82_point_native.cpp')
    return seen

def geometry_closure():
    roots={'utils':V/'utils','vascular_processing':V/'vascular_processing'}
    # Flow preprocessing/surface code authority is the active network workspace.
    special={'utils.cfd_preprocess':N/'utils/cfd_preprocess','utils.cfd_surface_prepare':N/'utils/cfd_surface_prepare'}
    def source(name):
        for key,root in sorted(special.items(),key=lambda x:-len(x[0])):
            if name==key or name.startswith(key+'.'):
                p=root.joinpath(*name[len(key):].strip('.').split('.')) if name!=key else root
                break
        else:
            first=name.split('.')[0]
            if first not in roots:return None
            p=roots[first].joinpath(*name.split('.')[1:])
        if p.is_dir():p=p/'__init__.py'
        else:p=p.with_suffix('.py')
        return p if p.is_file() else None
    stack=['utils.rodent_vasculature','utils.sampling.pipeline','utils.swc_roi_yaml_config',
           'utils.cfd_lumen.ultraliser_backend','utils.cfd_lumen','utils.cfd_preprocess','utils.cfd_surface_prepare']
    seen=set()
    while stack:
        name=stack.pop();p=source(name)
        if p is None or name in seen:continue
        seen.add(name)
        dest=name.replace('utils','geometry_tools',1) if name.startswith('utils') else name
        rel=Path(*dest.split('.'))
        rel=rel/'__init__.py' if p.name=='__init__.py' else rel.with_suffix('.py')
        add(p,rel)
        package=name if p.name=='__init__.py' else name.rpartition('.')[0]
        for a in ast.walk(ast.parse(p.read_text())):
            if isinstance(a,ast.ImportFrom):
                prefix='.'.join(package.split('.')[:len(package.split('.'))-a.level+1]) if a.level else ''
                target='.'.join(x for x in (prefix,a.module) if x)
                stack.append(target)
                stack.extend(target+'.'+x.name for x in a.names)
            elif isinstance(a,ast.Import):stack.extend(x.name for x in a.names)
        for i in range(1,len(name.split('.'))):stack.append('.'.join(name.split('.')[:i]))
    return seen

def relocate(text):
    replacements={
        'particle_3d':'mouse_vessel_flow.particles',
        'network_1d0d':'mouse_vessel_flow.network_models',
        'flow_solver_support':'mouse_vessel_flow.solver_support',
        'vascular_validation':'mouse_vessel_flow.mesh_tools',
        'utils':'mouse_vessel_flow.geometry_tools',
        'vascular_processing':'mouse_vessel_flow.vascular_processing',
    }
    # Rewrite import syntax only; never scientific metadata/string identifiers.
    for old,new in replacements.items():
        text=re.sub(r'(?m)^(\s*from )'+re.escape(old)+r'(?=[.\s])',r'\g<1>'+new,text)
        text=re.sub(r'(?m)^(\s*import )'+re.escape(old)+r'(?=[.\s])',r'\g<1>'+new,text)
    text=text.replace("Path(__file__).resolve().parents[3]", "PROJECT_ROOT")
    if 'PROJECT_ROOT' in text and 'from mouse_vessel_flow.paths import PROJECT_ROOT' not in text:
        # Insert after future imports to preserve syntax.
        lines=text.splitlines(True);tree=ast.parse(text)
        pos=tree.body[0].end_lineno if tree.body and isinstance(tree.body[0],ast.Expr) and isinstance(tree.body[0].value,ast.Constant) else 0
        for node in tree.body:
            if isinstance(node,ast.ImportFrom) and node.module=='__future__':pos=node.end_lineno
        lines.insert(pos,'from mouse_vessel_flow.paths import PROJECT_ROOT\n');text=''.join(lines)
    for old,new in [
        ('formal_3D_flow_solver/FEM_SimVascular','inputs/reference/fem'),
        ('particle_3d/reports/particle9a5_formal_trajectories','results/network_h0/particles/particle9a5_formal_trajectories'),
        ('particle_3d/reports/','inputs/reference/particles/'),
        ('particle_3d/contracts/','configs/particles/contracts/'),
        ('particle_3d/outputs/','results/development/'),
    ]:text=text.replace(old,new)
    text=text.replace("Path('/home/lzy/projects/sonovue_size_distribution_v0')", "PROJECT_ROOT/'inputs/reference/sonovue'")
    text=text.replace('Path(__file__).resolve().parents[2] / "contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json"',
                      'Path(__file__).resolve().parents[2] / "configs/particles/contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json"')
    return text

def run():
    particle_closure();geo=geometry_closure()
    for p in (N/'network_1d0d').glob('*.py'):add(p,Path('network_models')/p.name)
    for old,new in [(F/'solver_support/src/flow_solver_support','solver_support'),(F/'mesh_generate/src/vascular_validation','mesh_tools')]:
        for p in old.glob('*.py'):add(p,Path(new)/p.name)
    add(CE/'sonovue_size_distribution_v0/src/sonovue_sampler.py','particles/sonovue_sampler.py')
    source_zip=ROOT/'inputs/reference/source_identity/original_source_identity.zip'
    protected=json.loads((CE/'particle_workflow/particle_3d/reports/network_derived_flow_mb_validation_v1/data/code_snapshot_manifest.json').read_text())
    with zipfile.ZipFile(source_zip,'w',zipfile.ZIP_DEFLATED) as z:
        for relative,expected in protected.items():
            source=CE/'particle_workflow'/relative;data=source.read_bytes()
            assert hashlib.sha256(data).hexdigest()==expected,relative
            z.writestr('protected/'+relative,data)
        z.writestr('protected_manifest.json',json.dumps(protected,indent=2))
        for dst,src in selected.items():
            original=src.read_bytes();text=original.decode();new=relocate(text) if src.suffix=='.py' else text
            target=ROOT/'mouse_vessel_flow'/dst;target.parent.mkdir(parents=True,exist_ok=True);target.write_text(new)
            z.writestr('selected/'+str(dst),original)
            names=[]
            if src.suffix=='.py':
                names=[n.name for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef))]
            mapping.append(dict(source=str(src),path=str(target.relative_to(ROOT)),original_sha256=hashlib.sha256(original).hexdigest(),
                sha256=hashlib.sha256(target.read_bytes()).hexdigest(),symbols=names,changed=original!=target.read_bytes()))
            if original!=target.read_bytes():
                diff=''.join(difflib.unified_diff(text.splitlines(True),new.splitlines(True),fromfile=str(src),tofile=str(dst)))
                d=ROOT/'docs/audit/source_diffs'/dst.with_suffix('.diff');d.parent.mkdir(parents=True,exist_ok=True);d.write_text(diff)
    for d in (ROOT/'mouse_vessel_flow').rglob('*'):
        if d.is_dir() and not (d/'__init__.py').exists():(d/'__init__.py').write_text('"""Mouse vessel workflow components."""\n')
    (ROOT/'mouse_vessel_flow/__init__.py').write_text('"""Independent mouse vessel flow and microbubble workflow."""\n__version__ = "0.1.0"\n')
    (ROOT/'docs/audit/source_map.json').write_text(json.dumps(mapping,indent=2)+'\n')
    (ROOT/'docs/audit/geometry_dependency_closure.json').write_text(json.dumps(sorted(geo),indent=2)+'\n')
    print('SOURCE_MODULES',len(mapping),'geometry',len(geo),'protected_original',len(protected))

if __name__=='__main__':run()
