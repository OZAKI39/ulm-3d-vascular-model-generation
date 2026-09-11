"""Static provenance and dependency records; never imports project code."""
from scope import *

def cache_record(p):
    result=dict(path=str(p),exists=p.is_file())
    if p.is_file():
        result['sha256']=sha(p)
        result['values']={m[0]:m[1] for m in re.findall(r'^(CMAKE_HOME_DIRECTORY|CMAKE_PROJECT_NAME|CMAKE_BUILD_TYPE):[^=]*=(.*)$',p.read_text(),re.M)}
    return result

def native_sources(items):
    m=ROOTS['mirheo_starter'];h=ROOTS['hemocell_starter']
    build=json.loads((B/'isolated_build.json').read_text());spec=json.loads((RUN/'spec.json').read_text())
    lib=Path(build['library']);cache=cache_record(B/'native/build/CMakeCache.txt')
    confirmed=(build['source_directory']==str(B/'native/source')==cache['values']['CMAKE_HOME_DIRECTORY'] and
        spec['native_library']==str(lib) and spec['native_library_sha256']==build['library_sha256']==sha(lib) and
        build['patch_sha256']==sha(Path(build['patch'])))
    if not confirmed:raise RuntimeError('Latest native build/source evidence mismatch; preserve evidence and investigate')
    roots=[(m/'vendor/Mirheo',[m/'vendor/Mirheo/build',m/'metadata'],'Old installed WSL library; not the latest isolated repair'),
        (B/'native/source',[B/'native/build'],'CONFIRMED: latest explicit A6 run and isolated_build.json agree with actual library SHA-256 and CMake home'),
        (h/'vendor/HemoCell',[h/'build/upstream'],'HemoCell library source, confirmed CMake home'),
        (h/'vendor/HemoCell/palabos',[],'Extracted Palabos source; version is a reported archive revision, not the parent HemoCell HEAD')]
    rows=[]
    for p,builds,relationship in roots:
        prefix=p.relative_to(PROJECTS).as_posix()+'/'
        row=git_info(p);row.update(build_directories=[str(x) for x in builds],build_evidence=[cache_record(x/'CMakeCache.txt') for x in builds],
            included_files=sum(k.startswith(prefix) for k in items),included=True,current_runtime_relationship=relationship,
            source_copy_location='source/'+p.relative_to(PROJECTS).as_posix())
        row['contains_uncommitted_or_untracked_changes']=bool(row['status_porcelain']) if row['is_exact_git_root'] else None
        rows.append(row)
    rows[-1]['archive_source_record']=json.loads((h/'metadata/sources.json').read_text())
    subs=[]
    for sub in ['src/extern/backward-cpp','src/extern/cub','src/extern/cuda_variant','src/extern/pugixml','src/extern/pybind11','units/extern/googletest']:
        original=m/'vendor/Mirheo'/sub;isolated=B/'native/source'/sub
        def signatures(root):
            prefix=root.relative_to(PROJECTS).as_posix()+'/'
            return {k[len(prefix):]:sha(v['source_path']) for k,v in items.items() if k.startswith(prefix) and Path(v['source_path']).is_file()}
        a=signatures(original);b=signatures(isolated)
        subs.append(dict(relative=sub,original=git_info(original),original_files=len(a),isolated_files=len(b),
            actual_contents_copied=bool(a) and bool(b),only_original=sorted(set(a)-set(b)),only_isolated=sorted(set(b)-set(a)),
            differing_files=[k for k in sorted(set(a)&set(b)) if a[k]!=b[k]],
            note='Isolated submodule Git status may be uninitialized because Git metadata was omitted from the prior copy; actual source files are present. Two missing isolated pugixml Windows NuGet packaging files are preserved in original vendor only.'))
    return dict(status='CONFIRMED_LATEST_NATIVE_SOURCE',latest_source=str(B/'native/source'),
        latest_library=dict(path=str(lib),sha256=sha(lib),included=False,reason='Compiled WSL binary intentionally excluded'),
        installed_library=dict(path=build['installed_library'],sha256=sha(build['installed_library']),included=False),
        explicit_delivery_receipt=dict(path=str(B/'delivery_receipt.json'),sha256=sha(B/'delivery_receipt.json')),
        latest_spec=dict(path=str(RUN/'spec.json'),sha256=sha(RUN/'spec.json')),native_trees=rows,submodules=subs,
        ambiguity=[],external_sources_required=False,untracked_header=str(B/'native/source/src/mirheo/core/bouncers/repair_trace.h'),
        case_sources=[dict(path=str(h/'cases'/n),included=True,cache=cache_record(h/'build'/n/'CMakeCache.txt')) for n in ['single_rbc_shear_benchmark','pure_fluid_benchmark']])

def inputs(items,remote_dir):
    m=ROOTS['mirheo_starter'];h=ROOTS['hemocell_starter']
    spec=json.loads((RUN/'spec.json').read_text());old=m/'runs/single_rbc_benchmark'/CONFIG['legacy_campaign']/'solver/main_hemocell_1'
    mandatory=[RUN/'spec.json',Path(spec['mesh']),B/'native/local_before_halo_v2.patch',B/'native/source/src/mirheo/core/bouncers/repair_trace.h',
        B/'native/source/src/mirheo/core/simulation.cpp',B/'native/source/src/mirheo/core/bouncers/from_mesh.cu',
        m/'py_scripts/single_rbc_benchmark_repaired.yaml',m/'py_scripts/single_rbc_repair/continuous_worker.py',
        m/'py_scripts/single_rbc_repair/continuous_protocol.py',m/'py_scripts/single_rbc_repair/half_step_control.py',
        B/'frozen_benchmark_plan.json',B/'authorization.json',B/'material_matching.json',B/'A5_failure_step_geometry.json',
        B/'runtime_reviews/A5_continuous_preparation_30/review/review.json',h/'cases/single_rbc_shear_benchmark/benchmark.cpp',
        h/'cases/single_rbc_shear_benchmark/CMakeLists.txt',old/'config.xml',old/'RBC.xml',old/'RBC.pos',old/'reference.off']
    rows=[]
    for p in mandatory:
        key=p.relative_to(PROJECTS).as_posix()
        rows.append(dict(path=str(p),source_relative=key,status='INCLUDED' if key in items and p.is_file() and p.stat().st_size else 'MISSING',sha256=sha(p) if p.is_file() else None))
    if any(x['status']=='MISSING' for x in rows):raise RuntimeError('Missing necessary source/input: '+str([x['path'] for x in rows if x['status']=='MISSING']))
    unresolved=[
        dict(status='NOT_UPLOADED_BY_DESIGN',path=spec['native_library'],role='Mirheo compiled module',reason='Requires a separate cloud build and path/hash adaptation; no binary or build tree transfer'),
        dict(status='NOT_UPLOADED_BY_DESIGN',path=str(h/'build/upstream/libhemocell.a'),role='HemoCell compiled library',reason='CMake references this generated library; cloud build not authorized in this transfer'),
        dict(status='UNVERIFIED',path=str(m/'.venv'),role='Local Python, MPI, CUDA and library runtime',reason='Original .venv excluded. Existing /workspace/bloodflow/.venv preserved and fingerprinted, but no dependency import/install or runtime validation'),
        dict(status='UNVERIFIED',role='Cloud configuration and authorization',reason='WSL absolute paths and historical library/source protection hashes remain unchanged. Historical local GPU/build ledgers are evidence, not authorization to execute on Vast.ai'),
        dict(status='UNVERIFIED',role='Historical diagnostics/report regeneration',reason='Only selected run configs and provenance copied; raw profiles/frames/bounce traces under runs and derived analysis_cache omitted. Analysis entry points requiring these histories cannot be asserted runnable'),
        dict(status='UNVERIFIED',role='Dynamic input/output paths',reason='Static reading only; generated same_process_relaxed_geometry.off and runtime result paths are not pre-existing mandatory cold-start inputs. Arbitrary CLI choices and all historical provenance references were not expanded')]
    mapping=dict(original_content_rewritten=False,roots=[dict(local=str(root),staging='source/'+name,remote=remote_dir+'/source/'+name) for name,root in ROOTS.items()],
        existing_remote_environment='/workspace/bloodflow/.venv',external_sources=[],external_required_inputs=[],
        instruction='Mapping is documentation only. Do not use historical WSL absolute paths as verified cloud paths. Deployment adaptation is a separate task.')
    report=dict(status='TRANSFER_INPUTS_PRESENT_RUNTIME_UNVERIFIED',static_entries_reviewed=[str(m/'py_scripts/single_rbc_repair'/n) for n in ['continuous_worker.py','continuous_protocol.py','half_step_control.py','workflow.py']]+[str(old/'config.xml'),str(old/'RBC.xml'),str(h/'cases/single_rbc_shear_benchmark/CMakeLists.txt')],
        resolved=rows,unresolved=unresolved,required_external_files=[],missing_known_source_or_mesh_files=[],
        historical_hemocell_input_example=str(old),historical_example_is_new_cloud_parameter_freeze=False,
        file_transfer_completeness_and_cloud_runtime_readiness_are_distinct=True)
    return report,mapping,[x['source_relative'] for x in rows]
