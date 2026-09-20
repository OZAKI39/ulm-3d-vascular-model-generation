"""Copy an audited, bounded WSL source snapshot; never execute a solver.

Run only in a fresh handoff branch. Existing destination files must match.
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    target = Path(__file__).resolve().parents[2]
    repo = target.parents[1]
    git = lambda *a: subprocess.check_output(['git', '-C', str(repo), *a], text=True).strip()
    assert git('branch', '--show-current').startswith('sync/fem-simvascular-stage-q-particle-handoff-')
    assert git('rev-parse', 'HEAD') == git('rev-parse', 'origin/main')
    meta = target / 'sync_metadata'
    write(meta / 'main_base.json', dict(repository='OZAKI39/ulm-3d-vascular-model-generation',
          remote_url=git('remote', 'get-url', 'origin'), base_branch='main',
          base_commit=git('rev-parse', 'origin/main'), base_commit_date=git('show', '-s', '--format=%cI', 'origin/main'),
          sync_time=datetime.now(timezone.utc).isoformat(), source_FEM_path=str(source),
          branch=git('branch', '--show-current')))
    rows, omitted = [], []
    excluded_parts = {'.git', '__pycache__', '.pytest_cache', '.venv', 'CMakeFiles', 'plot_cache'}

    def copy(p, dest=None):
        relative = p.relative_to(source).as_posix()
        dest = dest or relative
        out = target / dest
        digest = sha(p)
        assert p.stat().st_size < 95 * 2**20, relative
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists():
            assert sha(out) == digest, dest
        else:
            shutil.copy2(p, out)
        assert sha(out) == digest
        rows.append(dict(source_path=relative, relative_path=dest, size_bytes=p.stat().st_size, sha256=digest))

    # All project code, tests, configurations, patches, and small input contracts.
    for top in ['src', 'scripts', 'tests', 'configs', 'patches', 'inputs', 'benchmarks']:
        for p in sorted((source / top).rglob('*')):
            if p.is_file() and not (excluded_parts & set(p.parts)) and p.suffix != '.pyc':
                copy(p)
    copy(source / 'README.md', 'reports/SOURCE_README.md')
    copy(source / 'pyproject.toml')
    for top in ['reports', 'logs', 'outputs']:
        for p in sorted((source / top).rglob('*')):
            if not p.is_file() or excluded_parts & set(p.parts):
                continue
            rel = p.relative_to(source).as_posix()
            size = p.stat().st_size
            selected = False
            if top == 'reports':
                selected = size < 2 * 2**20 or ('/sv1_3q/' in rel and size < 50 * 2**20)
                # These inventories describe machines/toolchains, not required run evidence.
                selected &= p.name not in {'final_environment.json', 'pre_install_environment.json'}
            elif top == 'logs':
                selected = size < 2**20 or ('/sv1_3q/remote/' in rel and p.suffix == '.log')
            else:
                selected = p.suffix in {'.xml', '.json', '.txt', '.dat', '.flow', '.yaml'} and size < 2 * 2**20
                selected &= '/native/' not in rel and '/plot_cache/' not in rel
                selected |= rel.startswith(('outputs/sv1/SV_MESH/', 'outputs/sv1/model/'))
                selected &= p.name not in {'CMakeCache.txt', 'controller_heartbeat'}
            if selected:
                copy(p)
            else:
                omitted.append(dict(source_path=rel, size_bytes=size, sha256=sha(p),
                    reason='Historical/machine inventory, dependency/build binary, or redundant intermediate output; retained in WSL'))

    # Preserve the complete accepted svMultiPhysics source tree, including licenses.
    upstream = source / 'external/sv13q/svMultiPhysics-reuse'
    expected = json.loads((source / 'reports/sv1_3q/source_patch.json').read_text())['after']
    attributes = {}
    pointers = []
    for rel, digest in expected.items():
        assert sha(upstream / rel) == digest, rel
        dest = 'vendor/svMultiPhysics_stage_q/' + rel
        data = (upstream/rel).read_bytes()
        if len(data)<1000 and data.startswith(b'version https://git-lfs.github.com/spec/v1\n'):
            pointers.append(dict(source_path=(upstream/rel).relative_to(source).as_posix(),relative_path=dest,
                size_bytes=len(data),sha256=digest,upstream_relative_path=rel,original_pointer_text=data.decode(),
                lfs_oid_sha256=re.search(rb'oid sha256:([a-f0-9]{64})',data)[1].decode(),
                unfetched_upstream_object_size_bytes=int(re.search(rb'\nsize (\d+)',data)[1]),
                classification='UNFETCHED_UPSTREAM_EXAMPLE_OR_INSTALLER_NOT_PARTICLE_INPUT'))
            continue
        if Path(rel).name == '.gitattributes':
            attributes[rel] = 'sync_metadata/upstream_gitattributes.txt' if rel == '.gitattributes' else 'sync_metadata/upstream_attributes/' + rel + '.txt'
            copy(upstream / rel, attributes[rel])
            (target/dest).parent.mkdir(parents=True, exist_ok=True)
            (target/dest).write_text('# Packaging only: upstream LFS rules are archived in sync_metadata.\n* -text -filter !diff !merge\n')
        else:
            copy(upstream / rel, dest)
    write(meta/'upstream_attribute_map.json', attributes)
    write(meta/'upstream_lfs_pointer_inventory.json',dict(status='ARCHIVED_AS_PROVENANCE_ONLY',files=pointers,
        count=len(pointers),LFS_enabled=False,downloads_performed=False,essential_frozen_artifacts_affected=False))
    upstream_commit = subprocess.check_output(['git', '-C', str(upstream), 'rev-parse', 'HEAD'], text=True).strip()
    patch = subprocess.check_output(['git', '-C', str(upstream), 'diff', '--binary', 'HEAD'])
    # git diff omits untracked files; explicitly include the accepted added header.
    for rel in sorted(set(expected) - set(subprocess.check_output(['git', '-C', str(upstream), 'ls-files'], text=True).splitlines())):
        result = subprocess.run(['git', 'diff', '--no-index', '--', '/dev/null', str(upstream / rel)], capture_output=True)
        assert result.returncode in (0, 1)
        data = result.stdout.replace(str(upstream).encode(), b'')
        patch += data
    (target / 'patches/upstream_to_frozen_stage_q.patch').write_bytes(patch)
    write(meta / 'upstream_source.json', dict(upstream_commit=upstream_commit,
          upstream_remote=subprocess.check_output(['git', '-C', str(upstream), 'remote', 'get-url', 'origin'], text=True).strip(),
          source_tree='vendor/svMultiPhysics_stage_q', verified_files=len(expected),
          unfetched_pointer_inventory='sync_metadata/upstream_lfs_pointer_inventory.json',
          upstream_pointer_files_archived=len(pointers),
          all_inherited_and_stage_q_changes_present=True,
          patch='patches/upstream_to_frozen_stage_q.patch', patch_sha256=sha(target/'patches/upstream_to_frozen_stage_q.patch')))

    run = 'outputs/sv1_3q/REAL_VASCULAR_GPU_ILU_REUSE_WINNER'
    aliases = {
        run + '/1-procs/result_071.vtu': 'frozen_reference/flow/steady_flow_stage_sv1_3q.vtu',
        run + '/1-procs/stFile_071.bin': 'frozen_reference/fem_checkpoint/stFile_071.bin',
        run + '/solver.xml': 'frozen_reference/run/solver.xml',
        run + '/STOP_SIM': 'frozen_reference/run/STOP_SIM',
        run + '/sv13q_stop_ack.json': 'frozen_reference/run/sv13q_stop_ack.json',
    }
    for p in (source / 'outputs/sv1/SV_MESH').rglob('*'):
        if p.is_file():
            aliases[p.relative_to(source).as_posix()] = 'frozen_reference/SV_MESH/' + p.relative_to(source/'outputs/sv1/SV_MESH').as_posix()
    for src, dest in aliases.items():
        copy(source / src, dest)
    write(meta/'copied_from_wsl.json', dict(files=rows, copied_files=len(rows), source=str(source),
          archived_lfs_pointer_files=len(pointers),original_files_verified=len(rows)+len(pointers),
          upstream_pointer_inventory='sync_metadata/upstream_lfs_pointer_inventory.json',
          byte_identical_before_secret_redaction=True, no_remote_download_required=True))
    write(meta/'omitted_artifacts.json', dict(files=omitted, count=len(omitted),
          rule='No dependency/toolkit/build trees, binary executables, caches, credentials; no file >=95 MiB; no LFS or size-hiding archive'))
    print(f'Copied {len(rows)} files; inventoried {len(omitted)} omitted evidence/output files; source tree unchanged.')


if __name__ == '__main__':
    main()
