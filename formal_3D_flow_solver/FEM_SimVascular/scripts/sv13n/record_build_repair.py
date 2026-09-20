"""Record the two known build-script corrections without rerunning any build.

The before-side is a minimal reconstruction of the rejected expressions from
the retained failure evidence, not a claimed byte-for-byte historical snapshot.
"""
import difflib
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / 'reports/sv1_3n'
name = 'scripts/sv13n/petsc_build_remote.py'
after = (ROOT / name).read_text()
before = after.replace("'CUDAFLAGS=-ccbin '", "'CUDAC_FLAGS=-ccbin '")
before = before.replace("'run successfully with cuda' in text(check).lower()",
                        "'run successfully with cuda' in text(check)")
assert before != after
patch = ''.join(difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                                   fromfile='a/' + name, tofile='b/' + name))
completion = 'scripts/sv13n/complete_gpu_build_remote.py'
patch += ''.join(difflib.unified_diff([], (ROOT / completion).read_text().splitlines(True),
                                    fromfile='/dev/null', tofile='b/' + completion))
target = ROOT / 'patches/sv1_3n/build_config_parser_repair.patch'
target.write_text(patch)
d = json.loads((R / 'compatibility_adapter.json').read_text())
repair = next(x for x in d['repairs'] if x['id'] == 'repair_02')
repair.update(
    patch=str(target.relative_to(ROOT)),
    patch_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
    added_lines=sum(x.startswith('+') and not x.startswith('+++') for x in patch.splitlines()),
    removed_lines=sum(x.startswith('-') and not x.startswith('---') for x in patch.splitlines()),
    diff_hunks=[x for x in patch.splitlines() if x.startswith('@@')],
    before_side_provenance='Minimal reconstruction of the two rejected expressions from retained configure/parser evidence; not a complete original-script snapshot.',
    current_file_sha256={n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest()
                         for n in (name, completion)},
)
(R / 'compatibility_adapter.json').write_text(json.dumps(d, indent=2) + '\n')
print('repair_02 patch and affected-line counts recorded; no build or solver run.')
