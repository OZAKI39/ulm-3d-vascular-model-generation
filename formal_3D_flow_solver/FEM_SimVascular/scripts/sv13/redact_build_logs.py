#!/usr/bin/env python3
"""Remove environment credentials from stage build logs before review/export."""
import json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,write_json
pattern=re.compile(r'(?:TOKEN|PASSWORD|PASSWD|SECRET|API_KEY|ACCESS_KEY|AUTHORIZATION|CREDENTIAL|OPEN_BUTTON)',re.I)
manifest=ROOT/'reports/sv1_3/log_redaction.json'
data=json.loads(manifest.read_text()) if manifest.exists() else {'reason':'PETSc configure environment dumps may contain credentials','files':[]}
for folder in (ROOT/'outputs/sv1_3').glob('remote_build*'):
    for p in folder.rglob('*'):
        if not p.is_file() or p.is_symlink() or p.suffix not in ('.log','.txt'):continue
        lines=p.read_text(errors='replace').splitlines(keepends=True)
        remove=lambda line: bool(pattern.search(line)) and ('=' in line or ':' in line)
        count=sum(remove(line) for line in lines)
        if count:
            before=sha256(p)
            p.write_text(''.join('[REDACTED environment credential line]\n' if remove(line) else line for line in lines))
            data['files'].append({'path':str(p.relative_to(ROOT)),'redacted_lines':count,
                                  'original_sha256':before,'redacted_sha256':sha256(p)})
write_json(manifest,data)
print('Build log credential redaction complete',flush=True)
