"""Distinguish original scientific source identity from relocated runtime source identity."""
from pathlib import Path
import hashlib, json, zipfile
from .paths import PROJECT_ROOT, project_path

def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def protected_originals():
    # Evidence is read as bytes only. Nothing in archive is imported or executed.
    with zipfile.ZipFile(project_path('inputs/reference/source_identity/original_source_identity.zip')) as z:
        manifest = json.loads(z.read('protected_manifest.json'))
        for name, expected in manifest.items():
            if hashlib.sha256(z.read('protected/'+name)).hexdigest() != expected:
                raise ValueError('Original protected source identity mismatch: '+name)
    return manifest

def runtime_identity():
    return {str(p.relative_to(PROJECT_ROOT)):sha256(p)
            for p in sorted((PROJECT_ROOT/'mouse_vessel_flow').rglob('*'))
            if p.suffix in ('.py','.cpp')}

def verify_source_integrity():
    manifest=json.loads(project_path('configs/runtime_source_identity.json').read_text())
    current=runtime_identity()
    if current!=manifest:
        raise ValueError('Runtime sources changed; audit and explicitly version the new run identity')
    return current
