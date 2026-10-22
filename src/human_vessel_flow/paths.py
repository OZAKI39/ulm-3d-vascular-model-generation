"""Explicit project files and exclusive output directories, independent of cwd."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]

def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def read_json(path):
    return json.loads(Path(path).read_text())

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def local_path(relative, root=ROOT):
    """Resolve a data/config path; old workspace references cannot silently load."""
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError('Project data path escapes this project: ' + str(relative))
    return path

def fresh_directory(path, root=ROOT):
    """Every new run gets a new directory; frozen input/result trees are protected."""
    path = Path(path).resolve()
    for rel in ('inputs', 'cases', 'diagnostics', 'configs', 'native', 'src',
                'results/population', 'results/reports', 'results/visualization'):
        if path.is_relative_to((Path(root)/rel).resolve()):
            raise ValueError('Refuse output inside a frozen or source directory: '+str(path))
    if path.exists():
        raise FileExistsError('A new output directory is required: '+str(path))
    path.mkdir(parents=True, exist_ok=False)
    return path

