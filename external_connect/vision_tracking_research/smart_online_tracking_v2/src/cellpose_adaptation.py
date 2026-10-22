"""Load an explicitly admitted microscope checkpoint; never overwrite cyto3."""
import json
from pathlib import Path
from annotation_store import sha


def checkpoint_for(root, force_original=False):
    root=Path(root);path=root/'data/cellpose_profile.json'
    if force_original or not path.exists():
        checkpoint=root/'models/cellpose/cyto3'
        return dict(id='ORIGINAL_CYTO3',checkpoint=str(checkpoint),checkpoint_sha256=sha(checkpoint),admitted=False)
    profile=json.loads(path.read_text())
    if profile.get('schema')!='MICROSCOPE_CELLPOSE_PROFILE_V1' or not profile.get('admitted'):
        raise ValueError('Cellpose profile is not admitted; original labels remain unchanged')
    checkpoint=Path(profile['checkpoint'])
    if not checkpoint.is_file() or sha(checkpoint)!=profile['checkpoint_sha256']:
        raise ValueError('Cellpose microscope checkpoint checksum mismatch')
    return profile
