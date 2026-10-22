"""Stream saved, checksum-verified full-video records without loading all frames."""
import json
from wcommon import OUT,sha
def full_records():
    for p in sorted((OUT/'results/full_video_chunks').glob('chunk_*.json')):
        if p.name.endswith('.sha.json'):continue
        assert sha(p)==json.loads(p.with_suffix('.sha.json').read_text())['sha256']
        yield from json.loads(p.read_text())
