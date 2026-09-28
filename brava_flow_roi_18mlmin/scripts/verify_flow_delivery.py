"""Verify transferred field identity and decode every delivered flow movie."""
from pathlib import Path
import hashlib,json
import imageio_ffmpeg

ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
active=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text())
assert active['status']=='PASS'
case=ROOT/active['case']
assert sha(case/'frozen_flow/steady_flow.vtu')==active['flow_sha256']
input_hashes=json.loads((case/'input_hashes.json').read_text())
for relative,digest in input_hashes.items():assert sha(case/relative)==digest,relative
movies=[]
for pose in ['0','15']:
    folder=ROOT/'visualization'/('candidate_'+pose)
    manifest=json.loads((folder/'MEDIA_VALIDATION.json').read_text())
    assert manifest['source_flow_sha256']==active['flow_sha256']
    assert len(manifest['media'])==5
    for row in manifest['media']:
        movie=folder/row['file']
        assert sha(movie)==row['sha256'],str(movie)
        reader=imageio_ffmpeg.read_frames(str(movie),pix_fmt='rgb24')
        header=next(reader);count=sum(1 for _ in reader)
        assert header['size']==(1920,1080) and header['fps']==24 and count==432
        movies.append(dict(pose=pose,kind=row['kind'],decoded_frames=count,sha256=row['sha256']))
        print('Verified',pose,row['kind'],count,flush=True)
record=dict(PASS=True,flow_sha256=active['flow_sha256'],input_files_verified=len(input_hashes),movies=movies,
    scope='Transferred inputs, flow hash and media decoding; numerical field checks are recorded separately.')
(ROOT/'reports/LOCAL_FLOW_VERIFICATION.json').write_text(json.dumps(record,indent=2)+'\n')
