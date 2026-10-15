"""Lossless byte-exact archival of completed, no-longer-read rejection NPZs.
Independent of frozen solver; never edits live checkpoints or accepted states.
"""
from common import *
import io,zlib,struct,time,sys,tempfile,signal
MAGIC=b'V3_REJECT_XOR_BYTEPLANES_V1\n'
def pause_owned_worker_away_from_directory_scan():
 # The frozen worker scans output files after every 400th accepted step.
 # Briefly suspend only our exact worker, away from that scan, so unlink/rename
 # cannot race is_file()/stat(). Physical time, state and integration are untouched.
 while True:
  pids=[]
  for q in Path('/proc').glob('[0-9]*/cmdline'):
   try:args=q.read_bytes().split(b'\0')
   except (FileNotFoundError,PermissionError):continue
   if str(A/'run_protocol.py').encode() in args:pids.append(int(q.parent.name))
  assert len(pids)<=1
  if not pids:return None
  pid=pids[0];status=Path(f'/proc/{pid}/status')
  assert '\nState:\tT' not in status.read_text(),'Worker already stopped by another owner'
  try:os.kill(pid,signal.SIGSTOP)
  except ProcessLookupError:continue
  try:
   while '\nState:\tT' not in status.read_text():time.sleep(.001)
   last=json.loads((O/'full/LAST_ACCEPTED.json').read_text());step=last['step']
   if step%400 and last['time_s']<END-1e-15 and not (O/'full/FAILURE_RECORD.json').exists():return pid
  except BaseException:
   try:os.kill(pid,signal.SIGCONT)
   except ProcessLookupError:pass
   raise
  os.kill(pid,signal.SIGCONT);time.sleep(.05)
def encode(paths,target):
 prev={};records=[]
 with target.open('wb') as out:
  out.write(MAGIC)
  for p in paths:
   data=dict(np.load(p));header=dict(name=p.name,original_sha256=sha(p),original_bytes=p.stat().st_size,arrays=[]);blobs=[]
   for key,a in data.items():
    assert not a.dtype.hasobject
    raw=np.frombuffer(a.tobytes(),np.uint8);base=prev.get(key);xor=base is not None and len(base)==len(raw);delta=raw^base if xor else raw
    blob=zlib.compress(delta.reshape(-1,a.dtype.itemsize).T.copy().tobytes(),6);blobs.append(blob);header['arrays'].append(dict(key=key,dtype=a.dtype.str,shape=list(a.shape),xor_previous_key=xor,packed_bytes=len(blob)));prev[key]=raw.copy()
   payload=json.dumps(header,separators=(',',':')).encode();out.write(struct.pack('<Q',len(payload)));out.write(payload)
   for blob in blobs:out.write(blob)
   records.append(header)
  out.flush();os.fsync(out.fileno())
 return records

def decode(path,restore_folder=None):
 prev={};records=[]
 with path.open('rb') as f:
  assert f.read(len(MAGIC))==MAGIC
  while b:=f.read(8):
   assert len(b)==8;size=struct.unpack('<Q',b)[0];header=json.loads(f.read(size));data={}
   for a in header['arrays']:
    dtype=np.dtype(a['dtype']);blob=f.read(a['packed_bytes']);assert len(blob)==a['packed_bytes'];raw=np.frombuffer(zlib.decompress(blob),np.uint8).reshape(dtype.itemsize,-1).T.copy().ravel()
    if a['xor_previous_key']:raw^=prev[a['key']]
    prev[a['key']]=raw.copy();data[a['key']]=raw.view(dtype).reshape(a['shape'])
   out=io.BytesIO();np.savez_compressed(out,**data);original=out.getvalue();assert len(original)==header['original_bytes'] and hashlib.sha256(original).hexdigest()==header['original_sha256']
   if restore_folder:
    dest=Path(restore_folder);dest.mkdir(exist_ok=True,parents=True);target=dest/header['name'];assert not target.exists();target.write_bytes(original)
   records.append(header)
 return records

if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='restore':
  print('RESTORED',len(decode(Path(sys.argv[2]),Path(sys.argv[3]))));sys.exit(0)
 if len(sys.argv)>1 and sys.argv[1]=='verify':
  print('VERIFIED',len(decode(Path(sys.argv[2]))));sys.exit(0)
 assert not (O/'PAYLOAD_MANIFEST.json').exists(),'Do not mutate a sealed archive; restore/verify modes remain available'
 dest=O/'full/rejected_archives';dest.mkdir(exist_ok=True);temp=Path('/workspace/y_clot/tmp/v3_rejection_archival_'+RUN);temp.mkdir(parents=True,exist_ok=True)
 paths=sorted([p for p in (O/'full/rejected').glob('*.npz') if time.time()-p.stat().st_mtime>45],key=lambda p:(int(p.stem.split('_')[-1]),p.name.startswith('trial')))
 if len(sys.argv)>1:paths=paths[:int(sys.argv[1])]
 for start in range(0,len(paths),200):
  batch=paths[start:start+200];serial=len(list(dest.glob('batch_*.bin')))+1;target=dest/f'batch_{serial:04d}.bin';scratch=temp/target.name;assert not scratch.exists() and not target.exists();tick=time.perf_counter();records=encode(batch,scratch);restored=decode(scratch);assert records==restored
  for p,record in zip(batch,records):assert sha(p)==record['original_sha256']
  archive_sha=sha(scratch)
  audit=dict(status='PASS_BYTE_EXACT_ROUNDTRIP',archive=target.name,archive_sha256=archive_sha,original_bytes=sum(q['original_bytes'] for q in records),archive_bytes=scratch.stat().st_size,files=records,wall_s=time.perf_counter()-tick,solver_runtime_unchanged=True,all_FP64_and_original_NPZ_bytes_recoverable=True,only_closed_rejection_snapshots=True,restore_command='source /workspace/y_clot/env.sh; python '+str(A/'archive_rejections.py')+' restore '+str(target)+' NEW_EMPTY_FOLDER')
  audit['archive_bytes']=scratch.stat().st_size
  pid=pause_owned_worker_away_from_directory_scan()
  try:
   os.replace(scratch,target)
   # Every original byte has been reconstructed and hashed before this deletion.
   for p in batch:p.unlink()
   audit['owned_worker_briefly_suspended_for_file_transaction']=pid is not None
   dump(target.with_suffix('.json'),audit)
  finally:
   if pid is not None:
    try:os.kill(pid,signal.SIGCONT)
    except ProcessLookupError:pass
  print('ARCHIVED',serial,len(records),audit['original_bytes'],audit['archive_bytes'],flush=True)
