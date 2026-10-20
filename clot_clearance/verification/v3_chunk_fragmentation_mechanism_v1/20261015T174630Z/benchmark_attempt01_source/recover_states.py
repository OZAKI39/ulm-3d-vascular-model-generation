"""Selectively decode ORIGINAL accepted-before states; never use rejected candidates.

Reuses the sealed XOR/byte-plane archive format, but only materializes the
requested records. Archive SHA covers all bytes; selected NPZ SHA must match
the original byte-exact hash. This adds selective retrieval, not new physics.
"""
from common import *
import io,struct,zlib,time
MAGIC=b'V3_REJECT_XOR_BYTEPLANES_V1\n'
def main():
 start=time.perf_counter();audit=json.loads((B/'full/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());trials=audit['trials']
 # Fixed recorded-event neighborhoods, not states chosen for desired outcomes.
 targets={1941,1943,1945,3280,3282};chosen=[x for x in trials if x['decision']=='REJECT' and x['from_step'] in targets]
 wanted={f'before_{x["trial"]:06d}.npz':x for x in chosen};dest=O/'derived/recovered_accepted';dest.mkdir(exist_ok=True);records=[]
 for meta in sorted((B/'full/rejected_archives').glob('batch_*.json')):
  md=json.loads(meta.read_text());matches={x['name'] for x in md['files']} & set(wanted)
  if not matches:continue
  path=meta.parent/md['archive'];assert sha(path)==md['archive_sha256'];prev={}
  with path.open('rb') as f:
   assert f.read(len(MAGIC))==MAGIC
   while length:=f.read(8):
    header=json.loads(f.read(struct.unpack('<Q',length)[0]));data={}
    for a in header['arrays']:
     dtype=np.dtype(a['dtype']);blob=f.read(a['packed_bytes']);assert len(blob)==a['packed_bytes']
     raw=np.frombuffer(zlib.decompress(blob),np.uint8).reshape(dtype.itemsize,-1).T.copy().ravel()
     if a['xor_previous_key']:raw^=prev[a['key']]
     prev[a['key']]=raw.copy()
     if header['name'] in matches:data[a['key']]=raw.view(dtype).reshape(a['shape'])
    if header['name'] in matches:
     buf=io.BytesIO();np.savez_compressed(buf,**data);original=buf.getvalue();assert len(original)==header['original_bytes'] and hashlib.sha256(original).hexdigest()==header['original_sha256']
     trial=wanted[header['name']];assert int(data['step'])==trial['from_step'] and float(data['time'])==trial['start_s']
     p=dest/header['name'];assert not p.exists();p.write_bytes(original)
     records.append(dict(file=p,original_name=header['name'],step=int(data['step']),time_s=float(data['time']),original_sha256=header['original_sha256'],archive=path,archive_sha256=md['archive_sha256'],semantics='ACCEPTED state before a later rejected candidate; not the rejected trial state'))
     print('RESTORED_ACCEPTED',header['name'],int(data['step']),float(data['time']),flush=True)
 assert {r['step'] for r in records}==targets
 dump(O/'RECOVERED_ACCEPTED_STATES.json',dict(status='PASS_BYTE_EXACT',records=records,wall_s=time.perf_counter()-start))
 p=O/'BASELINE_AUDIT.json';z=json.loads(p.read_text());z['event_1943_exact_node_frame_present']=True;z['event_1943_source']='Original accepted-before snapshot restored from rejection archive; see RECOVERED_ACCEPTED_STATES.json';dump(p,z)
if __name__=='__main__':main()
