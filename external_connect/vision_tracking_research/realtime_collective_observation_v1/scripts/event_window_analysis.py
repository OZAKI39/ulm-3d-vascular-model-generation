"""Read-only event comparison; never infers an acoustic effect from simulated labels."""
import argparse,json,gzip,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import dump

def compare(events,states,window_s=2):
 result=[]
 for event in events:
  t=event.get('source_pts');row={'event':event,'window_s':window_s,'acoustic_attribution':'NOT_VERIFIED','scope':'SOFTWARE_SEGMENTATION_ONLY' if event['source']=='SIMULATED_EVENT' else 'OBSERVATIONAL_NO_CAUSAL_CLAIM'}
  if t is None:row.update(status='NOT_OBSERVABLE',reason='NO_SOURCE_PTS_OR_TRUSTED_CLOCK_MAPPING');result.append(row);continue
  for name,left,right in [('before',t-window_s,t),('after',t,t+window_s)]:
   subset=[s for s in states if s.get('session_id')==event.get('session_id') and s.get('source_pts') is not None and left<=s['source_pts']<right and s.get('valid') and s.get('aggregation')]
   metrics={k:float(np.mean([s['aggregation']['regions'][0][k] for s in subset if s['aggregation']['regions'][0][k] is not None])) if any(s['aggregation']['regions'][0][k] is not None for s in subset) else None for k in ['count','nn_px','rg_px','radial_px_s','density_per_10k_px']}
   row[name]={'valid_samples':len(subset),'status':'OBSERVED_WINDOW' if len(subset)>=8 else 'NOT_OBSERVABLE','metrics':metrics}
  row['status']='SEGMENTATION_ONLY';result.append(row)
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('session');p.add_argument('--output',required=True);a=p.parse_args();d=Path(a.session)
 events=json.loads((d/'pzt_events.json').read_text())
 with gzip.open(d/'states.jsonl.gz','rt') as f:states=[json.loads(line) for line in f]
 dump(a.output,compare(events,states))
