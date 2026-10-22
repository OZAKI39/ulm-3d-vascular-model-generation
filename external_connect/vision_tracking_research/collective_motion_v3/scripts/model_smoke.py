"""Official checkpoint load/inference and known-shift smoke test, no fitting."""
import sys,time,argparse,traceback
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from motion_models import *
import av
def main():
    p=argparse.ArgumentParser();p.add_argument('model',choices=['LRAFT_ST','CoTracker3','UnLiteFlowNet']);a=p.parse_args();record=dict(identity=identity(),name=a.model)
    try:
        source=old.read_frames([10],VIDEO)[10];frames,meta=old.original_synthetic.build_sequence('translation',1100000,source,length=16);det=old.CandidateDetector();objects=[det(f)[0] for f in frames];times=np.arange(16)*meta['dt']
        model=CoTrackerOfficial() if a.model=='CoTracker3' else DenseOfficial(a.model)
        r=model.pair(frames,objects,times);r.pop('field',None)
        record.update(status='PASS_FINITE_INFERENCE' if r['estimate'].get('raw_estimate') is not None else 'RUN_WITH_REJECTION',result=r,weight_sha256=sha(model.path),truth=[18.,-12.,0.,0.,0.,0.])
    except Exception as e:record.update(status='FAILED',error=repr(e),traceback=traceback.format_exc())
    dump(OUT/'evidence'/f'model_smoke_{a.model}.json',record);print(record['status'],record.get('error'),flush=True)
if __name__=='__main__':main()
