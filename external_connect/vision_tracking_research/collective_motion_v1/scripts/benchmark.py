import sys,time,json,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
from synthetic import generate
from distribution import CandidateDetector,distribution
from methods import piv,cpd
from affine import parameters,empty

def run(split):
    folder=generate(split);detector=CandidateDetector();rows=[];start=time.perf_counter()
    selection=json.loads((ROOT/'results/selection.json').read_text()) if split=='test' else None
    for path in sorted(folder.glob('*.npz')):
        z=np.load(path);meta=json.loads(path.with_suffix('.json').read_text());ims=z['frames'];H=z['H'];pts=np.arange(len(ims))*meta['dt'];objects=[];dist=[];background=[]
        for k,im in enumerate(ims):
            o,g=detector(im);objects.append(o);dist.append(distribution(o,im.shape)[0]);background.append(np.median(ims[max(0,k-5):k+1],axis=0))
        for end in P['evaluation_end_frames']:
            for window in P['windows']:
                begin=end-window+1;dt=pts[end]-pts[begin];true=None
                if meta['global_affine_truth_available']:
                    HH=H[end]@np.linalg.inv(H[begin]);true=parameters(HH[:2,:2],HH[:2,2],dt,np.array([128.,128.])).tolist()
                modes=['raw','local','fixed_suppressed','candidate','multipass'] if split=='development' else [selection['mode']]
                results={}
                for mode in modes:
                    results['PIV_'+mode]=piv(ims[begin],ims[end],dt,mode='local' if mode=='multipass' else mode,objects_a=objects[begin],objects_b=objects[end],background_a=background[begin],background_b=background[end],multipass=mode=='multipass',bootstrap=P['bootstrap_replicates'])
                results['CPD']=cpd([o['xy'] for o in objects[begin]],[o['xy'] for o in objects[end]],dt,ims[0].shape)
                results['statistics']={'estimate':empty('MEAN_AND_COVARIANCE_DO_NOT_IDENTIFY_FULL_MATERIAL_AFFINE_MOTION','point_statistics'),'algorithm_ms':0.0,
                                       'distribution_before':dist[begin],'distribution_after':dist[end],
                                       'mean_displacement_px':(np.array(dist[end]['group_center_xy'])-np.array(dist[begin]['group_center_xy'])).tolist() if objects[begin] and objects[end] else None}
                for method,result in results.items():
                    result.pop('vectors',None)
                    row={'split':split,'case':meta['case'],'seed':meta['seed'],'method':method,'begin':begin,'end':end,'window':window,'dt_s':dt,'truth':true,'n_detected':[len(objects[begin]),len(objects[end])],**result}
                    rows.append(row)
        dump(ROOT/f'results/benchmark_{split}.json',{'identity':identity(),'rows':rows,'elapsed_s':time.perf_counter()-start,'complete':False})
        print(split,path.stem,'rows',len(rows),'elapsed',round(time.perf_counter()-start,1),flush=True)
        if time.perf_counter()-start>1800:raise RuntimeError('Preregistered batch walltime limit')
    dump(ROOT/f'results/benchmark_{split}.json',{'identity':identity(),'rows':rows,'elapsed_s':time.perf_counter()-start,'complete':True})
    if split=='development':
        scores={}
        for m in ['raw','local','fixed_suppressed','candidate','multipass']:
            rr=[r for r in rows if r['method']=='PIV_'+m and r['truth'] is not None];valid=[r for r in rr if all(r['estimate']['valid'].values())]
            err=[np.mean(np.abs(np.array(list(r['estimate']['values'].values()))-r['truth'])/P['parameter_error_limits']) for r in valid]
            scores[m]={'valid_fraction':len(valid)/len(rr),'mean_normalized_error':np.mean(err) if err else 1e6,'penalized_score':(np.mean(np.minimum(err,10)) if err else 10)+2*(1-len(valid)/len(rr))}
        best=min(scores,key=lambda m:scores[m]['penalized_score'])
        dump(ROOT/'results/selection.json',{'mode':best,'rule':'lowest development clipped normalized error + 2*(1-valid_fraction); independent frozen test never tunes','scores':scores,'identity':identity()})
        print('selection',best,scores,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('split',choices=['development','test']);run(p.parse_args().split)
