"""Original SAM via micro-sam, isolated from Qt. One request, bounded process lifetime."""
from pathlib import Path
import sys,json,time,traceback
ROOT=Path(__file__).resolve().parents[1]

def main():
    req=json.loads(Path(sys.argv[1]).read_text());start=time.monotonic()
    def status(text):print(json.dumps({'status':text}),flush=True)
    try:
        status('加载原版 SAM ViT-B；界面仍可操作，可随时取消')
        import numpy as np,imageio.v3 as iio,torch
        from micro_sam import util
        torch.set_num_threads(2)
        predictor=util.get_sam_model(model_type='vit_b',device='cuda' if torch.cuda.is_available() else 'cpu',
                                    checkpoint_path=str(ROOT/'models/sam_vit_b_01ec64.pth'))
        image=iio.imread(req['image']);status('读取或计算当前图片的嵌入')
        embeddings=util.precompute_image_embeddings(predictor,image,save_path=req['embedding'],ndim=2,verbose=False)
        util.set_precomputed(predictor,embeddings);status('根据正点 / 负点分割；输出仅为待审核候选')
        positive=np.asarray(req['positive_yx'],float).reshape(-1,2)
        negative=np.asarray(req['negative_yx'],float).reshape(-1,2)
        points=np.concatenate([positive,negative],axis=0)[:,::-1]
        labels=np.r_[np.ones(len(positive),np.int32),np.zeros(len(negative),np.int32)]
        masks,scores,_=predictor.predict(point_coords=points,point_labels=labels,multimask_output=False)
        np.savez_compressed(req['output'],mask=masks[0],score=scores,request_token=req['token'])
        status('DONE');print(json.dumps({'wall_s':time.monotonic()-start,'output':req['output']}),flush=True)
    except Exception:
        traceback.print_exc();sys.exit(1)
if __name__=='__main__':main()
