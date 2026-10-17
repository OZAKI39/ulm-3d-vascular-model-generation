"""Figures and audit for fixed-model test frames 49/50; no inference or tuning."""
from pathlib import Path
from datetime import datetime, timezone
import sys,json,shutil,difflib,re,argparse
import numpy as np,cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from annotation_store import sha
from vision_evaluation import evaluate_instances,aggregate
from render_supervised_delivery import contours,comparison,rgb

TAG='heldout_49_50_20261016_v1'
OUT=ROOT/'results'/TAG
FIG=ROOT/'figures'/TAG
FREEZE=ROOT/'results/expert_freezes/vision2_heldout_49_50_20261016_v1'
BASE=ROOT/'results/supervised_20261016_v1'
COLORS=['#147eab','#7550a3']


def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False))


def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',facecolor='white',dpi=230)
    plt.close(fig)


def main(figures_only=False):
    FIG.mkdir(exist_ok=figures_only)
    summary=json.loads((OUT/'summary.json').read_text())
    baseline=json.loads((BASE/'summary.json').read_text())
    frozen=json.loads((FREEZE/'frozen.json').read_text())
    rows=json.loads((OUT/'all_rows.json').read_text())
    selected=[r for r in frozen['rows'] if r['record']['split']=='test']
    truth_by_id={};images={};outside={};exclusion=[];recomputed=[]
    for r in selected:
        key=r['record']['id'];images[key]=cv2.imread(str(FREEZE/r['frozen_image_file']))
        with np.load(FREEZE/r['frozen_mask_file']) as z:truth_by_id[key]={k:z[k].copy() for k in z.files}
        full=truth_by_id[key]['SINGLE_BUBBLE'];coverage=np.zeros(full.shape,bool)
        for c in summary['coverage']:
            if c['image_id']!=key:continue
            for y0,x0,y1,x1 in c['boxes_yxyx']:coverage[y0:y1,x0:x1]=True
        outside[key]=[int(n) for n in np.unique(full) if n and not np.any((full==n)&coverage)]
    # Re-score saved predictions independently of model execution; no new inference.
    for r in rows:
        a=np.load(OUT/'predictions'/(r['artifact']+'.npz'))
        full=truth_by_id[r['image_id']]['SINGLE_BUBBLE'];y0,x0,y1,x1=r['box_yxyx']
        counts={int(n):int((full==n).sum()) for n in np.unique(full) if n}
        masks=truth_by_id[r['image_id']]
        check=evaluate_instances(a['truth'],a['labels'],r['objects'],
            (masks['UNCERTAIN'][y0:y1,x0:x1]>0)|(masks['OUT_OF_FOCUS'][y0:y1,x0:x1]>0),counts)
        assert check==r['metrics'],r['artifact']
        recomputed.append(r)
        if r['model']=='opencv' and r['domain']=='ROI':
            for n in check['excluded_gt_ids']:
                yy,xx=np.where(a['truth']==n);cx,cy=float(xx.mean()),float(yy.mean())
                exclusion.append(dict(image_id=r['image_id'],roi_index=r['roi_index'],id=n,
                    full_pixels=counts[n],crop_pixels=len(xx),local_center_xy=[cx,cy],
                    truncated=len(xx)!=counts[n],border_center=not(6<=cx<x1-x0-6 and 6<=cy<y1-y0-6)))
    assert len(exclusion)==13 and sum(len(v) for v in outside.values())==21
    dump(OUT/'coverage_audit.json',dict(excluded_truth=exclusion,unreviewed_saved_label_ids=outside,
                                      recomputed_prediction_rows=len(recomputed),metrics_match=True))
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    xx=np.arange(3)
    for route,name,color,delta in zip(['opencv','cellpose'],['OpenCV','Cellpose pretrained'],COLORS,[-.18,.18]):
        m=summary['metrics'][f'{route}/ROI']
        bars=axes[0,0].bar(xx+delta,[m['center'][k] for k in ['precision','recall','f1']],.36,label=name,color=color)
        axes[0,0].bar_label(bars,fmt='%.3f',fontsize=9,padding=3)
        axes[0,1].bar(np.arange(2)+delta,[m['center']['fp'],m['center']['fn']],.36,label=name,color=color)
        old=baseline['metrics'][('opencv' if route=='opencv' else 'cellpose_pretrained')+'/validation/ROI']
        axes[1,0].plot([0,1],[old['center']['f1'],m['center']['f1']],'o-',c=color,label=name+' center F1')
        axes[1,0].plot([0,1],[old['mask']['f1'],m['mask']['f1']],'s--',c=color,label=name+' mask F1')
        recall=[summary['metrics'][f'{route}/{d}']['center']['recall'] for d in ['ROI','FULL_FRAME']]
        axes[1,1].plot([0,1],recall,'o-',c=color,label=name)
    axes[0,0].set(xticks=xx,xticklabels=['Precision','Recall','Center F1'],ylim=(0,1.17),ylabel='Score',title='(a) Held-out reviewed ROI detection (149 eligible GT)')
    axes[0,0].legend(fontsize=9,loc='upper center',ncol=2)
    axes[0,1].set(xticks=[0,1],xticklabels=['False positives','False negatives'],ylabel='Count',title='(b) Unmatched predictions and missed human labels')
    axes[0,1].legend(fontsize=9)
    axes[1,0].set(xticks=[0,1],xticklabels=['Validation (25/26)','Test (49/50)'],ylim=(0,1.2),ylabel='F1 score',title='(c) Center accuracy does not imply accurate masks')
    axes[1,0].legend(fontsize=8,loc='upper center',ncol=2)
    axes[1,1].set(xticks=[0,1],xticklabels=['ROI inference','Full-frame inference'],ylim=(-.04,1),ylabel='Recall within reviewed test ROIs',title='(d) Full-frame Cellpose inference remains inadequate')
    axes[1,1].legend(fontsize=9)
    fig.suptitle('Fixed-model held-out evaluation | 2 source frames, 5 reviewed ROIs\n162 reviewed objects - 13 border/clipped exclusions = 149 eligible; no test-based tuning',fontsize=13)
    save(fig,'Fig_Heldout_1_fixed_model_accuracy')

    originals=[r for r in rows if r['model']=='opencv' and r['domain']=='ROI']
    fig,axes=plt.subplots(5,3,figsize=(12,16),layout='constrained')
    for i,r in enumerate(originals):
        y0,x0,y1,x1=r['box_yxyx'];crop=images[r['image_id']][y0:y1,x0:x1]
        gt=np.load(OUT/'predictions'/(r['artifact']+'.npz'))['truth']
        axes[i,0].imshow(rgb(crop));contours(axes[i,0],gt,np.unique(gt)[1:],'#9543b3')
        axes[i,0].set_title(f'Frame {int(r["image_id"].split("_")[-1])} / ROI {r["roi_index"]+1} | human masks')
        for j,route in enumerate(['opencv','cellpose'],1):
            rr=next(z for z in rows if z['model']==route and z['domain']=='ROI' and z['image_id']==r['image_id'] and z['roi_index']==r['roi_index'])
            arrays=np.load(OUT/'predictions'/(rr['artifact']+'.npz'))
            comparison(axes[i,j],crop,rr,arrays);m=rr['metrics']['center']
            axes[i,j].set_title(f'{"OpenCV" if j==1 else "Cellpose pretrained"} | TP {m["tp"]} / FP {m["fp"]} / FN {m["fn"]}')
        for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
    fig.legend(handles=[Line2D([],[],c='#9543b3',label='Human mask'),Line2D([],[],c='#05a2c2',label='Center-matched prediction'),
                        Line2D([],[],c='#d64938',label='False positive'),Line2D([],[],c='#e6a800',marker='+',ls='',label='False negative')],
               loc='outside lower center',ncol=4,frameon=False)
    fig.suptitle('Every reviewed test ROI | original image geometry\nNo test training, threshold search or best-example selection',fontsize=13)
    save(fig,'Fig_Heldout_2_all_reviewed_rois')

    fig,axes=plt.subplots(1,2,figsize=(12,5.8),layout='constrained')
    for ax,r,c in zip(axes,selected,summary['coverage']):
        key=r['record']['id'];full=truth_by_id[key]['SINGLE_BUBBLE']
        ax.imshow(rgb(images[key]));all_ids=[int(n) for n in np.unique(full) if n]
        contours(ax,full,[n for n in all_ids if n not in outside[key]],'#9543b3',.6)
        contours(ax,full,outside[key],'#d29314',.8)
        for i,(y0,x0,y1,x1) in enumerate(c['boxes_yxyx'],1):
            ax.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,color='#028275',lw=1.3))
            ax.text(x0+4,y0+18,str(i),color='white',bbox={'facecolor':'#028275','pad':2,'edgecolor':'none'},fontsize=9)
        ax.set_title(f'Image {49 if key.endswith("660") else 50} | {c["reviewed_counts"]["SINGLE_BUBBLE"]} reviewed, {len(outside[key])} outside review')
        ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Coverage audit | only confirmed rectangles supply ground truth\nPurple: masks in reviewed regions; amber: saved labels outside review (not scored)',fontsize=13)
    save(fig,'Fig_Heldout_3_review_coverage')

    if figures_only:
        render=json.loads((OUT/'render_identity.json').read_text())
        render.update(script_sha256=sha(__file__),outputs={str(p.relative_to(ROOT)):sha(p) for p in sorted(FIG.iterdir())})
        dump(OUT/'render_identity.json',render)
        return

    before=OUT/'documents_before';before.mkdir()
    docs=['FINAL_GATE.json','ANNOTATION_READY.md','OPENCV_VS_CELLPOSE.md','DATASET_SPLIT_AND_PROVENANCE.md',
          'SMART_TRACKING_V2_DEVELOPMENT.md','FINAL_RECOMMENDATION.md','GROUP_AND_BUBBLE_ADMISSION.md']
    for name in docs:shutil.copy2(ROOT/name,before/name)
    gate=json.loads((before/'FINAL_GATE.json').read_text())
    gate.update(utc=datetime.now(timezone.utc).isoformat(),status='HELDOUT_ROI_EVALUATION_COMPLETE; LIVE_GROUP_TRACK_ADMISSION_PENDING',
       latest_heldout_result_directory=str(OUT),expert_reviewed_source_frames=6,
       expert_labels_scope='REVIEWED_ROIS_ONLY: 2 train, 2 validation, 2 test source frames; no group/track truth',
       test_reviewed_source_frames=2,test_reviewed_rois=5,test_reviewed_single_masks=162,
       test_eligible_single_masks=149,test_unreviewed_saved_single_masks=21,
       HELDOUT_ROI_EVALUATION_COMPLETE=True,no_test_tuning=True,heldout_accuracy_metrics=summary['metrics'],
       heldout_scope='Within-video time-block test; not blind image exposure or cross-session validation',
       test_frames_now_consumed_for_evaluation=True)
    dump(OUT/'FINAL_GATE.json',gate);dump(ROOT/'FINAL_GATE.json',gate)
    lines=['| 方法 | TP / FP / FN | 精确率 | 召回率 | 中心 F1 | 掩码 F1 | 中心误差均值 / p95 (px) |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for route,name in [('opencv','OpenCV'),('cellpose','Cellpose cyto3 预训练')]:
        m=summary['metrics'][f'{route}/ROI'];c=m['center']
        lines.append(f'| {name} | {c["tp"]} / {c["fp"]} / {c["fn"]} | {c["precision"]:.4f} | {c["recall"]:.4f} | {c["f1"]:.4f} | {m["mask"]["f1"]:.4f} | {m["center_error_px_mean"]:.3f} / {m["center_error_px_p95"]:.3f} |')
    table='\n'.join(lines)
    report=f'''# 第49、50张：冻结参数后的留出区域测试

本次已实际使用这两张新标注完成测试。没有训练、调参、替换模型或连接相机；上一轮选择的OpenCV阈值10、cyto3预训练权重、直径参考11.113245963 px及min_size=5保持不变。原标注和旧结果均保留。

## 数据范围

| 图片序号 / 原始帧 | 已确认区域 | 已审核单泡 | 边缘/裁断排除 | 实际计分 | 另有框外已保存标签 |
|---|---:|---:|---:|---:|---:|
| 49 / 660 | 3 | 97 | 7 | 90 | 0 |
| 50 / 678 | 2 | 65 | 6 | 59 | 21 |
| 合计 | 5 | 162 | 13 | 149 | 21 |

第50张保存了86个标签，但只有65个位于两个有效确认框内；另外21个不是丢失，也不强行改为已审核。本次遵循已保存审核范围，只用已确认像素作真值。13个边缘/裁断排除沿用上一轮规则，逐ID原因见coverage_audit.json；没有为了成绩改动分母。框外不作背景，5个区域只算2张源图。

冻结输入位于`results/expert_freezes/vision2_heldout_49_50_20261016_v1/`，本轮冻结6张有审核区域的图，推理/评分仅取test的660、678两张。训练和验证仍按原时间块分开；自动SAM训练快照实际只包含训练90/108和验证360/374，不包含test标签。测试图由用户基于cyto3候选修正，非双人盲法独立标注；原始660帧曾用于早期无真值工程对照，因此本轮准确表述为**冻结参数的同录像时间块留出真值测试**，不称从未见过图像的完全盲测或跨实验泛化。

## 固定模型结果

{table}

中心以6 px一对一匹配计分，掩码以IoU≥0.5另行匹配。OpenCV在两张测试图的中心F1分别0.868/0.852，比cyto3的0.413/0.539高；本轮OpenCV的漏检及误检都较少。与上一轮验证F1相比，OpenCV由0.841到0.862、cyto3由0.644到0.460，是不同样本的表现差异，**不是本轮训练提升或退化**。

OpenCV仍仅有0.195的mask F1，cyto3为0.377。因此“能找到多数中心”和“准确恢复边界”必须分开：现在不能把OpenCV候选的面积、Rg、形状张量直接当成精确微泡形态。已匹配对象的中心误差也不能代表全部泡；漏检已单列。

完整1476×1120图推理后，仅在相同测试ROI计分：OpenCV TP129/FP20/FN20；cyto3 TP0/FP0/FN149，精确率无定义(null)，召回率0。没有将整帧外部未标区域当负例。cyto3全图输入问题保留为失败；未在test上尝试归一化修补、尺度扫描或阈值挑选。下一次修改只能在训练/验证样本开发，不能反复据这两张test的成绩优化后仍声称是未用过的测试集。

## 工程与交付

本轮推理、加载及核查墙钟共{summary['wall_s']:.3f} s；2张静态图不能给出新的实时FPS/p95。既有离线视频、软件计时和相机未连接状态保持其原报告身份，本轮不重复生成无新增运动信息的视频。没有团、短ID、急转弯或不确定/失焦人工真值，这些仍NOT_EVALUATED。BUBBLE_DETECTION_VALIDATED=false、group_shape=null、training_ready=false，无PZT/AFG硬件控制。

新增三组英文PNG/PDF：

1. `figures/{TAG}/Fig_Heldout_1_fixed_model_accuracy.*`：两模型固定参数测试及原验证对照。
2. `figures/{TAG}/Fig_Heldout_2_all_reviewed_rois.*`：全部5个区域，人工轮廓、匹配、误检和漏检。
3. `figures/{TAG}/Fig_Heldout_3_review_coverage.*`：完整原图及有效审核框，框外21个保存标签单独标示。

数据：preregistered.json在推理前写入，记录模型/config/源码/冻结数据SHA；all_rows.json、predictions/保存逐区域原始预测；summary.json、coverage_audit.json记录合并与排除。重新从20份保存预测独立计分，与运行结果一致。绘图脚本scripts/report_heldout.py，测试入口scripts/evaluate_heldout.py；旧结果不覆盖。

## 下一步

这两张已经成功完成当前单泡留出测试，无需重标。优先在既有训练/验证区域改善误检与边界恢复，并专门核查全图/ROI归一化差异，避免盲目增加训练轮次。短轨迹或团形态是另一个验收目标，只有推进这些目标时才需要相应的少量短ID或团掩码；目前不要求用户继续标整幅密集图。
'''
    (OUT/'REPORT.md').write_text(report)
    sections={
      'ANNOTATION_READY.md':'**第49、50张已完成本轮留出测试，无需重标。**第49张3个确认区、97泡；第50张2个确认区、65泡，另21个已保存标签在确认框外，仍保留。沿用边缘规则后149泡评分。现有单泡离线比较已经完成；后续若开展短轨迹/团形态验收，再补相应少量标签，不要求现在继续标完整幅图。',
      'DATASET_SPLIT_AND_PROVENANCE.md':'test新增第49/50张（源帧660/678），5个确认ROI、162个审核单泡，13个边缘/裁断排除，149个计分。21个框外标签未计分。冻结参数测试已完成；同录像、模型辅助人工修正，且660帧曾参与无真值工程测试，不是完全盲测。此后不能将这两张反复调参仍称未使用test。',
      'OPENCV_VS_CELLPOSE.md':table+'\n\n仅统计已确认测试ROI。OpenCV定位较好但mask F1低；cyto3全图输入在这些ROI漏掉149/149，失败保留。未进行test调参。',
      'SMART_TRACKING_V2_DEVELOPMENT.md':f'已完成第49、50张的冻结参数留出测试，共{summary["wall_s"]:.3f}秒，新增三组PNG/PDF和逐ROI证据。原引擎、评价规则、模型权重与配置均未修改，未训练，未连接相机。',
      'FINAL_RECOMMENDATION.md':'新test仍支持优先完善OpenCV小ROI定位：中心F1 0.862、误检20、漏检21。但mask F1仅0.195，不能视作形状准确；cyto3中心F1 0.460、mask F1 0.377。下一步在train/validation上改进边界与输入一致性，test不用于优化。',
      'GROUP_AND_BUBBLE_ADMISSION.md':'已有train/validation/test各2张源图的局部单泡真值，HELDOUT_ROI_EVALUATION_COMPLETE=true；科学检测、短轨迹、团形态准入仍false。掩码重合度低、场景有限，不能提升为任意实时识别或控制许可。'
    }
    for name,text in sections.items():
        old=(before/name).read_text()
        title,rest=old.split('\n',1)
        notice=f'\n\n> **最新：第49、50张留出测试已完成。** [本次结果](results/{TAG}/REPORT.md)；此前续研记录保留。\n'
        (ROOT/name).write_text(title+notice+rest+f'\n\n## 第49、50张留出测试（2026-10-16）\n\n{text}\n')
    progress=Path('/home/lzy/projects/microbubble_aggregation.md');old=progress.read_text()
    backup=progress.with_name(progress.name+'.bak_F10T_Vision2_test49_50_20261016')
    if backup.exists():raise FileExistsError(backup)
    shutil.copy2(progress,backup)
    block='''**第49、50张的新标注随后完成了冻结参数的留出测试。**第49张三个确认区域含97个单泡，第50张两个确认区域含65个，另21个已保存标签位于确认框外，保留但未计分。沿用原6 px边界/裁断规则排除13个后，149个单泡进入评价；没有使用test训练或调参。OpenCV的TP/FP/FN为128/20/21，精确率0.865、召回率0.859、中心F1 0.862；预训练cyto3为55/35/94、0.611/0.369/0.460。OpenCV仍较适合先做候选定位，但mask IoU≥0.5的F1只有0.195，cyto3为0.377，不能据中心成绩宣称边界、面积或形状准确。全图推理后只在同一审核区评价，cyto3漏掉149/149；失败保留，没有在test上修补参数。

本次实际计算与核查约16.60秒，新增三个英文PNG/PDF图组和所有区域对照，位于smart_online_tracking_v2/results/heldout_49_50_20261016_v1及对应figures目录。它是同录像时间块留出真值测试，人工使用了模型候选辅助，且660原始帧曾用于无真值工程对照；不声称完全盲测或跨实验日泛化。两张test已被评价，后续开发应在train/validation上改进，再保留新独立数据。当前单泡离线比较已完成，无需重标第49、50张；团/短ID仍缺真值，科学与控制准入仍未通过。本轮没有连接相机、没有新增实时性能测量或硬件控制，training_ready=false。

'''
    marker='---\n\n### 尾声：故事里还没讲完的部分';assert old.count(marker)==1
    new=old.replace(marker,block+marker);assert new.replace(block,'')==old
    headings=lambda s:re.findall(r'^#{1,6} .+$',s,re.M)
    assert headings(new)==headings(old)
    progress.write_text(new);p=OUT/'progress_document';p.mkdir()
    (p/'before.md').write_text(old);(p/'after.md').write_text(new)
    (p/'diff.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='before',tofile='after')))
    dump(p/'verification.json',dict(backup=str(backup),before_sha256=sha(backup),after_sha256=sha(progress),
                                  old_text_preserved=True,headings_unchanged=True))
    dump(OUT/'render_identity.json',dict(script_sha256=sha(__file__),inputs={str(p.relative_to(ROOT)):sha(p) for p in
        [OUT/'summary.json',OUT/'all_rows.json',BASE/'summary.json',FREEZE/'frozen.json']},
        outputs={str(p.relative_to(ROOT)):sha(p) for p in sorted(FIG.iterdir())}))
    print('Generated 3 figure pairs, checked 20 prediction records, updated reports and progress without changing labels.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--figures-only',action='store_true')
    main(parser.parse_args().figures_only)
