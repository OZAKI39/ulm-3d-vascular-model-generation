"""Render the completed detector trial without changing selection or labels."""
from pathlib import Path
import argparse, hashlib, json, sys
import cv2
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(out):
    s=json.loads((out/'summary.json').read_text());rows=json.loads((out/'rows.json').read_text())
    frozen=json.loads((out/'freeze/frozen.json').read_text());dataset=json.loads((out/'dataset.json').read_text())
    training=json.loads((out/'training.json').read_text());selection=s['selection']
    figures=out/'figures';figures.mkdir(exist_ok=True)
    arial=Path('/mnt/c/Windows/Fonts/arial.ttf')
    if arial.exists():font_manager.fontManager.addfont(str(arial));plt.rcParams['font.family']='Arial'
    plt.rcParams.update({'font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    names=['opencv_original','opencv_calibrated','cellpose_original','cellpose_tuned']
    titles=['OpenCV original','OpenCV calibrated','Cellpose original','Cellpose fine-tuned']
    colors=['#527b98','#bf854b','#7483b1','#629e85']
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for i,(kind,label) in enumerate([('precision','Precision'),('recall','Recall'),('f1','F1')]):
        vals=[s['metrics'][n+'/validation']['center'][kind] for n in names]
        axes[0].bar(np.arange(4)+(i-1)*.24,vals,.24,label=label,color=['#4b6e91','#71a69b','#c5a26c'][i])
    axes[0].set(xticks=np.arange(4),xticklabels=['OpenCV\noriginal','OpenCV\ncalibrated','Cellpose\noriginal','Cellpose\nfine-tuned'],ylim=(0,1.1),title='Validation: automatic center detection')
    axes[0].legend(frameon=False,ncols=3,loc='upper center',fontsize=9)
    for i,n in enumerate(names):
        c=s['metrics'][n+'/validation']['center']
        axes[1].bar(i-.17,c['fp'],.34,color='#c48a65',label='False positives' if i==0 else None)
        axes[1].bar(i+.17,c['fn'],.34,color='#6e87a4',label='Missed bubbles' if i==0 else None)
    axes[1].set(xticks=np.arange(4),xticklabels=['OpenCV\noriginal','OpenCV\ncalibrated','Cellpose\noriginal','Cellpose\nfine-tuned'],ylabel='Count',title='Validation errors (193 eligible bubbles)')
    axes[1].legend(frameon=False)
    fig.savefig(figures/'detector_validation.png',dpi=200);fig.savefig(figures/'detector_validation.pdf');plt.close(fig)
    ids=['microscope_20261019_train_085s','microscope_20261019_validation_150s']
    fig,axes=plt.subplots(2,5,figsize=(16,8),layout='constrained')
    for ri,fid in enumerate(ids):
        item=next(x for x in frozen['rows'] if x['record']['id']==fid)
        sample=next(r for r in rows if r['image_id']==fid and r['model']=='opencv_original')
        y0,x0,y1,x1=sample['box_yxyx']
        image=cv2.imread(str(out/'freeze'/item['frozen_image_file']))[y0:y1,x0:x1]
        with np.load(out/'freeze'/item['frozen_mask_file']) as z:truth=z['SINGLE_BUBBLE'][y0:y1,x0:x1]
        for ci,name in enumerate(['truth']+names):
            ax=axes[ri,ci];ax.imshow(cv2.cvtColor(image,cv2.COLOR_BGR2RGB))
            if name=='truth':labels=truth;subtitle='Human-reviewed labels';color='#3bb894'
            else:
                r=next(r for r in rows if r['model']==name and r['image_id']==fid)
                with np.load(out/'predictions'/(r['artifact']+'.npz')) as z:labels=z['labels']
                c=r['metrics']['center'];subtitle=f"TP {c['tp']} | FP {c['fp']} | FN {c['fn']}";color='#ecb062'
            for n in np.unique(labels):
                if not n:continue
                contours,_=cv2.findContours(np.uint8(labels==n),cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
                for c in contours:
                    q=c[:,0,:];ax.plot(np.r_[q[:,0],q[0,0]],np.r_[q[:,1],q[0,1]],color=color,lw=.65)
            ax.set_title(('Reference' if name=='truth' else titles[ci-1])+'\n'+subtitle,fontsize=10)
            ax.set(xticks=[],yticks=[])
            if ci==0:ax.set_ylabel('Train: 85 s' if ri==0 else 'Validation: 150 s')
    fig.savefig(figures/'reviewed_roi_comparison.png',dpi=180);fig.savefig(figures/'reviewed_roi_comparison.pdf');plt.close(fig)
    # Candidates remain explicitly rejected; these are evidence, not active profiles.
    search=json.loads((out/'opencv_search.json').read_text())
    ocv=dict(schema='MICROSCOPE_OPENCV_PROFILE_V1',id='new_microscope_20261019_candidate',
        admitted=s['opencv_update_admitted'],camera_image_wh=[1280,1024],parameters=search['configs'][search['selected_index']],
        reason='VALIDATION_RECALL_DECREASED',study=str(out),summary_sha256=sha(out/'summary.json'),
        scientific_scope='Reviewed regions from one compressed recording; not physical calibration')
    cp=dict(schema='MICROSCOPE_CELLPOSE_PROFILE_V1',id='new_microscope_20261019_candidate',
        admitted=s['cellpose_update_admitted'],checkpoint=training.get('weights'),checkpoint_sha256=training.get('weights_sha256'),
        reason='VALIDATION_F1_AND_RECALL_DECREASED',study=str(out),summary_sha256=sha(out/'summary.json'))
    for name,value in [('opencv_candidate_profile.json',ocv),('cellpose_candidate_profile.json',cp)]:
        (out/name).write_text(json.dumps(value,indent=2))
    report='''# 新显微镜检测适配：首轮结果

本轮已完成 Cellpose 微调、OpenCV 小范围参数校准和保留测试评价；两个更新候选均未通过预登记替换条件，默认检测和原始权重保持原样。已实现新权重和独立显微镜检测配置的加载接口。

## 数据与评价范围

- 已确认训练：3 张 / 198 个单泡；验证：5 张 / 203 个；测试：第 19、22 张 / 87 个。
- 第 8 张仍为草稿，排除；只使用已确认区域。排除截断和边缘对象后，训练 197 个、验证 193 个、测试 87 个用于中心检测评价。
- 所有数据来自同一段压缩录像，结果不能称跨实验验证。自动检测不使用人工提示点。
- 中心匹配半径 6 px；轮廓匹配 IoU ≥ 0.5。不同方法共用相同标签和评分规则，大面积错误预测不静默排除。
- 微调只用训练裁片，100 epochs，最终权重一次评价；验证决定选择，选择 JSON 保存之后才评价测试图。测试图不参与优化。

## 验证结果

| 方法 | 精确率 | 召回率 | 中心 F1 | 误检 | 漏检 |
|---|---:|---:|---:|---:|---:|
'''
    for n,title in zip(names,titles):
        c=s['metrics'][n+'/validation']['center'];report+=f"| {title} | {c['precision']:.1%} | {c['recall']:.1%} | {c['f1']:.1%} | {c['fp']} | {c['fn']} |\n"
    report+='''
OpenCV 的候选把阈值从 10 改为 12，误检下降，但漏检增加，不满足本轮“召回率不得下降”的预登记要求。精度和召回率的取舍是真实存在的，不能在看完测试图后放宽门槛以宣布成功。

Cellpose 微调候选在训练集和验证集的召回率均降低。本轮证据只能确定未改善自动候选质量，不能仅凭训练损失下降归因为已适应；少量场景、增强/尺度和优化设置的影响尚未分离。

![验证比较](figures/detector_validation.png)
![实际标签与预测](figures/reviewed_roi_comparison.png)

## 保留测试

候选未通过验证，所以保留测试只评估已选择的原版，不评价后再重新选择。
'''
    for name in ['opencv_original','cellpose_original']:
        c=s['metrics'][name+'/test']['center'];report+=f"\n- {name}: TP {c['tp']}，FP {c['fp']}，FN {c['fn']}；精确率 {c['precision']:.1%}，召回率 {c['recall']:.1%}，F1 {c['f1']:.1%}。"
    report+='\n\n原版 Cellpose 在这两张测试图上漏检很多，不能称已适配。OpenCV 的中心检测表现不等于准确的泡轮廓或真实粒子身份。\n'
    report+=f"\n## 耗时及恢复\n\nCellpose 训练 {training['wall_s']:.2f} s；恢复后的本次计算进程 {s['current_process_wall_s']:.2f} s。异常终止前的完整进程耗时没有可靠记录，未补造总耗时。复用先前保存的 OpenCV 搜索和 17 份预测，未重新训练已完成模型。\n"
    report+=f"\n1280×984 图像上的 OpenCV 检测核 p95 为 {s['benchmark']['opencv_full_image_ms']['p95']:.2f} ms；仅检测核，不包含相机传输、PIV、界面和录像，不能换算成真实相机整套系统帧率。\n"
    report+='''
## 当前可用状态与后续

- SAM 保留前一轮已验证模型；这与自动 Cellpose 检测分别评价。
- 新显微镜工作台保留原版 cyto3。没有把失败候选写入活动模型配置。
- 实时检测接口支持独立配置、相机分辨率检查和逐帧配置身份；无配置时沿用历史行为。当前未启用此次候选。
- 整条管腔未绘制：全通道回放和真实相机精度验证均为 NOT_RUN。启动实时界面后仍需 Draw LIVE mask；没有推断或代画人工确认掩码。
- 下一轮先比较训练区域中的小/淡泡与大/深泡的漏检分布，在训练/验证集上定位 Cellpose 失败原因，再决定有限的优化设置调整。第 19、22 张已作为本轮保留测试使用，后续若反复根据它们修改方案，就不能再称它们为未见测试集。
- 原始标签、失败权重、输入快照、源码、日志和 SHA-256 均保留。
'''
    (out/'REPORT.md').write_text(report)
    manifest={str(p.relative_to(out)):sha(p) for p in [out/'REPORT.md',out/'summary.json',out/'protocol.json',out/'selection.json',out/'training.json',Path(__file__)] if p.is_relative_to(out)}
    manifest.update({str(p.relative_to(out)):sha(p) for p in figures.iterdir() if p.is_file()})
    (out/'report_identity.json').write_text(json.dumps(dict(script=str(Path(__file__).resolve()),script_sha256=sha(__file__),artifacts=manifest),indent=2))
    print(out/'REPORT.md')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();main(a.output.resolve())
