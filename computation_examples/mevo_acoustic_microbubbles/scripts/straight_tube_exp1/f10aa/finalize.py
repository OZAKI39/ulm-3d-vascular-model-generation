"""Final artifact validation, progress-document preservation and hash manifest."""
import re,difflib,json,datetime
import numpy as np
from PIL import Image,ImageDraw
import fitz,av
from .common import *

def main():
    verify();figs=read(BASE/'figures/MANIFEST.json')['figures'];videos=read(BASE/'videos/MANIFEST.json')['videos']
    assert len(figs)==8 and len(videos)==2
    reviewed=read(BASE/'evidence/VISUAL_REVIEW.json');assert reviewed['status']=='PASS'
    for r in figs:
        p=BASE/'figures'/(r['name']+'.png');assert sha(p)==r['png_sha256'];im=Image.open(p);im.verify()
        p=p.with_suffix('.pdf');assert sha(p)==r['pdf_sha256']
        with fitz.open(p) as pdf:assert len(pdf)==1 and len(pdf[0].get_text())>50
    for r in videos:
        p=BASE/'videos'/r['file'];assert sha(p)==r['sha256'];times=[]
        with av.open(str(p)) as reader:
            for f in reader.decode(video=0):times.append(float(f.pts*f.time_base))
        assert len(times)==r['frames'] and np.allclose(np.diff(times),1/r['fps'],atol=1e-9)
    gate=read(BASE/'FINAL_GATE.json');gate['gates']['FIGURES_VALIDATED']='PASS';gate['gates']['VIDEOS_VALIDATED']='PASS';gate['visual_review_sha256']=sha(BASE/'evidence/VISUAL_REVIEW.json');write(BASE/'FINAL_GATE.json',gate)
    summary=read(BASE/'reports/SUMMARY.json');ledger=read(BASE/'evidence/EXECUTION_LEDGER.json')
    # Preserve original text and heading sequence; insert only a clearly bounded new episode.
    doc=Path('/home/lzy/projects/microbubble_aggregation.md');backup=doc.with_name(doc.name+'.bak_F10AA_20261018')
    evidence=BASE/'evidence/progress_document';evidence.mkdir(exist_ok=True)
    if not backup.exists():backup.write_bytes(doc.read_bytes())
    original=backup.read_text();marker='【F10AA-Simulation，2026-10-18】';anchor='### 尾声：故事里还没讲完的部分'
    text=f'''{marker}这一轮把主要问题转为“在物理和数值门槛不变时能快多少”。先冻结126条历史轨迹及16个合法初态、32个ON/OFF力学状态，再实测原代码瓶颈：128泡反复生成8128×384的稠密接触矩阵最贵；进入接触后，对全部约束列的线性求解进一步放大成本。新增矩阵自由全粒子对约束、向量化润滑阻力与仅按需解活跃列；APECSS、F2、2.5MHz均匀水场、零流、真实球面Sequential QP、壁面和全部原门槛保留。制造的FCC/立方近接触状态仅用于力学核查，最多532个同时活跃约束，不称声学聚集结果。

在相同初态、动作、步长和逐步位置/力/反力输出下，原单线程与新单线程的完整墙钟分别为：32泡2ms，7.787→2.646s（2.944倍）；64泡2ms，34.632→5.760s（6.012倍）；128泡1ms，44.374→6.561s（6.763倍）。128泡在进程启动前固定4个BLAS线程可降至5.293s（8.384倍），运行中切换到4线程曾出现native崩溃，失败日志保留，不能混称所有四线程均不工作。CuPy FP64通过数值门槛但比优化CPU慢；多个独立任务每进程1线程，4进程暖轮总吞吐为单进程2.706倍，不是单条轨迹快2.706倍。服务器有其他CPU任务，不能称整机独占基准。

本轮完整保存59条轨迹尝试，58条完成、1条因预算中止；65,259个接受步独立检查原阻力和最终真实法向的力平衡、间隙互补及原门槛。正式32/64/128窗口的固定ID、COM、Rg、形状张量、长短比、邻近图和事件时序一致，dt/cache通过。优化128泡2ms已完成，Rg缩小3.7435%、硬接触峰13对，半步长固定ID最大差2.5446nm、细缓存0.03285nm；原CPU参考在1.659ms达到240s预算后停止，只能核查公共前缀，不能宣称完整2ms原实现对照通过。256泡仅做合法无接触短步预检，约29.4ms/步，未运行长轨迹或密集状态扩展。

9个Supervisor尝试（含导入/native/API失败）的累计墙钟为{ledger['total_supervisor_wall_s']:.3f}s，不含环境恢复、下载、少量独立诊断及本地绘图/审查。已生成8组PNG/PDF和2段完整解码MP4，固定物理帧间隔10µs、播放30fps，不把6.7s播放时间当物理时间。代码分支dev/f10aa-solver-acceleration-20261018，实际运行提交、原始HDF5、哈希和失败分别保存于新F10AA工程results/straight_tube_exp1/f10aa。后续优先64/128泡的有限更长窗口及多初态/动作，并继续收敛与接触事件检查。加速不证明真实PDMS、电声映射、接触和融合已校准，保留EXPLORATORY、PHYSICS_UNCALIBRATED、HOMOGENEOUS_WATER_CONTROL_PROXY、ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED、training_ready=false。

'''
    assert marker not in original and original.count(anchor)==1
    after=original.replace(anchor,text+anchor)
    current=doc.read_text();assert current==original or marker in current,'Concurrent unrelated document edit; preserve and merge manually'
    heads=lambda s:[v for v in s.splitlines() if re.match(r'^#{1,6}\s',v)]
    assert heads(original)==heads(after)
    before=evidence/'before.md';aft=evidence/'after.md';before.write_text(original);aft.write_text(after)
    (evidence/'change.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),after.splitlines(True),fromfile='before.md',tofile='after.md')))
    doc.write_text(after);write(evidence/'SHA256.json',dict(document=str(doc),backup=str(backup),before_sha256=sha(before),backup_sha256=sha(backup),after_sha256=sha(aft),document_sha256=sha(doc),heading_sequence_unchanged=True,all_original_text_preserved=after.replace(text,'')==original))
    (BASE/'README.md').write_text('# F10AA-Simulation\n\n[中文报告](reports/F10AA_DEVELOPMENT_REPORT.md) · [活动代码](reports/ACTIVE_CODE_INDEX.md) · [推荐配置](reports/FINAL_RECOMMENDATION.md) · [数值Gate](FINAL_GATE.json)\n\n8组科研图：`figures/`。2段MP4：`videos/`。原始HDF5及失败：`full/`、`extension/`、`threaded/`、`parallel/`。每次执行身份和原始日志：`logs/`。完整SHA：`INPUT_OUTPUT_SHA256.json`。\n\n原始大数据本地与远端保留；Git仅保存代码、输入和可审阅的报告/数值证据及小型图表，不包含环境或完整轨迹大文件。详见reports/README.md。\n')
    manifest={}
    for p in BASE.rglob('*'):
        if p.is_file() and p.name not in ['INPUT_OUTPUT_SHA256.json']:
            manifest[str(p.relative_to(BASE))]=dict(sha256=sha(p),bytes=p.stat().st_size)
    write(BASE/'INPUT_OUTPUT_SHA256.json',dict(identity=identity(),files=manifest,source_sha256=identity()['sha256'],scope='F10AA delivery; excludes this self-referential manifest itself, historical external inputs separately frozen in PREREGISTRATION'))
    print('Validated',len(figs),'figures,',len(videos),'videos,',len(manifest),'files; progress document updated')
if __name__=='__main__':main()
