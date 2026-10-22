"""Minimal progress-document update and portable delivery seal."""
import shutil,difflib,re,subprocess
from .common import *

def main():
    sci=read(BASE/'SCIENTIFIC_RESULTS.json');rows=sci['rows'];ref=max([r for r in rows if r['status']=='PASS' and r['task'].get('action')=='X_PAIR' and not r['task'].get('variant')],key=lambda r:r['completed_s']);cost=read(BASE/'COMPUTE_COST.json');fig=next(r for r in read(BASE/'figures/MANIFEST.json')['figures'] if r['number']==11);vp=BASE/'videos/Video_F10AF_1_stripe_and_pattern_formation.mp4';project=Path('/home/lzy/projects');p=project/'microbubble_aggregation.md';backup=Path(str(p)+'.bak_F10AF_20261019');e=BASE/'evidence';e.mkdir(exist_ok=True)
    heading='### 2.8 可观测尺度的空间排列与慢时间积分（F10AF）'
    if heading not in p.read_text():
        old=p.read_text();assert not backup.exists(),'Existing backup must be inspected, not overwritten';shutil.copy2(p,backup);(e/'progress_before.md').write_text(old)
        assets=project/'microbubble_aggregation_assets';assets.mkdir(exist_ok=True);fp=assets/'f10af_actual_pattern_evolution.png';vv=assets/'f10af_pattern_attempt.mp4';shutil.copy2(BASE/'figures'/fig['files']['png']['path'],fp);shutil.copy2(vp,vv)
        oldroot=str(ROOT).replace('mevo_observable_f10af_20261019','mevo_target_f10ae_20261019');new=old.replace('[F10AE 工程][sim-root]；[指定形态与恢复报告][ae-report]','[F10AF 工程][sim-root]；[可观测尺度研究报告][af-report]');new=new.replace('```text\n'+oldroot+'\n```','```text\n'+str(ROOT)+'\n```',1)
        paragraph='当前分支为 `dev/f10af-observable-multifidelity-20261019`，基于已核查的 F10AE `306d899fa4a9be18b8a926688a6e25d389193138`；每条实验的实际执行提交见对应 `EXECUTION.json`。声学聚集历史依据来自 F10Z，完整数值引擎沿用 F10AA；F10AF 新增受限慢时间积分及大范围空间指标。旧的 `computation_examples/mevo_acoustic_microbubbles` 目录不能直接当作最新基线。'
        new=re.sub(r'^当前分支为.*$',paragraph,new,count=1,flags=re.M);new=re.sub(r'^\[sim-root\]:.*$', '[sim-root]: '+str(ROOT),new,count=1,flags=re.M)
        new=re.sub(r'^\| 本批逐泡参数 \|.*$', '| 本批逐泡参数 | 半径 1.78 µm；32、64、128、200 泡均有有限窗口测试，具体准入按算例区分 |',new,count=1,flags=re.M)
        new=re.sub(r'^\| 当前时间尺度 \|.*$',f'| 当前时间尺度 | F10AF 最长完成 {ref["completed_s"]*1e3:g} ms 完整逐泡参考；128/200 泡本轮限 2 ms，未完成秒级仿真 |',new,count=1,flags=re.M)
        new=new.replace('当前标准步长为 1 µs。','完整参考步长为 1 µs；F10AF 仅在已验证的适用范围内采用自适应慢时间积分，接近接触时返回完整求解器。',1)
        block=f'''{heading}

研究对象从紧凑泡团扩展到约 1.0 × 0.8 mm 的分散微泡云，并复用 F10Y 的二维群体指标。X 双源在 2.5 MHz 下存在约 **304.49 µm** 间距的 X 向稳定位置，但静态周期不等于微泡已形成条带。最长 X 双源试验完成 **{ref['completed_s']*1e3:g} ms**，最大逐泡位移约 **{ref['max_fixed_ID_displacement_m']*1e6:.3f} µm**；尚未确认条带、二维阵列或持续整体旋转，也未证明这些变化能被真实 Chronos 可靠观测。

新增受限 SciPy 慢时间积分，保留原 F1/F2；远离接触时利用原阻力矩阵等于单位阵的性质，接近接触时返回原完整求解。短窗口可明显加速，50 ms 接触增多后优势缩小；紧凑初态在约92 µs触发的原反力平衡失败被如实复现，未放宽门槛。100 ms 的必要精细复算因预测成本超过单条预算而未运行。

本轮 **{sci['HF_attempts']} 次完整物理尝试、{sci['HF_complete']} 条完成**；科学计算阶段墙钟合计约 **{cost['serial_scientific_wall_total_s']/60:.2f} 分钟**，审计、分析和媒体生成另计。SMAC3 与 BoTorch 官方接口已通过冻结实测数据回放验证；由于没有高保真阳性模式，正式多保真筛选未准入。另完成8条有界全物理Sobol对照，不能据此宣称AI方法优于历史CMA-MAE。

![实际空间排列与结构因子](microbubble_aggregation_assets/f10af_actual_pattern_evolution.png)

[28项中文结论与完整报告][af-report] · [实际条带尝试视频](microbubble_aggregation_assets/f10af_pattern_attempt.mp4) · [代码与复现入口][af-code]。12组PNG/PDF和3段60 fps真实数据视频已保存；视频是毫秒物理过程的慢放。声流、自旋、真实PDMS声场及电压—声压映射未校准，`training_ready=false`。下一步优先测量真实声场与群体运动机制，再建立可重复的阳性模式和观测标定。

'''
        assert '\n## 3. 人工标注' in new;new=new.replace('\n## 3. 人工标注','\n'+block+'## 3. 人工标注',1);new+=f'\n[af-report]: {BASE}/F10AF_OBSERVABLE_DEVELOPMENT_REPORT.md\n[af-code]: {BASE}/ACTIVE_CODE_INDEX.md\n'
        before_heads=[s for s in old.splitlines() if re.match(r'^#{1,6} ',s)];after_heads=[s for s in new.splitlines() if re.match(r'^#{1,6} ',s)];assert all(s in after_heads for s in before_heads)
        p.write_text(new);(e/'progress_after.md').write_text(new);(e/'progress.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(backup),tofile=str(p))))
        write(e/'PROGRESS_DOCUMENT_AUDIT.json',dict(backup=str(backup),before_sha256=sha(backup),after_sha256=sha(p),all_old_headings_preserved=True,old_headings=before_heads,new_headings=after_heads,assets={str(z):sha(z) for z in [fp,vv]},identity=identity()))
    # Explicit state-construction provenance supplements the immutable state index.
    write(BASE/'INITIAL_STATE_PROVENANCE.json',dict(wide_seed_rule='462010 + N',independent_seeds=[462101,462102,462103],loose=dict(seed=462030,scale=.2),contact=dict(seed=462040,manual_first_pair='Centers +/-1.4 R on X, gap 0.8 R; remaining positions sampled; not acoustically formed'),prescribed_top=dict(seed=462041,construction='Y scaled by 0.2; Z set to sqrt((0.8 mm - R - 5 um)^2 - Y^2); not buoyancy formed'),compact='Historical F10AE S2_r0 demonstration state',identity=identity()))
    original=Path('/home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1/src/aggregation.py');vendor=ROOT/'scripts/straight_tube_exp1/f10af/vendor/f10y_aggregation.py'
    write(e/'F10Y_PROVENANCE.json',dict(original=str(original.resolve()),original_sha256=sha(original),vendor=str(vendor.relative_to(ROOT)),vendor_sha256=sha(vendor),bytewise_identical=original.read_bytes()==vendor.read_bytes(),units='Numerical input micrometres; no camera pixel calibration',identity=identity()))
    unchanged=['src','scripts/straight_tube_exp1/f10aa','scripts/straight_tube_exp1/f10ab','scripts/straight_tube_exp1/f10ac_qd','scripts/straight_tube_exp1/f10ad_multifrequency','scripts/straight_tube_exp1/f10ae'];diff=subprocess.check_output(['git','diff','306d899fa4a9be18b8a926688a6e25d389193138','HEAD','--',*unchanged],cwd=ROOT,text=True);assert not diff;write(e/'SOURCE_PRESERVATION.json',dict(baseline='306d899fa4a9be18b8a926688a6e25d389193138',unchanged_paths=unchanged,diff_empty=True,identity=identity()))
    raw={r['id']:dict(path=r['trajectory_path'],sha256=r['trajectory_sha256']) for r in read(BASE/'ALL_HF_RUNS.json') if r.get('trajectory_path')}
    for r in read(BASE/'AUDIT.json')['rows']:assert r['raw_sha256']==raw[r['id']]['sha256']
    write(BASE/'REMOTE_RAW_TRAJECTORIES.json',dict(server='/workspace/mevo_acoustic_microbubbles/f10af_observable_20261019',files=raw,scope='Every-step authoritative HF data retained on server; corresponding exact-frame local datasets separately hashed'))
    excluded={'INPUT_OUTPUT_SHA256.json','LOCAL_RECEIPT.json'};files={str(p.relative_to(BASE)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*') if p.is_file() and p.name not in excluded and not p.name.endswith('.tar.gz') and (p.suffix!='.lock' or p.name.endswith('_environment.lock'))}
    source={str(p.relative_to(ROOT)):sha(p) for sub in ['scripts/straight_tube_exp1/f10af','configs/straight_tube_exp1/f10af','tests/f10af'] for p in (ROOT/sub).rglob('*') if p.is_file() and '__pycache__' not in str(p)}
    write(BASE/'INPUT_OUTPUT_SHA256.json',dict(files=files,source_files=source,upstream_protected=read(BASE/'PREREGISTRATION.json')['protected_sha256'],remote_raw='REMOTE_RAW_TRAJECTORIES.json',excluded=['self manifest','local transport receipt','redundant source archives','lock files'],identity=identity()))
if __name__=='__main__':main()
