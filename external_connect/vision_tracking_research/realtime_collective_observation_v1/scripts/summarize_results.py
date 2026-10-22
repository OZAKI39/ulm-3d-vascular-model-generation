import json,gzip,sys,collections
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,stats

def rows(path):
 with gzip.open(path,'rt') as f:return [json.loads(line) for line in f]
def main():
 study=ROOT/'results/formal_study';r=json.loads((study/'replay_1x/summary.json').read_text());data=rows(study/'replay_1x/states.jsonl.gz');offline=rows(study/'offline.jsonl.gz');synth=json.loads((ROOT/'results/synthetic_frozen/summary.json').read_text());lab=next((ROOT/'results/local_gui_formal').glob('session_*/summary.json'),None);local=json.loads(lab.read_text()) if lab else None
 span=r['wall_s'];ttl100=sum(v['valid'] and v['timings']['end_to_end_software_ms']<=100 for v in data)/span
 summary=dict(replay=r,local_gui=local,matched_100ms_fresh_motion_hz=ttl100,offline=dict(frames=len(offline),median_candidates=float(np.median([len(v['candidates']) for v in offline])),mean_valid_windows=float(np.mean([v['valid_windows'] for v in offline])),aggregation_statuses=dict(collections.Counter(v['aggregation']['status'] for v in offline))),synthetic=synth)
 dump(ROOT/'results/delivery_metrics.json',summary)
 table='| 阶段 | p50 ms | p95 ms |\n|---|---:|---:|\n'+''.join(f"| {k} | {r['stages'][k]['p50']:.3f} | {r['stages'][k]['p95']:.3f} |\n" for k in ['decode_ms','scan_ms','piv_ms','aggregation_ms','analysis_ms','complete_processing_ms','queue_age_ms','end_to_end_software_ms','publish_ms'])
 resources=r['resources'];cpu=np.mean([x['cpu_percent'] for x in resources]);gpu=np.mean([x['gpu_percent'] for x in resources if x.get('gpu_percent') is not None]);mem=max(x['rss_mib'] for x in resources)
 (ROOT/'REALTIME_PERFORMANCE.md').write_text(f'''# 实测性能

正式执行快照 `{r['identity']['git_head']}`，远程 RTX4090、原 CuPy B4。下表是完整937帧按真实PTS释放，分析允许最新帧覆盖的真实运行数据。

{table}

接收 {r['source_sent']} 帧，分析 {r['counters']['processed']} 帧，分析跳过 {r['analysis_skipped']}，过期入口拒绝 {r['counters'].get('analysis_expired',0)}，新鲜有效运动 {r['counters']['fresh_valid']} 次，{r['fresh_valid_hz']:.3f} Hz；总墙钟 {span:.3f} s。100 ms 同门槛运动更新为 {ttl100:.3f} Hz；它是每帧区域运动状态，不是旧F10X的带多个六参数组状态，不能直接当作算法精度提高。

decode+analysis 是完整计算时间，release→completion 是软件接收后年龄，包含排队；UI呈现另测。两者都不是相机曝光到显示的延迟。目标33.3ms的p95尚未达到。不能声称逐帧30fps或真实60fps分析。

CPU 平均 {cpu:.1f}%（100%=一个逻辑核），GPU设备平均 {gpu:.2f}%，RSS峰值 {mem:.1f} MiB。低GPU占用不是任务失败；全通道检测、聚集指标、Python对象和显示都有CPU成本，不为刷利用率增加负荷。

优化证据：pilot ui_v2 scan/PIV/aggregation p95 20.91/29.19/17.08ms，analysis p95 87.18ms。配置OpenCV/解码线程、独立TTL与资源采样、紧凑异步记录后，正式各阶段如表。运行负载含UI/编码与headless不同，故不把pilot与正式比值视为严格单因素加速倍数。

UI、1×/2×/0.5×、长稳态有独立 summary；最终报告同时列出这些差异。profile滚动窗口上限6000，长稳态完整日志保留，可重算。
''')
 gui_rows=[]
 for folder,label in [('gpu_ui_final','远程 RTX4090 / CuPy / GUI + 录屏'),('local_gui_formal','本机 CPU / Rust / WSLg'),('local_gpu_formal','本机 RTX4060 / CuPy / 冷缓存'),('local_gpu_warm','本机 RTX4060 / CuPy / 热缓存')]:
  p=next((ROOT/'results'/folder).glob('session_*/summary.json'))
  measured=json.loads(p.read_text())
  stage=measured['stages']['complete_processing_ms']
  gui_rows.append(f"| {label} | {measured['fresh_valid_hz']:.2f} | {measured['ui_refresh_hz']:.2f} | {stage['p50']:.2f} / {stage['p95']:.2f} | {measured['analysis_skipped']} |\n")
 with (ROOT/'REALTIME_PERFORMANCE.md').open('a') as f:
  f.write('\n## 实际界面与本机后端\n\n每次均完整接收937帧；新鲜有效更新和UI刷新是不同指标。各运行保留自己的代码哈希、配置、资源与summary，不能将下列运行当作严格单因素基准。\n\n| 运行环境 | 新鲜有效 Hz | UI刷新 Hz | 计算 p50 / p95 ms | 跳过分析帧 |\n|---|---:|---:|---:|---:|\n'+''.join(gui_rows)+'\n首次CuPy启动包含CUDA编译，冷缓存运行的单次计算最大约7.38s；旧状态会失效，接收与最新原图显示仍独立继续。本机启动脚本默认使用已安装的GPU后端，可通过 `--backend B2` 显式使用Rust。最终远程界面的Qt画布paint p50/p95为4.78/10.22ms；这是绘制调用耗时，不是曝光到屏幕的光学延迟。\n')
 (ROOT/'REPLAY_SIMULATION_REPORT.md').write_text('# 真实录像回放\n\n源SHA `3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c`；1476×1120，937帧，PTS 0–31.1688s，步长0.0333s，呈现时长31.2021s，容器约31.232s。完整核查见 results/software_validation/source_audit.json。没有修改原视频。\n\n'+''.join(f"- {n}: {json.loads((study/n/'summary.json').read_text())['counters']}；跳过 {json.loads((study/n/'summary.json').read_text())['analysis_skipped']}。\n" for n in ['replay_1x','stress_2x','debug_half'])+'\n1×真实节奏；2×只是约60.06Hz软件加速输入；0.5×约15.02Hz调试。循环会话重建、暂停/恢复、EOF与停止已有实际测试，不复制帧冒充新运动。\n')
 sr=synth['rows'];table='| case / seed | math | measured image | first image convergence s |\n|---|---|---|---|\n'+''.join(f"| {v['case']} / {v['seed']} | {v['math_status']} | {v['image_status']} | {v['image_response_s']} |\n" for v in sr)
 (ROOT/'SYNTHETIC_AGGREGATION_VALIDATION.md').write_text(f'''# SYNTHETIC_AGGREGATION_TEST

开发种子94001用于工程调试。最终冻结种子95001/95002在运行前登记，阈值未按最终结果调整。10类场景×2种子；先用已知点集/速度验证数学趋势，再用真实录像图块、冻结OpenCV检测器、原CuPy PIV估计图像趋势。没有PZT声学模拟、没有训练大模型。

数学方向/拒判正确 {synth['math_correct']}/{len(sr)}；图像方向/拒判正确 {synth['image_correct']}/{len(sr)}；图像误报 {synth['image_false_positive']}；图像拒判 {synth['image_rejected']}。这里将预期不可观测的漏检/失焦拒判算作正确防护，不能将这个比例称为真实聚集精度或独立跨视频泛化率。素材仍来自同一条真实录像，数学种子独立，光学域不独立。

{table}

完整每帧记录、输入图块SHA、clips与失败均保留。Onset出现后需时间窗积累；图像处理和数学窗口的首次提示时间分别记录。真实录像聚集准确率 NOT_VERIFIED。
''')
 (ROOT/'CAMERA_FAILURE_SIMULATION.md').write_text('''# 故障模拟与恢复

results/software_validation/fault_events.json 保存无新帧、短/长断流、重新连接、重复ID、PTS跳变、错尺寸、空帧、抖动/源丢帧、事件跨断流、在途停止、坏二进制的实际事件与观察结果。unit tests另含慢处理拥塞、source切换、灰度/BGR、会话重置和Qt意外关闭。

results/bridge_mock/summary.json 用真实OS管道测试继承的WindowsLatestBridge循环、多个有限子进程会话与EOF回收；results/windows_ipc_mock/summary.json 则是实际Windows Python→WSL的1280×1024合成BGR传输。后者没有相机，也不是公网带宽测量。

故障测试不表示真实Chronos发生过相同错误。重复帧/坏尺寸拒绝，长无新帧TTL失效；恢复需要合法新帧，跨session不衔接运动。stop/UI close都回收线程和子进程。真实相机端到端/重连/曝光时钟仍NOT_MEASURED。

已保留开发失败：Qt event命名冲突、录屏初版编码PTS冲突、SciPy整通道速度不足；相应修复有单独提交与测试。v4/v5任务在等待旧稳定性任务时被本任务主动停止，用最终v6快照替代，不算正式科研运行。
''')
if __name__=='__main__':main()
