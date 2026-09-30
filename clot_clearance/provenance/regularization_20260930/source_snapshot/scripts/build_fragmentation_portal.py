"""Build an offline entry point for the existing prototype's new results."""
import json
from pathlib import Path

root=Path(__file__).resolve().parents[1]
history=json.loads((root/'results/streaming_fragmentation_demo/history.json').read_text())
data=[{key:row[key] for key in ['represented_cycles','broken_bond_fraction','detached_volume_fraction','cleared_volume_fraction','detached_components']} for row in history]
html='''<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>clot_clearance · 断裂与输运验证</title>
<style>
:root{color-scheme:dark}body{margin:0;background:#07101e;color:#e8f0fa;font:17px/1.7 system-ui,sans-serif}main{max-width:1280px;margin:auto;padding:32px 24px 60px}h1{font-size:34px;line-height:1.3}h2{margin:36px 0 12px;font-size:24px}a{color:#66cfff}p{max-width:1000px}.tag{color:#89a9c7;font-size:14px;letter-spacing:.08em}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}.card,section{border:1px solid #26374b;border-radius:12px;padding:18px;background:#0c1728}.card strong{font-size:26px;color:#66cfff;display:block}.muted{color:#aabbd0;font-size:15px}video,img{width:100%;height:auto;border-radius:6px}input{width:100%;accent-color:#66cfff}button{background:#1b3954;color:#fff;border:1px solid #47647b;border-radius:7px;padding:8px 13px;font:inherit;cursor:pointer}nav{display:flex;gap:18px;flex-wrap:wrap}.stats{font-variant-numeric:tabular-nums}ul{padding-left:24px}summary{cursor:pointer}
</style><main>
<div class="tag">EXISTING CLOT_CLEARANCE / STRAIGHT-PIPE WORKFLOW</div>
<h1>从损伤、断键到脱离与输运</h1>
<p>既有 NOSB-PD 原型的独立扩展。载荷与经验循环损伤产生实际断键，脱离材料保留内部活动键并继续接受规定流体力。没有人工开裂命令，所有 360 个粒子均保留，动画使用真实位移。</p>
<p><strong>Qualitative forced-failure verification · UNCALIBRATED</strong><br>仅验证数值流程；不预测临床清除率、真实超声响应或实验碎片尺寸。177 个脱离分量中有 172 个单颗粒，低维支撑使用未经标定的延续近似。</p>
<nav><a href="FRAGMENTATION_REPORT_ZH.md">完整中文报告</a><a href="README_FRAGMENTATION.md">复现说明</a><a href="FRAGMENTATION_MODEL.md">模型与限制</a><a href="OPEN_RESULTS.html">保留的原轻度损伤结果</a></nav>
<div class="cards" style="margin-top:24px">
<div class="card">首次断键<strong>N = 2,000</strong></div><div class="card">首次脱离<strong>N = 11,000</strong></div>
<div class="card">最终脱离体积<strong>60.555556%</strong></div><div class="card">最终清除体积<strong>43.333333%</strong></div>
</div>
<h2>逐状态查看</h2><section>
<img id="frame" src="visualization/streaming_fragmentation_demo/frames/damage/frame_0000.png" alt="实际粒子位置、损伤和稳定碎片编号的双面板图">
<label for="step">代表循环状态</label><input id="step" type="range" min="0" max="25" value="0" step="1">
<div style="display:flex;gap:10px;flex-wrap:wrap"><button type="button" id="play">播放</button><button type="button" data-step="0">初始</button><button type="button" data-step="2">首次断键</button><button type="button" data-step="11">首次脱离</button><button type="button" data-step="25">末态</button></div>
<p class="stats" id="stats" aria-live="polite"></p><p class="muted">滑块读取保存状态，不插值或改变科学数据；首次脱离帧尚未积累明显输运距离。ID 按分裂谱系更新，显式标签显示最大的脱离分量。</p>
</section>
<h2>完整动画</h2><section>
<p>损伤主面板 / 碎片副面板</p><video controls preload="metadata" poster="visualization/streaming_fragmentation_demo/frames/damage/frame_0025.png" src="visualization/streaming_fragmentation_demo/fragmentation_damage.mp4"></video>
<p><a href="visualization/streaming_fragmentation_demo/fragmentation_damage.mp4">下载 MP4</a> · <a href="visualization/streaming_fragmentation_demo/fragmentation_damage.gif">GIF</a></p>
<p>碎片主面板 / 损伤副面板</p><video controls preload="metadata" poster="visualization/streaming_fragmentation_demo/frames/fragments/frame_0025.png" src="visualization/streaming_fragmentation_demo/fragmentation_fragments.mp4"></video>
<p><a href="visualization/streaming_fragmentation_demo/fragmentation_fragments.mp4">下载 MP4</a> · <a href="visualization/streaming_fragmentation_demo/fragmentation_fragments.gif">GIF</a></p>
</section>
<h2>四阶段与时间历程</h2><section>
<a href="visualization/streaming_fragmentation_demo/figures/four_panel_summary.png"><img src="visualization/streaming_fragmentation_demo/figures/four_panel_summary.png" alt="初始、首次断键、首次脱离和后期输运四阶段图"></a>
<a href="visualization/streaming_fragmentation_demo/figures/four_panel_summary.pdf">四阶段 PDF</a>
<a href="visualization/streaming_fragmentation_demo/figures/time_histories.png"><img src="visualization/streaming_fragmentation_demo/figures/time_histories.png" alt="损伤、断键、附着、脱离、清除和碎片数量时间历程"></a>
<a href="visualization/streaming_fragmentation_demo/figures/time_histories.pdf">时间历程 PDF</a>
</section>
<h2>参数、原始数据与证据</h2><section><ul>
<li><a href="configs/streaming_fragmentation_demo.json">完整新配置</a> · <a href="provenance/fragmentation_extension/CONFIG_DIFF.json">相对原温和算例的全部参数变化</a></li>
<li><a href="results/streaming_fragmentation_demo/SUMMARY.json">计算摘要</a> · <a href="results/streaming_fragmentation_demo/history.csv">统计 CSV</a> · <a href="results/streaming_fragmentation_demo/states.npz">完整状态 NPZ</a></li>
<li><a href="results/streaming_fragmentation_demo/particles.pvd">粒子 PVD</a> · <a href="results/streaming_fragmentation_demo/bonds.pvd">键 PVD</a>（本地 ParaView 打开）</li>
<li><a href="results/streaming_fragmentation_demo/components.json">分量与质心</a> · <a href="results/streaming_fragmentation_demo/lineage.json">分裂谱系</a> · <a href="results/streaming_fragmentation_demo/clearance_ledger.json">清除记账</a> · <a href="results/streaming_fragmentation_demo/failure_ledger.jsonl">逐键失效证据</a></li>
<li><a href="results/streaming_fragmentation_demo/INDEPENDENT_AUDIT.json">独立数值核验</a> · <a href="results/streaming_fragmentation_demo/TRANSPORT_FORCE_ATTRIBUTION.json">流体力开关对照</a> · <a href="verification/fragmentation_extension/MILD_REPRODUCTION.json">旧算例逐位回归</a></li>
<li><a href="logs/fragmentation/tests_final.log">26 项测试日志</a> · <a href="provenance/fragmentation_extension/FINAL_QC.json">最终交付核验</a> · <a href="results/streaming_fragmentation_demo/IDENTITY.json">配置与源代码身份</a></li>
</ul><p class="muted">机械代理时间 0.052 s 与 N/1 MHz = 0.025 s 的周期账本不同，均非已标定治疗时间。清除面 x = 1 mm 仅计账，不删除材料。原型仍需真实微泡流场、载荷耦合、材料与损伤参数的科学验证。</p></section>
<script>
const states=__DATA__;
const step=document.getElementById('step'), frame=document.getElementById('frame'), stats=document.getElementById('stats'), play=document.getElementById('play');let timer=null;
function show(){const i=Number(step.value),s=states[i];frame.src='visualization/streaming_fragmentation_demo/frames/damage/frame_'+String(i).padStart(4,'0')+'.png';stats.textContent='N = '+s.represented_cycles.toLocaleString()+' | 断键 '+(100*s.broken_bond_fraction).toFixed(4)+'% | 脱离 '+(100*s.detached_volume_fraction).toFixed(4)+'% | 清除 '+(100*s.cleared_volume_fraction).toFixed(4)+'% | 脱离分量 '+s.detached_components;}
function stop(){clearInterval(timer);timer=null;play.textContent='播放';}
step.addEventListener('input',()=>{stop();show();});
document.querySelectorAll('[data-step]').forEach(b=>b.addEventListener('click',()=>{stop();step.value=b.dataset.step;show();}));
play.addEventListener('click',()=>{if(timer){stop();return;}if(Number(step.value)===25){step.value=0;show();}play.textContent='暂停';timer=setInterval(()=>{step.value=Number(step.value)+1;show();if(Number(step.value)===25)stop();},500);});show();
</script></main></html>
'''
(root/'OPEN_FRAGMENTATION.html').write_text(html.replace('__DATA__',json.dumps(data)),encoding='utf8')
print('Wrote OPEN_FRAGMENTATION.html')
