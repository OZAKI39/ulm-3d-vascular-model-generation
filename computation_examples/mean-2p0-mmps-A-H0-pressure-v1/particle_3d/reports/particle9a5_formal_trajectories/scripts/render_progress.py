"""Human-facing figures for each completed development/benchmark stage."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
ROOT=Path(__file__).resolve().parents[4];R=ROOT/'particle_3d/reports/particle9a5_formal_trajectories'
font=R/'assets/NotoSansCJKsc-Regular.otf'
if font.is_file():
 from matplotlib import font_manager
 font_manager.fontManager.addfont(str(font));plt.rcParams['font.family']=FontProperties(fname=font).get_name()
plt.rcParams.update({'figure.facecolor':'#f5f7fb','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'font.size':11})
events=json.loads((R/'data/CORE500_COHORT.json').read_text())['events'];d=np.array([e['diameter_um'] for e in events]);positions=np.array([e['birth_center_m'] for e in events])*1e6
fig,ax=plt.subplots(1,3,figsize=(14,4.5));fig.suptitle('阶段审核：dt = 1.0 ms；CORE500 已冻结，正式轨迹生产尚未完成',fontsize=17)
ax[0].hist(d,bins=24,color='#3676ad',edgecolor='white');ax[0].set(xlabel='Diameter (µm)',ylabel='Accepted births',title='前500个 accepted births；没有按出口挑选')
p=ax[1].scatter(positions[:,0],positions[:,1],c=d,s=14,cmap='viridis',alpha=.75);ax[1].set(xlabel='Inlet x (µm)',ylabel='Inlet y (µm)',title='原始出生位置（投影）');fig.colorbar(p,ax=ax[1],label='Diameter (µm)')
ax[2].plot([0,1,2,3],[1.5,3,6,12],'o-',color='#b7612f');ax[2].set(xticks=[0,1,2,3],xticklabels=['历史','初始','按需延长','按需延长'],ylabel='Trajectory age limit (s)',title='3 → 6 → 12：只延长，不缩短')
fig.tight_layout(rect=(0,0,1,.90))
for ext in ['png','pdf']:fig.savefig(R/'figures'/('Stage_A_Cohort_and_Horizon.'+ext),dpi=300)
plt.close(fig)
(R/'OPEN_RESULTS.html').write_text('''<!doctype html><html lang="zh"><meta charset="utf-8"><title>P9-A.5 正式轨迹审核</title><style>body{font-family:sans-serif;max-width:1250px;margin:32px auto;background:#f5f7fb;color:#203044}img{width:100%}section{padding:24px;background:white;border-radius:12px;margin:20px 0}</style><h1>P9-A.5 正式轨迹：dt = 1.0 ms 验证中</h1><section><h2>已完成</h2><p>CORE500 在轨迹执行前冻结；NEW流场核验通过；旧回归165项通过；历史正式上限核实为1.5秒。</p><p>用户已将名义步长改为 1.0 ms；旧 0.25 ms 结果单独归档，不计入新步长正式样本。当前阶段尚不宣称完成正式500条。下一步比较固定24条轨迹的1/2/4/6/8/12/16 workers，再启动连续扩样。</p></section><section><img src="figures/Stage_A_Cohort_and_Horizon.png"><p><a href="figures/Stage_A_Cohort_and_Horizon.pdf">PDF</a> · <a href="data/CORE500_COHORT.json">CORE500</a> · <a href="data/historical_horizon_audit.json">历史上限核验</a> · <a href="logs/baseline_tests.txt">旧回归日志</a> · <a href="logs/server/server_resource_inventory.txt">服务器资源</a></p></section></html>''')
print('STAGE_A_FIGURES_HTML_READY')
