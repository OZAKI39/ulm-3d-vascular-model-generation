"""WSL-only guarded notebook/entrypoint update after terminal evidence is synced."""
from pathlib import Path
import json,hashlib,datetime,difflib,os
R=Path(__file__).resolve().parents[3];TASK='clot_irregular_lps_discretization_v4';RUN='20261019T221000Z';O=R/'verification'/TASK/RUN;D=R/'analysis'/TASK
def sha(b):return hashlib.sha256(b).hexdigest()
def guarded(path,new,audit):
 before=path.read_bytes();old=before.decode();after=new(old).encode();audit.mkdir(parents=True,exist_ok=False);(audit/'before.txt').write_bytes(before);(audit/'after.txt').write_bytes(after)
 assert path.read_bytes()==before,'Document changed concurrently; do not overwrite'
 tmp=path.with_name(path.name+'.v4.tmp');tmp.write_bytes(after);os.replace(tmp,path)
 (audit/'update.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),after.decode().splitlines(True),fromfile=str(path)+' before',tofile=str(path)+' after')))
 (audit/'HASHES.json').write_text(json.dumps(dict(path=str(path),before_sha256=sha(before),after_sha256=sha(after)),indent=2)+'\n')
def main():
 q=json.loads((O/'NATIVE_V4_RESULT.json').read_text());v=json.loads((O/'DELIVERY_VALIDATION.json').read_text());assert q['status'] in ['COMPLETED','NUMERICAL_FAILURE'] and v['status']=='PASS';s=json.loads((O/'SCIENTIFIC_OUTCOME.json').read_text());e=json.loads((O/'ELASTIC_PATCH_CHECK.json').read_text());geo=json.loads((O/'IRREGULAR_REFERENCE_GEOMETRY.json').read_text());stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');audit=O/'provenance/document_update'/stamp
 base='temp_storage/clot_github_sync_20261015/clot_clearance';report=f'{base}/analysis/{TASK}/REVIEW_CLOT_IRREGULAR_LPS_V4.md';index=f'{base}/results/{TASK}/{RUN}/OPEN_RESULTS.html'
 section=f'''## Irregular Variational LPS Discretization V4

**不规则参考、完整变分求积与单次原生轨迹（{RUN}）**：本轮在8125材料点、1875锚点的标准试样上开发。保留原10强区、每个原始材料ID相区、G/K/ρ、强弱阈值、150.5µm horizon、预缺口物理几何和1000Pa表面剪切协议。按用户确认，守恒V3积分体积1.015625mm³及材料质量1.015625mg；1.2×1.2×0.6mm是材料点中心包络，0.864mm³，两个体积口径分开。

种子20261019，球内各向同性扰动提案，单点模长≤12.5µm；固定外表面、缺口近邻和锚点，5037材料点实际移动。通过微小几何投影保持原强相积分比例33.32923077%，实际扰动分布受边界和投影条件化。最小间距{geo['minimum_point_spacing_m']*1e6:.5f}µm。共享映射控制体由6个辅助四面体积分，体积全部为正；辅助单元只用于几何求积与表面支撑。参考KD-tree重建材料图515058→476895条原始键，端点ID对决定身份；断后不补键。

完整endpoint-star系数改为逐节点Vi、邻点Vj和重新计算的mi；原生力/势统一N/J，质量逐ID取ρVi。接触保留原排斥规律和50µm起始距离，以对称平均体积平方延伸到非等体积；等体积极限核对通过。L1使用真实不规则节点几何；L2使用映射控制体面、全键跨面列表及当前六四面体遮挡，仍为部分近似受力面。合力、力矩、虚功、历史提交、静态试算回退和restart续步均通过代表性验证。

729点仿射小试样体积响应基本不变；三个剪切响应相对变化{100*e['relative_modulus_changes']['xy']:.3f}%/{100*e['relative_modulus_changes']['xz']:.3f}%/{100*e['relative_modulus_changes']['yz']:.3f}%，方向spread约{100*e['shear_directional_spread']['irregular']:.3f}%。这些是相对离散比较；规则试样自身表观剪切约1317Pa而非材料卡940Pa，不能把该小型全场仿射检查当作宏观标定。

V3真实损伤态8步计时发现主要热点为谱系与诊断20.080s、磁盘/事件输出9.079s；原生积分及材料1.150s。原成员ID直接查新标签等保持语义优化，使窗口44.787→28.829s（1.554倍）；状态、全标量CSV和父子事件一致。完整态窗口没有加速，不把短窗提升当作11.7小时全程实测提升。

唯一H_IRREGULAR_LPS_V4状态 **{q['status']}**，合法接受到 **{q['accepted_time_s']*1e3:.8f}ms**，{q['accepted_steps']}接受步、{q['rejected_predictions']}预测拒绝。2ms为数值试运行门槛，后段从同一合法checkpoint续算；BASE=0.25µs、最小事件步0.03125µs，不设4GiB/1800s预算。最大已接受归一化能量残差{100*q['max_relative_R']:.7f}%，R独立、不加入D；wall {q['wall_s']/3600:.4f}h。匹配接受时间的20–200节点游离分量V3/V4={s['final_target_count']['V3']}/{s['final_target_count']['V4']}，V4时程最多{q['max_N20_200']}个。成员由存活原始材料图定义，质量、体积和等效名义点数分别计算。

一次受约束随机点云同时改变参考坐标、连接和积分，不能据此宣称统计稳健性、真实血栓微结构或声流机制已经验证。正式标签EXPLORATORY_UNCALIBRATED、MECHANICAL_SHEAR_SURROGATE；数值停止或未观察窗口均在报告中单独标明，未自动追加物理工况。见[八问科研报告]({report})及[十图、两个真实1080p视频]({index})。视频为固定Z竖直相机、真实接受态和位移×1；球形节点仅表示积分体积，不伪造连续碎块外壳。

'''
 def notebook(old):
  assert '## Irregular Variational LPS Discretization V4' not in old,'V4 notebook section already present'
  marker='# ⑥ 外加剪切';at=old.index(marker);out=old[:at]+section+old[at:]
  marker='# ③ 切成点、连上弹簧\n';note='\n**V4补充**：以下保留原Y血管案例的规则离散说明；最新标准材料试样已采用逐ID非均匀控制体积和不规则内部参考点云，详见下方“ Irregular Variational LPS Discretization V4 ”。\n'
  return out.replace(marker,marker+note,1)
 guarded(Path('/home/lzy/projects/clot_rupture.md'),notebook,audit/'notebook')
 prefix=f'''# 当前：不规则变分LPS参考与单次轨迹 V4（{RUN}）

任务 `{TASK}`；原生状态 **{q['status']}**；合法接受 **{q['accepted_time_s']*1e3:.8f}ms / 3.6ms**。8125材料点+1875锚点，原10强区与阈值保持；非均匀体积及质量、完整变分力、参考KD-tree原始图、映射控制体L1/L2。积分体积1.015625mm³与中心包络0.864mm³分开。

先读 `analysis/{TASK}/FINAL_STATUS.json`、`REVIEW_CLOT_IRREGULAR_LPS_V4.md`、`REPRODUCE.md`。图与视频 `results/{TASK}/{RUN}/OPEN_RESULTS.html`；正式证据 `verification/{TASK}/{RUN}`。独立分支 `clot-clearance/irregular-lps-v4-20261020`。只有1次正式物理轨迹；禁止覆盖冻结输出、自动改阈值或追加案例。EXPLORATORY_UNCALIBRATED、MECHANICAL_SHEAR_SURROGATE。

## 历史入口原文（保留）

'''
 guarded(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',lambda old:prefix+old,audit/'entrypoint');(D/'ENTRYPOINT_PREFIX.txt').write_text(prefix)
 (D/'NOTEBOOK_UPDATE.json').write_text(json.dumps(dict(status='UPDATED_WITH_HASH_GUARD_AND_BACKUP',audit_path=str(audit.relative_to(R))),indent=2)+'\n');print('LOCAL_DOCUMENTS_UPDATED',audit)
if __name__=='__main__':main()
