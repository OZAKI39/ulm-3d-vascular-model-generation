"""Run locally after final server evidence is copied and verified."""
from pathlib import Path
import json,hashlib,difflib,os,shutil
from datetime import datetime,timezone
ROOT=Path('/home/lzy/projects/temp_storage/clot_github_sync_20261015/clot_clearance')
TASK='v3_chunk_fragmentation_load_direction_control_v1';RUN='20261016T105642Z';O=ROOT/'verification'/TASK/RUN
def read(name):return json.loads((O/name).read_text())
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 status=read('FINAL_STATUS.json');assert status['NUMERICAL_STATUS']=='FULL_PROTOCOL_NUMERICAL_PASS';assert not (O/'WSL_DOCUMENT_UPDATE.json').exists()
 p=Path('/home/lzy/projects/clot_rupture.md');before=p.read_bytes();assert sha(before)=='2d6de693154a3b83a7951a46e9fc402c8813ca415fc6ea72ac9f50d3784c4d9a','Document changed; reconcile before insertion'
 text=before.decode();marker='# ⑥ 外加剪切';assert text.count(marker)==1;assert 'V3参考外法向载荷对照（20261016T105642Z' not in text
 top=read('V3_DIRECTION_CONTROL_TOPOLOGY.json')['reference_normal'];fin=top['final'];maximum=top['history_maxima'];cand=read('V3_CHUNK_CANDIDATES.json');energy=read('DYNAMIC_ENERGY_LEDGER.json');ne=energy['runs']['reference_normal']['final'];oe=energy['runs']['historical_tangent']['final'];early=read('EARLY_WINDOW_COMPARISON.json')['comparison'];cost=read('RUN_RESOURCES.json');audit=read('NUMERICAL_AUDIT.json');interpret=read('SCIENTIFIC_INTERPRETATION.json')
 block=f'''\n**V3参考外法向载荷对照（{RUN}，2026-10-16）**：完成唯一一组真实V3法向合成牵引工况。453个原始加载ID、面积权重与原始1927材料点/1561固定锚点、157637连接均保持。法向来自原V3真实暴露表面4918个有向三角形，依原面片最近加载点归属进行面积加权、归一化；只用暴露面，排除附着面。453点全部通过“向外1µm位于血栓外且血管腔内、向内1µm位于血栓内”的闭合几何检查，无待定方向。原面积逐点重建差为0；曲面法向与历史切向约成90°，保留历史切向模长因子。

历史Gaussian中心、scale=780.2322247780305Pa、sigma=0.28mm、1.1605ms脉冲和当前位置计算规则完全复用；只替换参考方向向量，同坐标同时间逐点力模长核对通过。没有新增加载点、重复叠加切向力或提高幅值。原C++变分软化库、完整endpoint-star能量、本体/界面材料、原始m、事件分辨VV/拒绝回滚、接触、壁碰撞、背景COM阻力和附着数值阻尼均未改，未重编译，历史C1没有重跑。

0.8ms实际接受态：法向本体/界面完全失效{int(early['bulk_failed']['normal'])}/{int(early['interface_failed']['normal'])}条、脱离{int(early['detached_nodes']['normal'])}点、最大自由分量{int(early['largest_detached_nodes']['normal'])}点；相同物理时间附近的历史切向接受态分别为{int(early['bulk_failed']['historical_tangent'])}/{int(early['interface_failed']['historical_tangent'])}条、{int(early['detached_nodes']['historical_tangent'])}点、{int(early['largest_detached_nodes']['historical_tangent'])}点。保存永久原生检查点，实际重新载入并核验x/v/力/κ/d/拓扑/物理标量，保留累计外功和余额，从同一接受态续算至2.33ms；没有重算前缀。

法向全程接受{status['accepted_steps']}步，末态本体/界面断键{audit['last_accepted']['bulk_failed']}/{audit['last_accepted']['interface_failed']}条，实际脱离{fin['detached_nodes']}点、{fin['detached_components']}分量，质量{fin['detached_mass_kg']:.9g}kg；历史最大/最终最大自由分量{maximum['largest_detached_nodes']}/{fin['largest_detached_nodes']}点。≥14/34/114节点历史最大同时分量数{maximum['N_ge14']}/{maximum['N_ge34']}/{maximum['N_ge114']}，末态1–2节点未解析碎屑{fin['tiny_nodes']}点。上一轮原定义筛查得到≥14/34/114节点候选唯一ID集合数{cand['unique_identity_sets']['14']}/{cand['unique_identity_sets']['34']}/{cand['unique_identity_sets']['114']}；诊断候选不等于真实脱落，筛查非穷尽。{interpret['COHESION_COMPARISON']} {interpret['SECONDARY_BREAKUP']}

{interpret['EARLIEST_DAMAGE_COMPARISON']} 失效空间组织按原参考键中点、原自由表面距离和原加载区分区量化；没有将损伤云绘制为连续物理裂面。

法向人工加载功{ne['W_load_J']*1e9:.9g}nJ，历史切向{oe['W_load_J']*1e9:.9g}nJ，比值{energy['normal_over_tangent_load_work']:.9g}；法向末态D={ne['D_J']*1e9:.9g}nJ、K={ne['K_J']*1e9:.9g}nJ。不同轨迹改变Gaussian取值和力—速度关系，因此相同标量规则不是相同输入功。本轮未执行等能量或第二幅值对照，不能认证法向加载的能量效率更高。

数值状态`{status['NUMERICAL_STATUS']}`；最大归一化独立能量余额{100*audit['maxima']['relative_R']['value']:.8g}%，R未计为材料D。原生/独立Python力、能量、不可逆历史、原始身份、锚点与VV核验通过。单rank单线程，主段墙钟{cost['main_wall_s']:.3f}s、CPU {cost['main_CPU_s']:.3f}s、峰值RSS {cost['peak_RSS_bytes']/1024**3:.6f}GiB，拒绝{audit['rejected_damage_trials']}次；未新增墙钟/RSS停止门槛。节点图使用真实接受态、位移×1、共同尺度/相机；历史采样不一致处标注实际匹配时间，不插值制造断键。

`DIRECTION_CONTROL_STATUS={status['DIRECTION_CONTROL_STATUS']}`；`CHUNK_FORMATION_STATUS={status['CHUNK_FORMATION_STATUS']}`（≥34节点/≥200µm参考体积尺度）；`COHESION_RETENTION_STATUS={status['COHESION_RETENTION_STATUS']}`；`ENERGY_COMPARABILITY_STATUS={status['ENERGY_COMPARABILITY_STATUS']}`。{interpret['SCIENTIFIC_CONCLUSION']} 本体与界面Gc仍未标定，法向牵引为合成力学对照，不是已确认的微泡力；200/300µm并非医学安全目标或已认证空间分辨率。下一轮建议：{interpret['NEXT_RECOMMENDATION']} 本轮完成后停止。

见[载荷方向对照图表](temp_storage/clot_github_sync_20261015/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)、[28项科学审查](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{TASK}/{RUN}/REVIEW_V3_CHUNK_FRAGMENTATION_LOAD_DIRECTION_CONTROL_V1.md)。\n\n'''
 index=text.index(marker);after=(text[:index]+block+text[index:]).encode();assert after.decode().replace(block,'',1)==text
 out=O/'wsl_document';out.mkdir(exist_ok=True);(out/'clot_rupture.before.md').write_bytes(before);(out/'clot_rupture.after.md').write_bytes(after);(out/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),after.decode().splitlines(True),fromfile=str(p)+'.before',tofile=str(p)+'.after')))
 assert sha(p.read_bytes())==sha(before);tmp=p.with_name(p.name+'.direction.tmp');shutil.copy2(p,tmp)
 with tmp.open('wb') as f:f.write(after);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p);assert p.read_bytes()==after
 record=dict(status='PASS',actual_WSL_path=str(p),utc=datetime.now(timezone.utc).isoformat(),before_sha256=sha(before),after_sha256=sha(after),readback_sha256=sha(p.read_bytes()),all_historical_text_preserved=True,insert_start_line=text[:index].count('\n')+1,backup='wsl_document/clot_rupture.before.md',diff='wsl_document/clot_rupture.diff')
 (O/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n');print(json.dumps(record,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
