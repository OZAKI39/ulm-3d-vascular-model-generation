"""Update the actual WSL scientific record by insertion; preserve prior bytes."""
from pathlib import Path
import json,hashlib,difflib,os,shutil
from datetime import datetime,timezone
root=Path('/home/lzy/projects/temp_storage/clot_github_sync_20261015/clot_clearance')
out=root/'verification/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z'
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 p=Path('/home/lzy/projects/clot_rupture.md');status=json.loads((out/'FINAL_STATUS.json').read_text());small=json.loads((out/'SMALL_SPECIMEN_VALIDATION.json').read_text());assert status['status']=='COMPLETED_P0_P4_NO_NEW_V3'
 assert not (out/'WSL_DOCUMENT_UPDATE.json').exists(),'Already recorded; do not duplicate insertion'
 before=p.read_bytes();text=before.decode();assert sha(before)=='4a6b49d5d5e8a3908643bef1873ed7d7b7f102d12a84feb24c36aeacbc0fe6e7','Document changed: reconcile before editing'
 marker='# ⑥ 外加剪切';assert text.count(marker)==1;assert 'V3块体形成机制诊断（20261015T174630Z' not in text
 c0,c1=small['cases'];insert=f'''\n**V3块体形成机制诊断（20261015T174630Z，2026-10-15）**：完成P0代码依赖审计、1414个冻结基线文件哈希核对、全部10503个含初态接受拓扑复核、17个完整损伤关键帧的候选筛查、客观混合模式与完整endpoint-star虚拟断键能量诊断，以及两组实际原生三维预裂纹小试样。未删除冻结复现/构建源码；新工具共享图分析实现，未改生产C++、材料律、界面卡或接触刚度，未重跑真实V3。

真实V3中≥34/114节点（离散体积尺度≥200/300µm）的合格弱连接凝聚候选均为0，`NO_COHERENT_CHUNK_CANDIDATE`限于明确记录的关键帧与图搜索，不声称穷尽所有子集。仅1ms的三点[1625,1632,1649]满足候选条件，由bulk键90571拴住；κ=0.1698109341、d=0.6999879333，该键1.00425ms自然失效，三点末态又分开。它不是潜在大块。真实游离分量历史最大13点于0.554125ms形成，42条原始内部键仅13条存活，0.5745ms分为9/3/1点；支持局部凝聚性在脱离前已丧失并发生二次裂解，不等于主体整体粉碎。

0.376ms接受态从原无损归档按原SHA精确恢复：574/1429均仍附着，间距28.1995µm、重叠0.4360095dx、排斥力12.6385µN；邻域破坏早于接触峰值，不能凭时间相关性认定接触是主因。首bulk失效在自由表面附近；早期断键端点95.85%位于一个δ内，但表面受载、几何与邻域不完整混杂，未据此引入表面修正。原始V3几何/键与所用exposed_surface源一致。

客观局部极分解诊断发现低于历史伸长起点但剪切指标≥0.15的存活连接，关键帧峰值3365条；0.15仅为比较尺度，未标定critical skew。三点桥键局部星秩不足，剪切为NOT MEASURED，且已受当前拉伸判据软化；未证明大块存在剪切漏检桥。单键ΔU_virtual通过同坐标下完整endpoint-star能量重计算与独立整模型核验，不是Gc或可直接累加的单键耗散。用户提供的Karmakar2025全文已核验，NOSB/Ogden、独立界面弹簧及1.38/2.27长度比不可直接套入当前LPS。

小试样为1792材料点+336固定锚点，2736初始预切；两组仅顶面规定方向z/x不同，其余参数与240µm净位移、2ms半余弦、2ms夹具释放后0.5ms自由积分均相同。夹具释放不删材料键，释放前含夹具点的分量不算自由块。法向/切向分别完成{c0['accepted_steps']}/{c1['accepted_steps']}接受步，末态bulk自然断键{c0['last_row']['bulk_failed']}/{c1['last_row']['bulk_failed']}、界面断键{c0['last_row']['interface_failed']}/{c1['last_row']['interface_failed']}，历史最大自由分量{c0['maximum_free_nodes']}/{c1['maximum_free_nodes']}节点。最大归一化能量残差{100*c0['max_relative_energy_residual']:.6g}%/{100*c1['max_relative_energy_residual']:.6g}%，原生与独立Python力/势能/耗散/κ/d、VV自由DOF及最终真实拓扑核验通过。首次夹具释放接口失败和不满足精确重放要求的批量重放均保留；改为单NVE及原逐步调用重放后核验续算，没有调材料或强行造块。

`ROOT_CAUSE_STATUS=SUPPORTED_HYPOTHESIS`；真实V3的`CHUNK_FORMATION_STATUS=NOT_OBSERVED`（≥200µm尺度）；独立合成试样为`{status['SYNTHETIC_CHUNK_FORMATION_STATUS']}`。13点刚脱离时内部相对动能约0.475nJ而储能约0.0137nJ，COM拖曳几乎不直接耗散相对运动，是二次裂解的待验证机制。尚无唯一因果结论或经标定改动满足新V3准入，故本轮不追加正式V3。下一轮唯一建议：{status['NEXT_RECOMMENDED_CHANGE']}全局时间步/空间收敛、材料Gc及生理碎片尺寸仍未认证；旧G0和全部历史记录保留。

见[本轮图表索引](temp_storage/clot_github_sync_20261015/clot_clearance/results/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z/index.html)、[十项科学问题与实际运行报告](temp_storage/clot_github_sync_20261015/clot_clearance/verification/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z/REVIEW_V3_CHUNK_FRAGMENTATION_V1.md)。\n\n'''
 index=text.index(marker);after=(text[:index]+insert+text[index:]).encode();assert after.decode().replace(insert,'',1)==text
 folder=out/'wsl_document';folder.mkdir(exist_ok=True);(folder/'clot_rupture.before.md').write_bytes(before);(folder/'clot_rupture.after.md').write_bytes(after);(folder/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),after.decode().splitlines(True),fromfile=str(p)+'.before',tofile=str(p)+'.after')))
 assert sha(p.read_bytes())==sha(before),'Concurrent edit detected before write'
 tmp=p.with_name(p.name+'.v3_chunk.tmp');shutil.copy2(p,tmp)
 with tmp.open('wb') as f:f.write(after);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p);actual=p.read_bytes();assert actual==after
 record=dict(status='PASS',actual_WSL_path=str(p),utc=datetime.now(timezone.utc).isoformat(),before_sha256=sha(before),after_sha256=sha(after),readback_sha256=sha(actual),all_historical_text_preserved=True,insert_start_line=text[:index].count('\n')+1,backup='wsl_document/clot_rupture.before.md',diff='wsl_document/clot_rupture.diff')
 (out/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n');print(json.dumps(record,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
