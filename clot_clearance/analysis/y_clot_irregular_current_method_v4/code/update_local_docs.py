"""Guarded WSL notebook update; preserve every pre-existing entry verbatim."""
from pathlib import Path
import hashlib,json,datetime,difflib,os
R=Path(__file__).resolve().parents[3];TASK='y_clot_irregular_current_method_v4';RUN='20261020T041100Z';D=R/'analysis'/TASK;O=R/'verification'/TASK/RUN
def guarded(path,transform,audit):
 before=path.read_bytes();after=transform(before.decode()).encode();audit.mkdir(parents=True,exist_ok=False)
 (audit/'before.txt').write_bytes(before);(audit/'after.txt').write_bytes(after)
 assert path.read_bytes()==before,'Concurrent document edit; stop without overwriting'
 tmp=path.with_name(path.name+'.y_v4.tmp');tmp.write_bytes(after);os.replace(tmp,path)
 (audit/'update.diff').write_text(''.join(difflib.unified_diff(before.decode().splitlines(True),after.decode().splitlines(True),fromfile=str(path)+' before',tofile=str(path)+' after')))
 (audit/'HASHES.json').write_text(json.dumps(dict(path=str(path),before_sha256=hashlib.sha256(before).hexdigest(),after_sha256=hashlib.sha256(after).hexdigest()),indent=2)+'\n')
def main():
 q=json.loads((O/'SCIENTIFIC_OUTCOME.json').read_text());v=json.loads((O/'DELIVERY_VALIDATION.json').read_text());assert v['status']=='PASS'
 s=q['statistics']['Y V4'];old=q['statistics']['Y regular'];stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');audit=O/'provenance/document_update'/stamp
 base=f'temp_storage/clot_github_sync_20261015/clot_clearance';title='## Irregular V4 Application to the Original Y-Vessel Clot'
 section=f'''{title}

**原Y形血管三瓣壁附血块的V4实际应用（{RUN}）**：这是对已有工程构造几何的复用，不是患者重建。保留1927材料点、1561锚点、原相区ID与原附着面；1515边界/界面材料点不施加初始扰动，只扰动412内部点；正式动力学中材料点均可运动，仅1561锚点固定。PD积分体积0.240875mm³、质量0.240875mg和45.7706279%强相积分比例保持；连续形状约0.271577mm³与PD积分体积明确分开。没有把长方体V4的质量带入Y案例。

种子20261019、最大扰动12.5µm、最小间距28.4943µm；逐点Vi由共享映射控制体积分，参考键157637→153986。完整非均匀变分力、ρVi质量、不可逆软化和原始图碎块身份沿用V4；L1/L2附着面/开放间隙和原Y管腔检查保留。壁接触刚度50N/m、有效半径20µm、阻尼比0.5，阻尼系数按实际质量计算。有限Y样本仿射剪切响应变化−0.2654%/+0.03025%/+1.41078%，仅是相对离散核验。

唯一新物理轨迹状态 **{q['status']}**，实际接受 **{q['accepted_time_s']*1e3:.8f}ms**；累计{q['accepted_steps']}步。1000Pa指定当前表面剪切协议，2ms上升、0.4ms平台、0.4ms撤载，之后被动拖曳；步长0.25µs，最小事件步0.03125µs。最大归一化能量残差{100*q['maximum_normalized_energy_residual']:.7f}%，壁外功单列，不重复计墙面势能。没有新求解CFD、声流、微泡或双向FSI。

当前末态游离{s['free_nodes']}点/{s['free_components']}个分量，最大{s['largest_nodes']}点；20–200节点游离块末态/时程峰值{s['target_20_200']}/{s['peak_target_count']}，历史Y规则点阵对应{old['target_20_200']}/{old['peak_target_count']}。实际质量按Vi求和；单次离散应用不证明真实微结构、碎块稳定性或声学可实现性。图形仍显示规则的固定外边界，内部才有扰动。

见[开发与科研报告]({base}/analysis/{TASK}/REVIEW_Y_CLOT_IRREGULAR_CURRENT_METHOD_V4.md)及[实际时程、载荷场与Z轴旋转]({base}/results/{TASK}/{RUN}/OPEN_RESULTS.html)。红色点云、世界Z竖直、位移×1；18s旋转视频仅环绕t=0初态，另有实际接受态碎裂视频，二者不混淆。球体只表示Vi，没有生成连续壳层。

'''
 def notebook(text):
  assert title not in text
  at=text.index('# ⑥ 外加剪切');return text[:at]+section+text[at:]
 guarded(Path('/home/lzy/projects/clot_rupture.md'),notebook,audit/'notebook')
 prefix=f'''# 当前：原 Y 血管血块上的 V4 不规则方法应用（{RUN}）

任务 `{TASK}`；状态 **{q['status']}**，实际接受 **{q['accepted_time_s']*1e3:.8f} ms**。1927材料点+1561锚点，412内部点扰动、153986参考键，原外边界/附着面及材料相区保持。积分体积0.240875mm³、质量0.240875mg。

先读 `analysis/{TASK}/FINAL_STATUS.json`、`REVIEW_Y_CLOT_IRREGULAR_CURRENT_METHOD_V4.md`、`REPRODUCE.md`。实际时程、载荷场、4K/1080p Z轴旋转视频在 `results/{TASK}/{RUN}/OPEN_RESULTS.html`；原生证据在 `verification/{TASK}/{RUN}`。仅1条新物理轨迹，历史Y不重跑。EXPLORATORY_UNCALIBRATED、MECHANICAL_SHEAR_SURROGATE；没有新CFD/声流/双向FSI。不得覆盖接受态或自动追加工况。

## 历史入口原文（保留）

'''
 guarded(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',lambda text:prefix+text,audit/'entrypoint')
 (D/'ENTRYPOINT_PREFIX.txt').write_text(prefix)
 (D/'NOTEBOOK_UPDATE.json').write_text(json.dumps(dict(status='UPDATED_WITH_HASH_GUARD_AND_BACKUP',audit_path=str(audit.relative_to(R))),indent=2)+'\n')
 print('LOCAL_DOCS_UPDATED',audit)
if __name__=='__main__':main()
