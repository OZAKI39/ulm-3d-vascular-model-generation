"""WSL-only guarded scientific notebook / current-entry update after delivery."""
from common import *
import difflib

def guarded(path,new,backup):
    before=path.read_bytes();before_sha=sha(path);backup.parent.mkdir(parents=True,exist_ok=True)
    assert not backup.exists();backup.write_bytes(before)
    assert sha(path)==before_sha
    tmp=path.with_name(path.name+'.y_current_method.tmp');tmp.write_text(new);os.replace(tmp,path)
    difference=''.join(difflib.unified_diff(before.decode().splitlines(True),new.splitlines(True),fromfile=str(path)+'.before',tofile=str(path)))
    backup.with_suffix('.diff').write_text(difference)
    assert path.read_text()==new
    return dict(path=str(path),before_sha256=before_sha,after_sha256=sha(path),backup=str(backup),diff=str(backup.with_suffix('.diff')))

def main():
    assert str(R).startswith('/home/lzy/projects/')
    z=json.loads((O/'APPLICATION_RESULT.json').read_text());v=json.loads((O/'VISUALIZATION_VALIDATION.json').read_text())
    report=Path('temp_storage/clot_github_sync_20261015/clot_clearance/verification')/TASK/RUN/'REVIEW_Y_CLOT_CURRENT_METHOD_V2.md'
    index=Path('temp_storage/clot_github_sync_20261015/clot_clearance/results')/TASK/RUN/'OPEN_RESULTS.html'
    marker=f'**当前方法移入原 Y 血管不规则血块（{RUN}）**'
    addition=f'''\n\n{marker}：按用户要求执行一次更正后的正式工况；前一附着面误载荷适配运行已主动停止并排除科学结论。保留原1927材料点、1561固定锚点、157637原始键、全部参考体积与加权体积，以及无预制缺口初态。沿用G=940Pa、K=100000/3Pa和强/弱变分软化；将当前材料场按物理尺度平移裁剪到原血块，得到882强点、8个原强区标签，比例45.7706%，不重新调参凑回规则试样的33.3292%。这不是患者影像几何，也不是只有几何改变的严格单因素对照。

采用当前L1/L2部分表面求积、+X切向相对速度代理、1000Pa名义时程和原Y管腔检查。原Y壁20µm有效半径、50N/m罚接触和0.5阻尼比参与；壁功与流体功分列。没有叠加旧Gaussian、COM背景阻力或静态CFD力，没有新CFD、微泡声流或双向FSI。原图含被曲面空隙排除的近邻，表面支撑只保留节点和六条原始边均存在的辅助四面体，不补键。

原生状态 **{z['status']}**，接受到 **{z['accepted_time_ms']:.6f}ms**。末态游离 **{z['detached_nodes']}/1927点（{z['detached_fraction']:.3%}）**、**{z['free_components']}个分量**，最大 **{z['largest_free_nodes']}点**；20–200点分量末态/历史最多 **{z['final_20_200']}/{z['maximum_20_200']}**。本体/界面失效 **{z['failed_bulk']}/{z['failed_interface']}**，最大归一化能量余额 **{z['max_relative_energy_residual']:.6%}**，历史最大有效L2面 **{z['maximum_L2_patches']}**。数值与可视化检查不等于材料标定或完整水动力认证。

红色、Z轴竖直视频为1920×1080、24fps，真实位移×1，仅使用接受态；另有Y血管全景、初态/末态、实际载荷场与动力学曲线。详见[本轮开发报告]({report.as_posix()})和[静态图与视频入口]({index.as_posix()})。\n'''
    p=Path('/home/lzy/projects/clot_rupture.md');old=p.read_text();assert marker not in old
    anchor='**V3 附加红色可视化（2026-10-19）**';i=old.index(anchor);end=old.index('\n',i)
    record=guarded(p,old[:end]+addition+old[end:],O/'wsl_document/clot_rupture.before.md')
    entry=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';old=entry.read_text();head=f'''# 当前：原 Y 血管不规则血块上的 V3 方法应用（{RUN}）

任务 `{TASK}`；状态 **{z['status']}**，接受时间 {z['accepted_time_ms']:.6f}ms。仅1次正式工况；原1927材料点+1561锚点、无缺口、原材料图保留。原当前强弱场物理平移裁剪，45.7706%强相；1000Pa当前表面剪切代理、L1/L2和原Y壁接触。无新CFD、微泡声流或双向FSI。

末态游离{z['detached_nodes']}点/{z['free_components']}分量；最大{z['largest_free_nodes']}点；20–200点块末态/历史最多{z['final_20_200']}/{z['maximum_20_200']}。材料与完整湿润面未标定，数值验收不替代物理认证。

源码 `analysis/{TASK}/`；先读 `verification/{TASK}/{RUN}/FINAL_STATUS.json`、`REVIEW_Y_CLOT_CURRENT_METHOD_V2.md`、`DELIVERY_VALIDATION.json`。图与1080p红色视频：`results/{TASK}/{RUN}/OPEN_RESULTS.html`。新实验使用新目录；保留本轮全部接受数据与冻结哈希，不自动重跑或推送GitHub。

## 历史入口原文（保留）

'''
    er=guarded(entry,head+old,O/'wsl_document/ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
    dump(O/'WSL_DOCUMENT_UPDATE.json',dict(status='PASS',notebook=record,entrypoint=er,physics_recomputed=False))
    print('GUARDED_WSL_RECORDS_UPDATED',record,er,flush=True)

if __name__=='__main__':main()
