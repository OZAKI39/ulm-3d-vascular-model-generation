"""Update the required WSL notebook and active entry without losing prior text."""
from common import *
import difflib,datetime

def guarded_write(path,old,new,prefix):
    before=hashlib.sha256(old.encode()).hexdigest();folder=O/'document_updates';folder.mkdir(exist_ok=True)
    (folder/(prefix+'.before.md')).write_text(old)
    assert sha(path)==before,'Document changed after it was read; inspect before merging'
    path.write_text(new);assert path.read_text()==new
    (folder/(prefix+'.after.md')).write_text(new)
    diff=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(path)+'.before',tofile=str(path)))
    (folder/(prefix+'.diff')).write_text(diff)
    return dict(path=str(path),before_sha256=before,after_sha256=sha(path),diff=str((folder/(prefix+'.diff')).relative_to(R)))

def main():
    notebook=Path('/home/lzy/projects/clot_rupture.md');assert notebook.exists(),'Run this update in WSL with Y_CLOT_PROJECT_ROOT set'
    assert json.loads((O/'ARTIFACT_VALIDATION.json').read_text())['status']=='PASS'
    result=json.loads((O/'NATIVE_V3_RESULT.json').read_text());q=result['objective'];native=result['native'];last=result['last_row'];status=json.loads((O/'FINAL_STATUS.json').read_text())
    old=notebook.read_text();heading='## Fragment Surface Coupling and Physical Heterogeneity V3';assert heading not in old
    linkroot=os.path.relpath(R,notebook.parent);rel=f'{linkroot}/verification/{TASK}/{RUN}';figrel=f'{linkroot}/results/{TASK}/{RUN}'
    pc=lambda x:'未定义' if x is None else f'{100*x:.6g}%'
    block=f'''{heading}

**真实碎块表面求积与相对流体作用（{RUN}）**

先回放 V2 的七个末态目标块，没有重跑旧动力学。五个零加载块中，U1164、U1804、H3809、H3846 没有完整存活、正体积的辅助四面体；U1848 有5个四面体，但16个候选面未通过开放判定。U1884 仅2个加载面，H3769 仅1个，后者另有1个开放但不可达面。没有发现这七块的零加载来自原始ID所有权或旧公式节点映射错误。

保留原有受支持三角面 L1；新增仅用于真正游离分量的 L2 原始控制体积暴露面，检查完整存活非局部键的跨面桥接、当前开口、局部仿射变形、L1 重叠、对偶单元和真实四面体障碍。不能只由近邻断键认定裂面；无法确认的面保留 UNRESOLVED。七块均得到非零的部分有效求积，但不代表完整物理湿润面积。没有凸包外壳、恢复断键或按预期尺寸切分材料。

同一组局部仿射权重配对面位置、面速度和节点力。L1权重非负；L2外侧面采用受限带符号权重（最小−0.25、绝对和≤2），以保持真实面中心的力矩和虚功关系。离线最大合力/力矩/功率误差约5.42×10⁻²⁰N、5.82×10⁻²²N·m、2.71×10⁻²⁰W；这些是离散恒等式检查，不能认证局部牵引分布精度。原生少量微步及重启续步通过，κ/d和接受态提交语义保留。

正式仅一个 H：{native['status']}，接受到{native['time_s']*1e3:.7f}ms；同一8125材料点、1875锚点、原始图、缺口、10强区和材料卡，同一1000Pa均匀+X驱动。β=40Pa·s/m，使用当前切面相对滑移；2.8ms后外部驱动为零但被动阻力保留。唯一主要模型变化为表面表示及随之产生的流体作用。原生库未改，无新Uniform、幅值或随机种子扫描。

末态20–200节点游离块{q['n_fragments_20_200']}个，时程最多{q['max_target_count']}个；{'节点数条件达到至少5个，其他质量与存活指标仍分别报告' if q['n_fragments_20_200']>=5 else '未达到至少5个的目标，不改变成功定义'}。最大游离块{max(result['final_free_sizes'],default=0)}点，目标质量/初始质量{pc(q['mass_fraction_target_initial'])}，细屑1–2点/游离质量{pc(q['fine_debris_mass_fraction_detached'])}，总游离质量/初始质量{pc(q['total_detached_mass_fraction'])}。撤载时目标块的原成员集合完整存活{q['exact_target_sets_surviving_drive_off']}个，末态目标区间后代{q['target_descendants_of_drive_off_targets']}个。计数、形貌、掉点及二次分裂由存活原始图与原始成员ID分别判断。

Wdrive/Wdrag={last['W_drive_J']*1e9:.7g}/{last['W_drag_J']*1e9:.7g}nJ；最大相对数值余额{native['max_relative_energy_residual']*100:.7g}%，R不并入D。新裂面与L2外功是可能重叠的诊断子集，不重复加账。所有功均沿实际节点位移积分。正式接受{native['accepted_steps']}步、预测拒绝{native['predictive_rejections']}次，CPU{native['cpu_s']/3600:.5g}h、wall{native['wall_s']/3600:.5g}h、峰值RSS{native['peak_RSS_bytes']/2**30:.5g}GiB。

物理尺度仍 MICROSTRUCTURE_SCALE_UNCALIBRATED：Liu的0.1–1mm为断裂过程区；microCT研究中的数百微米RBC聚集与数十微米边界不等于当前强纤维蛋白相的实测相关长度。正式保留冻结10强区；米制相关长度生成器仅离线检查100/200/400µm假设及dx=25/50/100µm取样，没有用20–200目标反推强区大小。原始跨不同强区直接存活连接3004条，强区不等于预设碎块。

科学标签 EXPLORATORY_UNCALIBRATED、PRESCRIBED_SHEAR_MECHANICAL_SURROGATE、CURRENT_SURFACE_APPROXIMATE_COUPLING。小尺度碎屑和部分凝聚体的表面仍未解析；零数值载荷不是零物理阻力。单次新旧轨迹对比不支持把某次二次裂解归因于阻力单独作用，也不构成已验证声流双向FSI。

详见[七问报告]({rel}/REVIEW_CLOT_FRAGMENT_SURFACE_COUPLING_V3.md)、[八组科研图与两个1080p真实视频]({figrel}/OPEN_RESULTS.html)。实际接受态坐标位移×1，固定全景范围，原始ID载荷与帧来源可追踪；详细后续研究排序以本轮报告中的实际证据为准。

'''
    marker='# ⑥ 外加剪切';assert old.count(marker)==1;new=old.replace(marker,block+marker,1)
    notebook_record=guarded_write(notebook,old,new,'clot_rupture')
    entry=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';previous=entry.read_text();assert '# 当前碎块表面耦合 V3' not in previous
    lead=f'''# 当前碎块表面耦合 V3（{RUN}）

当前任务 {TASK}，analysis/{TASK}/；新L2控制体积求积补充冻结V2受支持三角面，完整存活原始图定义碎块，近似相对流体作用。

状态={status['status']}；正式H={native['status']}；实际接受{native['time_s']*1e3:.7f}ms，末态/历史最多20–200节点游离块={q['n_fragments_20_200']}/{q['max_target_count']}。

正式V2几何、10强区、材料、1000Pa时程保持；本轮唯一正式H已使用，不追加重跑。MICROSTRUCTURE_SCALE_UNCALIBRATED、PRESCRIBED_SHEAR_MECHANICAL_SURROGATE、CURRENT_SURFACE_APPROXIMATE_COUPLING。节点受力支撑率不等于连续表面完整性；未解析水动力不能写成物理零阻力。

先读 verification/{TASK}/{RUN}/FINAL_STATUS.json、REVIEW_CLOT_FRAGMENT_SURFACE_COUPLING_V3.md、DELIVERY_VALIDATION.json、ARTIFACT_VALIDATION.json。图及两个MP4：results/{TASK}/{RUN}/OPEN_RESULTS.html。代码与紧凑证据在本轮独立Git提交，完整接受态及大二进制留在归档，不自动推送。

## 历史入口原文（保留）

'''
    entry_record=guarded_write(entry,previous,lead+previous,'ACTIVE_DEVELOPMENT_ENTRYPOINTS')
    dump(O/'WSL_DOCUMENT_UPDATE.json',dict(status='PASS',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),notebook=notebook_record,active_entry=entry_record,prior_sections_preserved=True))
    print('WSL_DOCUMENTS_UPDATED_WITH_BACKUPS',flush=True)

if __name__=='__main__':main()
