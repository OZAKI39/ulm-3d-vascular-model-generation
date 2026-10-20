"""Guarded local note: distinguish visual markers from the physical discretization."""
from pathlib import Path
import json, hashlib, difflib, datetime, os

R=Path(__file__).resolve().parents[3]
TASK='y_clot_v4_reference_style_v1'
RUN='20261020T084000Z'
O=R/'verification'/TASK/RUN


def main():
    validation=json.loads((O/'VALIDATION.json').read_text())
    assert validation['status']=='PASS'
    path=Path('/home/lzy/projects/clot_rupture.md')
    before=path.read_bytes();text=before.decode()
    title='## Y-Clot V4 Reference-Style Visualization and Point-Cloud Scope'
    assert title not in text
    base=f'temp_storage/clot_github_sync_20261015/clot_clearance'
    section=f"""{title}

**按指定旋转视频风格重绘（{RUN}）**：本轮仅使用已完成Y血块V4的0–3.6ms接受态，未启动物理求解。白底、红色强弱相节点、蓝灰锚点、细存活键、六面网格、世界Z竖直、低角度正交相机与单一图例沿用20261019T205403Z参考。4K材料标记14px、锚点8px；本次标记仅表达点的位置与相区，不再以球半径表达Vi。非均匀Vi、质量及所有物理状态未改。

**离散方式名称澄清**：1927材料点中只有412个内部点被扰动；1515个边界/界面材料点保留初始规则排列，1561固定锚点亦规则。因此应称“保留规则边界的内部扰动点云”，不能称为整个血块的完整不规则离散。边界材料点只在初始几何构造中不扰动，后续动力学仍可运动，只有锚点固定。原V4任务4.1节允许边界固定，当前实现采用此分支；78.6%的材料点未扰动，参考风格的短键加深又突出规则连接。单纯重绘不改变点位。

旋转视频为初态相机环绕，18s/540帧；实际碎裂视频读取149个接受态、501视频帧，按真实时刻标注，无坐标插值或位移放大。两类均导出4K/1080p。载荷图直接读取已核验的实际面力，时刻1.95187ms，统一箭头比例并标注2µN；最大牵引1001.819Pa。所有输入哈希保持，四个视频完整解码通过。

见[本轮结果]({base}/results/{TASK}/{RUN}/OPEN_RESULTS.html)、[重绘说明]({base}/analysis/{TASK}/REVIEW.md)及[点云范围核对]({base}/verification/{TASK}/{RUN}/POINT_CLOUD_REGULARITY.json)。此前Vi球体版本仍为历史输出，不覆盖。

"""
    at=text.index('# ⑥ 外加剪切')
    after=(text[:at]+section+text[at:]).encode()
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    audit=O/'notebook_update'/stamp;audit.mkdir(parents=True,exist_ok=False)
    (audit/'before.txt').write_bytes(before);(audit/'after.txt').write_bytes(after)
    (audit/'update.diff').write_text(''.join(difflib.unified_diff(
        text.splitlines(True),after.decode().splitlines(True),
        fromfile=str(path)+' before',tofile=str(path)+' after')))
    assert path.read_bytes()==before,'Concurrent document change'
    temp=path.with_name(path.name+'.reference_style.tmp');temp.write_bytes(after);os.replace(temp,path)
    metadata=dict(status='UPDATED',path=str(path),backup=str(audit),
        before_sha256=hashlib.sha256(before).hexdigest(),after_sha256=hashlib.sha256(after).hexdigest())
    (audit/'HASHES.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (O/'NOTEBOOK_UPDATE.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('NOTEBOOK_UPDATED',audit)


if __name__=='__main__':
    main()
