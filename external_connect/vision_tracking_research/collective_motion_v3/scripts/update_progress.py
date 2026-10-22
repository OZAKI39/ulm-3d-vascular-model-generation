"""Authorized minimal Chinese progress addition, backups and heading checks."""
import sys,re,shutil,difflib,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
def main():
    assert json.loads((OUT/'FINAL_GATE.json').read_text())['ENGINEERING_DELIVERY_PASS']
    path=Path('/home/lzy/projects/microbubble_aggregation.md');backup=path.with_name('microbubble_aggregation.md.bak_F10W_20261017');before=path.read_text();assert '**2026年10月17日的 F10W-Vision' not in before,'ALREADY_UPDATED_REFUSE_DUPLICATE';assert not backup.exists(),'REFUSE_BACKUP_OVERWRITE'
    section='''**2026年10月17日的 F10W-Vision 实际补齐了整段Chronos录像，并完成开源运动模型迁移对照。**原片SHA与历史一致，1476×1120、937帧，视频PTS 0–31.1688 s，画面呈现31.2021 s；此前F10U/F10V各只有120帧局部分析。本轮固定历史ROI、冻结单遍候选增强PIV，937/937逐帧有记录，没有省略中间时段；887帧通过质量门槛、50帧拒判，其中2帧warmup、26帧向量不一致、22帧支持不足。94.664%是可观测质量通过率，不是真实准确率。候选数量24–81、首尾72→31，末段支持下降和六参数突变保留，不能把可见计数、中心或二维散度当成已验证的材料聚集。

使用官方LRAFT-ST学生权重、CoTracker3在线权重、UnLiteFlowNet-PIV权重以及POT Partial/Unbalanced OT，未重写骨干或训练新六参数CNN。新开发46、验证92、校准46、最终138条已知视觉运动序列采用不同源帧、图块区域和种子，但仍只有一段录像；历史F10U/F10V另做238次诊断，不再称独立盲测。最终冻结风险门槛后，PIV接受98/138且3个错误、LRAFT55且1错误、CoTracker105且4错误、Unbalanced OT87且0错误；UnLite未找到合法门槛，0接受。共同89例中CoTracker和Unbalanced OT部分分量精度优于PIV，但有限样本零错误不等于真实安全。UnLite官方8例可运行，迁移微泡后幅值和假接受失败仍保留；本轮自监督适配NOT_RUN，没有伪造训练收益。

完整原片的四路因果定速回放均真实输入937帧。修复逐帧重复权重哈希后，PIV处理860帧、丢77帧，完整管线p95约45.88 ms，新鲜有效六参数21.40 Hz；LRAFT约55.63 ms、5.55 Hz；CoTracker约85.56 ms、1.60 Hz；Unbalanced OT约12.14 ms、2.12 Hz。OT很快的拒判不代表更高有效吞吐；CoTracker当前适配器需要重建16帧历史，不能把该性能当作官方持续流式实现的极限。原慢版本、数值错误、Schema失败及全部科学误判均保留。单槽最新帧、断帧清空、100 ms TTL、本机消费者和v2→v1兼容测试通过，真实相机与端到端延迟未测。

全视频离线五配置及全画幅扫描146.744 s，新数据生成7.109 s，正式模型对照63.626 s，四路修正后回放125.776 s；GPU实际使用RTX4090，基准CUDA峰值分配约1516.8 MiB。11组英文PNG/PDF、两段完整937帧视频和一段明确selected clip的64帧比较视频均已编码、完整解码、保存PTS和SHA并视觉检查。交付位于独立归档工作树mevo_vision_f10w_20261017的collective_motion_v3，分支dev/f10w-pretrained-full-video-20261017；5214项历史保护文件不变。

本轮结论仍是保留PIV观测主线，CoTracker3和Unbalanced OT保留为有潜力的研究对照。优先增加独立显微录像及已知载物台运动真值，同时研究CoTracker持续缓存和OT稳定求解，再决定有限自监督；不因合成成绩改善而进入自动PZT。保留EXPLORATORY、PHYSICS_UNCALIBRATED、training_ready=false；真实六参数准确率和跨实验泛化未验证，没有连接相机、控制PZT/AFG或训练物理动作世界模型。

---

'''
    marker='### 尾声：故事里还没讲完的部分';assert before.count(marker)==1;after=before.replace(marker,section+marker,1);heads=lambda s:re.findall(r'^#{1,6} .*$',s,re.M);assert heads(before)==heads(after);assert after.replace(section,'',1)==before
    folder=OUT/'evidence/progress';folder.mkdir(parents=True,exist_ok=True);(folder/'before.md').write_text(before);shutil.copy2(path,backup);assert path.read_text()==before
    tmp=path.with_suffix('.md.F10W.tmp');tmp.write_text(after);tmp.replace(path);(folder/'after.md').write_text(after);(folder/'change.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='before.md',tofile='after.md')))
    dump(folder/'audit.json',dict(path=str(path),backup=str(backup),before_sha256=sha(backup),after_sha256=sha(path),heading_structure_equal=True,headings=heads(before),original_content_exactly_preserved=True,identity=identity()))
    print('Updated',path,'backup',backup)
if __name__=='__main__':main()
