"""Insert one Chinese F10R paragraph; verify all existing text/headings survive."""
from pathlib import Path
from datetime import datetime, timezone
import difflib
import re
from .common import BASE, read, write, sha, identity


def main():
    doc=Path('/home/lzy/projects/microbubble_aggregation.md')
    backup=doc.with_name(doc.name+'.bak_F10R_20261014')
    out=BASE/'progress';out.mkdir(exist_ok=True)
    gate=read(BASE/'analysis/FINAL_GATE.json');cost=read(BASE/'analysis/COSTS.json')
    assert gate['ENGINEERING_PASS']=='PASS'
    assert read(BASE/'figures/VISUAL_REVIEW.json')['status']=='PASS'
    before=doc.read_bytes();text=before.decode('utf-8')
    assert '**F10R ' not in text,'Already updated; do not duplicate the historical paragraph'
    assert sha(doc)=='355f53f0f89c7f5d53392cd53dec2e7c007902015f40447c49290d8355988b7a','Concurrent document edit: review before touching'
    if backup.exists():assert backup.read_bytes()==before,'Never overwrite a different backup'
    else:backup.write_bytes(before)
    paragraph=f'''**F10R 将连续四路幅相扩展成面向未来实验操作的离线仿真序列。**从最新公开F10Q同步HEAD核查950个科学源码/配置等文件与完整本地科研版本字节相同，在独立worktree开发；有界引用审计仍无安全废弃对象可删，原M0/B2/硬球、失败和历史结果均保留。按用户决定暂不开发AFG1062、功放、USB等真实仪器通信，只实现ONE_HOT_PZT、AXIS_PAIR_CONTROL及FOUR_CHANNEL_INDEPENDENT三种共同物理入口；每段记录编号、起止、2.5 MHz、四路启停、Pa-equivalent幅值和相位，1500源上限不是Vpp或安全电压。采用IDEAL_SCHEDULED_SEGMENTS与公共相位时钟，非网格事件缩短步长对齐，实际两段控制在1 ms切换并保存前后逐泡声压，没有虚构实测USB时延。保留原32个固定ID合法接触初态、充水无流、R3封闭X壁面、缓存、APECSS/F2及原五次Sequential QP门槛；新入口重放原X的全部2001帧坐标逐位相同，最大差0 m。X绝对目标(0.2,0,0) µm和0.075 µm容差未动：OFF严格静止，原固定X、同轴X对、两段X的2 ms位移分别约0.200896、0.131844、0.200675 µm，误差约0.000896、0.068156、0.000675 µm，基础、COM、形状、连通与ID均通过；单片X−/X+各1 ms的对称序列近乎相消，未到目标。两段只比原固定控制改善约0.222 nm误差，不能概括为普遍更优。±5°初始团旋转均观察到约0.200894 µm，但未对旋转另做完整半步/缓存准入，不当独立实验样本。仅改X+相位的场诊断给出约3.73%/度而非粗算6%/度；−6/+6/+10°完整2 ms分别移动0.155457/0.244999/0.273417 µm，观察均在容差，+10°完整数值对照通过但余量仅1.583 nm，其他场扫描点未伪造轨迹。

新Y开跑前单片诊断确认7 mm、2.5 MHz圆片的轴上零点距4.504167 mm，而Y片距中心4.5 mm；同1500 Pa-equivalent源的团心Y/X声压约13.778/2785.553 Pa，效率约1/202，故所有Y结论标HOMOGENEOUS_WATER_PROXY_SPECIFIC，不外推真实PDMS。旧Y在1.9 ms的二次力/接触反力约数百pN，末步活动接触对64/64/62；固定几何只换缓存力差很小，而不同末态同缓存的二次力差达91–134 pN，显示接触重排、几何反馈和步长/场差异共同敏感，不能归为唯一原因。四个0.5 ms候选后选择原Y控制X源降为0.25倍的新弱Y，实际泡处声压约416–471 Pa，独立2 ms/半步/细场移动约0.012944 µm、目标误差0.187056 µm，未达Y目标，但COM、形状、连通和ID均通过，最大ID步长/缓存差仅约0.00000420/0.00097581 nm；旧Y的COM通过、形状及ID失败结论原样保留。新弱Y的Rg只收缩约0.295 nm、X/Y比近1，虽超过所测差异，仍只是演示接触初态的微小松弛，不是指定形态或聚合融合。等幅同轴对只算场的连续压力峰相位斜率X约−0.843 µm/度，Y外侧强峰约−0.774 µm/度，符号按φplus−φminus；原Y近中心选峰混淆分支的错误拟合已保留并作废，压力峰尚未证明是三维稳定陷阱，也未做长时陷阱轨迹。本轮25条有限轨迹完成各自窗口，最长2 ms；已计时的科学执行、分析及必要核验约{cost['total_scientific_s']:.3f} s，其中远端四个独立Supervisor进程约{cost['remote_campaign_wall_s']:.3f} s，纯GPU核时未单独剖析，传输/绘图/开发不混入求解耗时。170个远端文件逐项SHA回传；归档标签的输出循环覆盖问题只修派生元数据和OFF零参考判定，原始轨迹及SHA不改。四组英文PNG/PDF实际打开检查，完整输入、源码快照、逐步反力和失败/NOT_RUN证据归档。下一轮优先验证移动X压力峰的三维力稳定性并有限扩大X距离；Y先确认真实声场和电压声压标定，再研究低风险干涉控制，不急于任意团形状。持续 **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED，training_ready=false**；详见独立F10R工作区docs/FAST_PZT_F10R_DEVELOPMENT.md及results/straight_tube_exp1/f10r/。

'''
    marker='---\n\n### 尾声：故事里还没讲完的部分'
    assert text.count(marker)==1
    after=text.replace(marker,paragraph+marker).encode('utf-8')
    assert after.replace(paragraph.encode('utf-8'),b'',1)==before
    headings=lambda s:re.findall(r'^#{1,6} .*$',s,re.M)
    assert headings(text)==headings(after.decode('utf-8'))
    (out/'before.md').write_bytes(before);(out/'after.md').write_bytes(after)
    (out/'change.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),after.decode('utf-8').splitlines(True),fromfile='before',tofile='after')),encoding='utf-8')
    doc.write_bytes(after)
    write(out/'UPDATE_AUDIT.json',dict(status='PASS',identity=identity(),updated_utc=datetime.now(timezone.utc).isoformat(),
          path=str(doc),backup=str(backup),before_sha256=sha(out/'before.md'),backup_sha256=sha(backup),after_sha256=sha(doc),
          previous_bytes_preserved=True,headings_before=headings(text),headings_after=headings(after.decode('utf-8')),
          diff_sha256=sha(out/'change.diff'),insertion_only=True))
    print('PASS: progress backup + insertion only; all prior bytes and headings preserved')


if __name__=='__main__':main()
