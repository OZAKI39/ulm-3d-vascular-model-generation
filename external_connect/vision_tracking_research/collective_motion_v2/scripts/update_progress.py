"""Minimal append before existing epilogue; exact backup/diff/headings audit."""
import sys,json,re,shutil,difflib,datetime
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from vcommon import ROOT,dump,sha
p=Path('/home/lzy/projects/microbubble_aggregation.md');backup=p.with_name(p.name+'.bak_F10V_20261017')
if backup.exists():raise RuntimeError('BACKUP_ALREADY_EXISTS: never overwrite history')
before=p.read_bytes();shutil.copy2(p,backup);old=before.decode('utf-8');marker='### 尾声：故事里还没讲完的部分';assert old.count(marker)==1
addition='''**2026年10月17日的 F10V-Vision 实际训练了学习型群体运动观测器，但没有通过替换PIV的科研准入。**本轮在独立工作树、独立服务器环境中复用F10U的OpenCV检测、PIV/CPD、六参数矩阵对数和观测接口，使用同一Chronos录像的四个不相交素材区域与帧段，生成8192条训练、1024条验证、1024条校准和1536条最终测试三帧序列。图块、固定背景、漏检、失焦、启停及真正非仿射运动均保留；热图来自实际检测，不是真值中心。只有一段录像，因此仍不能声称跨实验泛化。

Direct6灰度、灰度加热图、纯热图、抗热图扰动和GeoMotion分别进行了两种随机种子训练，另有两帧和去局部监督对照，共12项、258个epoch、11094个优化步，正式训练墙钟约11.75分钟。最佳权重和拒判阈值按验证集冻结，校准使用独立分组，最终测试后没有重新调参。学习模型在换素材后出现小运动先验/输出塌缩，没有找到满足原定5%联合错误风险且至少接受20样本的验证阈值，因此正式AI输出全部拒判；不能把零接受说成零错误成功。

同一最终测试上，单遍PIV接受1111/1536，其中84个全局运动样本超容差、1个非仿射样本误接受，共85个失败。相同25%诊断覆盖率下，PIV风险约0.52%，融合Direct6约86.98%，GeoMotion约86.46%，不支持用AI替换PIV。独立校准区间较宽，较高覆盖率并不代表有用预测或真实实验校准。真实录像因果回放中，融合AI完整管线p95约21.13ms、GeoMotion约23.56ms，均比同机单遍PIV约44.29ms快；但AI新鲜有效六参数为0Hz，PIV约25.24Hz。真实六参数准确率仍NOT_VERIFIED，相机端到端延迟NOT_MEASURED。

七组英文PNG/PDF、三段按实际PTS播放且完整解码的MP4、全部失败记录、best/last权重及SHA位于独立归档工作树mevo_vision_f10v_20261017的collective_motion_v2。新M4-collective-2.0与旧1.0显式适配、断帧及TTL测试通过；历史4715项保护文件未改变。保留EXPLORATORY、PHYSICS_UNCALIBRATED、training_ready=false；没有连接相机、操作PZT/AFG或训练物理动作世界模型。下一轮先在新的开发素材上解决局部位移学习和固定背景捷径，再增加独立录像及已知载物台运动真值，并采用新的最终留出测试；本轮结论是继续保留PIV观测主线。

---

'''
new=old.replace(marker,addition+marker);p.write_text(new);e=ROOT/'evidence/progress_document';e.mkdir(exist_ok=True);(e/'before.md').write_bytes(before);(e/'after.md').write_text(new);(e/'change.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(backup),tofile=str(p))))
head=lambda s:re.findall(r'^#{1,6}\s+.*$',s,re.M)
assert head(old)==head(new);assert new.replace(addition,'',1)==old
record=dict(path=str(p),backup=str(backup),before_sha256=sha(backup),after_sha256=sha(p),all_original_text_preserved=True,heading_structure_unchanged=True,original_heading_count=len(head(old)),added_characters=len(addition));dump(e/'audit.json',record);print(record)
