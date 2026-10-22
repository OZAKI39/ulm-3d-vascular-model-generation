"""Minimal append inside the established narrative; preserve every old byte."""
from pathlib import Path
import hashlib,difflib,json,re
ROOT=Path(__file__).resolve().parents[1]
p=Path('/home/lzy/projects/microbubble_aggregation.md')
backup=p.with_name('microbubble_aggregation.md.bak_F10T_Vision2_20261016')
out=ROOT/'evidence/progress_document';out.mkdir(exist_ok=True)
current=p.read_bytes()
before=backup.read_bytes() if backup.exists() else current;text=before.decode()
marker='**F10T-Vision-2 把实验视觉推进到人工标注与真实在线候选测试。**'
assert marker not in text, 'Original backup already contains this section'
audit=json.loads((ROOT/'results/annotation_audit_current.json').read_text())
addition=f'''**F10T-Vision-2 把实验视觉推进到人工标注与真实在线候选测试。**复用已经能在Windows正常工作的Chronos接收器，通过独立WSL视觉进程和有界最新帧接口处理真实1280×1024预览；没有重写RTSP，也没有修改原相机、Tektronix工程。准备72张原始录像候选和12张实时帧，共84张，按独立时间块划分训练、验证和测试，直播帧单独隔离；真实曝光时间、光学倍率和像素到微米标定仍未知，不能把像素轨迹直接当成物理速度。

**标注工具的实际卡顿也保留了修复证据。**初版逐帧重建窗口、同步加载SAM及复杂图层切换妨碍操作，改用napari原生手工轮廓、同窗口换帧、中文步骤、版本化草稿和自动保存；原版SAM只在需要时后台运行，可取消且有超时。真实窗口中的多对象描画、保存重开、跨分辨率换帧及后台取消检查通过，换帧约0.11秒。修复时已有两张草稿原样保留；截至2026年10月16日15:18（苏黎世）报告快照，用户已正式确认1张训练帧中的39个单泡，另有1张草稿，验证／测试帧尚未确认。软件测试、自动候选均不是真值，独立评价继续EXPERT_GROUND_TRUTH_PENDING。

**速度和识别是否正确分别验收。**OpenCV＋Norfair在固定原始256×384 ROI的20秒真实在线测量中处理约60.08 Hz，新鲜图像候选约59.93 Hz，完整软件处理p95约8.37 ms；仍有跳帧和大于16.7 ms的个别帧，不是逐帧硬实时。唯一AI检测路线Cellpose cyto3＋Norfair恢复后实际处理约11.48 Hz，p95约107.63 ms，230帧中151帧的处理结果超过软件新鲜度上界，当前场景无可用候选。两次直播先后运行且后台负载不同，不能据此比较识别准确率；另完成同一录像同帧同ROI的12帧全图／ROI对照和90帧因果离线回放，候选数差异很大，但没有人工真值就不判断谁漏检或误检更少，也没有微调。

**真实微泡与团形状仍未准入。**Norfair仅用短时观测关联，丢失不发送旧坐标，跨session和时间间隙重置；原始图像区域的面积、Rg和形状张量可以计算，但正式group_shape仍为null，单泡识别、短轨迹精度、团分割和团附近单泡分布尚待人工确认。已验证微泡／团更新率均为0，不能把60 Hz候选称为闭环视觉。四组英文PNG/PDF及两段实际H264视频已保存，视频保留真实源时间顺序，不插帧；标注说明、失败、中断、源码身份、哈希与性能日志位于external_connect/vision_tracking_research/smart_online_tracking_v2。未进行LED延迟测量，CAMERA_TO_CONTROLLER_END_TO_END_LATENCY=NOT_MEASURED；仅mock，没有PZT/AFG自动控制。下一步先完成少量完整人工帧和8帧短轨迹的确认，再冻结评价和决定是否有必要开展一次有限微调，training_ready=false。

'''
needle='---\n\n### 尾声：故事里还没讲完的部分'
assert text.count(needle)==1
if current!=before:
 current_text=current.decode();start=current_text.index(marker);end=current_text.index(needle,start)
 assert current_text[:start]+current_text[end:]==text, 'Unrelated user edits detected; do not overwrite'
 (out/'after_initial_delivery.md').write_bytes(current) if not (out/'after_initial_delivery.md').exists() else None
after=text.replace(needle,addition+needle)
head=lambda s:re.findall(r'^#{1,6}\s.*$',s,re.M)
assert head(text)==head(after)
assert after.replace(addition,'',1)==text
if backup.exists():assert backup.read_bytes()==before,'Existing backup differs; do not overwrite'
else:backup.write_bytes(before)
(out/'before.md').write_bytes(before);(out/'after.md').write_bytes(after.encode())
(out/'change.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),after.splitlines(True),fromfile='before',tofile='after')))
p.write_bytes(after.encode())
h=lambda b:hashlib.sha256(b).hexdigest()
(out/'audit.json').write_text(json.dumps({'before_sha256':h(before),'after_sha256':h(after.encode()),'backup':str(backup),
                                       'headings_preserved':True,'old_content_preserved':True,'added_paragraphs':4},indent=2))
print('Updated narrative; original content and headings preserved.')
