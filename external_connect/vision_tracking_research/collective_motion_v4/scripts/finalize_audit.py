"""History/progress/artifact audit. Run after figures/videos/report review is complete."""
import sys,json,time,re,difflib,shutil,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
def main():
 before=json.loads((ROOT/'evidence/history_before.json').read_text());changed=[]
 for rel,digest in before['files'].items():
  p=Path(before['root'])/rel
  if not p.exists() or sha(p)!=digest:changed.append(rel)
 assert not changed,changed
 dump(ROOT/'evidence/history_after.json',dict(checked=len(before['files']),changed=changed,baseline=before['git_head'],utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
 progress=Path('/home/lzy/projects/microbubble_aggregation.md');backup=progress.with_name(progress.name+'.bak_F10X_20261017');folder=ROOT/'evidence/progress_document';folder.mkdir(exist_ok=True)
 if not backup.exists():shutil.copy2(progress,backup)
 old=backup.read_text();current=progress.read_text();assert 'F10X-PIV：全通道观测' not in old
 a=json.loads((ROOT/'results/delivery_analysis.json').read_text());off=next(r for r in a['performance'] if r['name']=='verified_offline_regular');on=next(r for r in a['performance'] if r['name']=='verified_replay_regular')
 addition=f'''**F10X-PIV：全通道观测、成熟后端加速与仍未解决的弱对比问题（2026-10-17）**

本轮在F10W基线上独立开发collective_motion_v4，实际编译OpenPIV官方Rust、运行RTX4090上的TorchPIV及CuPy官方PIV，没有重写FFT/峰值核心或训练新AI。新版Rust与SciPy同条件位移差约3.55e-15；CuPy全范围规则批量最快，但FFT边界/有效支持不同，不能称为完全相同算法。通道可观察mask为805495px、1283个全窗安全位置，是历史ROI面积8.19倍；初版中央纹理排除遮蔽了真实泡，审计后修正，旧失败保留。人工部分审核503/505中心被活动窗口覆盖，而合成低对比仅10/120，不能声称全通道所有泡都可靠发现。

完整937/937帧离线扫描及四种调度对照已完成；正式窗口支持足迹并集平均占观察mask {100*a['coverage_mean']:.2f}%，不是逐像素水流真值。因果背景gap3、99.9ms区间在开发/验证后冻结；独立程序种子final44条接受30、错接4（对称旋转和局部相反运动），局部增强对照错接11。多遍final错接1但更慢，保留作下一轮预登记参考，不用final调参。仍仅一段真实录像，没有跨实验真实速度验证。

正式离线完整处理p50/p95为{off['p50']:.2f}/{off['p95']:.2f}ms；实际PTS定速规则回放处理{on['frames']}/937、丢{on['dropped']}，p95 {on['p95']:.2f}ms，新鲜有效群帧{on['fresh_group_frame_hz']:.2f}Hz。队列最大1，断帧重新预热，TTL过期为空值；33.3ms与60fps均未准入。局部群是空间候选集合，不是物理融合或跨帧身份，六参数仍按F10U矩阵对数约定。解析Jacobian试验在稀疏Huber问题出现不同优化路径，发现后对鲁棒活跃情况恢复原求解器；整片PIV和六参数已严格核对，失败版本不用于正式图。

已交付9组英文PNG/PDF、两段各937帧全长视频及一段60帧selected-clips后端比较，保存实际PTS、全解码、源码身份、SHA、视觉审查和失败证据。目录为独立归档mevo_vision_f10x_20261017的collective_motion_v4，分支dev/f10x-full-channel-piv-20261017。主要后续工作是独立显微运动真值、弱对比召回、非仿射拒判、CPU局部拟合与序列化，不是扩大模型或自动控制。EXPLORATORY、PHYSICS_UNCALIBRATED、training_ready=false；没有相机连接、PZT/AFG动作或世界模型训练。

---

'''
 marker='### 尾声：故事里还没讲完的部分';assert old.count(marker)==1
 new=old.replace(marker,addition+marker)
 if 'F10X-PIV：全通道观测' not in current:assert current==old;progress.write_text(new)
 else:assert current==new,'Refuse silently replace a changed progress update'
 headings=lambda s:re.findall(r'^#{1,6} .*$',s,re.M)
 assert headings(old)==headings(new)
 (folder/'before.md').write_text(old);(folder/'after.md').write_text(new);(folder/'change.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='before.md',tofile='after.md')))
 dump(folder/'audit.json',dict(path=str(progress),backup=str(backup),before_sha256=sha(backup),after_sha256=sha(progress),headings_preserved=True,headings=headings(new),historical_text_preserved=new.replace(addition,'')==old))
 # Manifest excludes itself and volatile ignored caches/environments; large artifacts remain addressable.
 hashes={}
 for directory in ['src','scripts','configs','tests','results','figures','videos','logs','evidence']:
  for p in (ROOT/directory).rglob('*'):
   if p.is_file() and '__pycache__' not in str(p) and '.pytest_cache' not in str(p) and p.suffix not in ['.tar','.pyc'] and not p.name.startswith('.'):
    hashes[str(p.relative_to(ROOT))]=dict(sha256=sha(p),bytes=p.stat().st_size)
 for p in ROOT.glob('*.md'):hashes[p.name]=dict(sha256=sha(p),bytes=p.stat().st_size)
 for name in ['FINAL_GATE.json','RECOVERY_MANIFEST.json']:
  p=ROOT/name
  if p.exists():hashes[name]=dict(sha256=sha(p),bytes=p.stat().st_size)
 dump(ROOT/'INPUT_OUTPUT_SHA256.json',dict(source_video=dict(path=str(VIDEO),sha256=sha(VIDEO)),files=hashes,scope='Self and final Git identity excluded to avoid recursive hash; bytewise hashes of actual delivered files',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
 print('history preserved',len(before['files']),'manifest files',len(hashes),'progress headings',len(headings(new)),flush=True)
if __name__=='__main__':main()
