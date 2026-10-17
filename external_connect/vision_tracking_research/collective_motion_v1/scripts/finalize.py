"""Integrity/progress archive. Idempotent, no alteration of research inputs."""
import sys,json,re,difflib,subprocess,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
import jsonschema

before=json.loads((ROOT/'evidence/protected_inputs_before.json').read_text())
after={p:sha(p) for p in before};changed=[p for p in before if before[p]!=after[p]]
dump(ROOT/'evidence/protected_inputs_after.json',after)
assert not changed,changed
schema=json.loads((ROOT/'configs/collective_state.schema.json').read_text());validator=jsonschema.Draft202012Validator(schema);count=0
for fn in ['causal_replay.json','causal_replay_fast.json']:
    for row in json.loads((ROOT/'results'/fn).read_text())['records']:
        validator.validate(row['state']);count+=1
dump(ROOT/'evidence/schema_validation.json',{'PASS':True,'actual_message_count':count,'schema_sha256':sha(ROOT/'configs/collective_state.schema.json')})

# Explicit causal changes of visible distribution, never interpreted as material velocity.
rows=json.loads((ROOT/'results/real_observation.json').read_text())['rows'];changes=[]
for a,b in zip(rows[:-1],rows[1:]):
    dt=b['pts']-a['pts'];da,db=a['distribution'],b['distribution']
    changes.append({'frame_id':b['frame_id'],'source_pts':b['pts'],'dt_s':dt,'visible_count_delta':db['visible_count']-da['visible_count'],
      'density_change_per_px2_s':(db['density_summary']['mean_per_px2']-da['density_summary']['mean_per_px2'])/dt,
      'visible_mean_change_px_s':(np.array(db['group_center_xy'])-da['group_center_xy'])/dt,
      'material_translation_inferred':False,'ingress_egress_verified':False})
dump(ROOT/'results/distribution_changes.json',{'scope':'visible candidates; changing support/detections may dominate','rows':changes})

progress=Path('/home/lzy/projects/microbubble_aggregation.md');backup=progress.with_name(progress.name+'.bak_F10U_Vision_20261017');ev=ROOT/'evidence/progress_update';ev.mkdir(exist_ok=True)
tag='**2026年10月17日的 F10U-Vision 改为不依赖单泡长期ID的群体观测。**'
text=progress.read_bytes().decode('utf-8')
if tag not in text:
    assert not backup.exists(),'Backup already exists: inspect instead of overwriting'
    orig=text;backup.write_bytes(progress.read_bytes());(ev/'before.md').write_bytes(progress.read_bytes())
    addition='''**2026年10月17日的 F10U-Vision 改为不依赖单泡长期ID的群体观测。**沿用原OpenCV候选检测，在真实Chronos录像120帧中输出可见中心、数量、密度、协方差和主轴，采用OpenPIV互相关与pycpd仿射配准估计平移、二维散度、涡量和两项剪切。所有速率使用实际PTS，有限变换经矩阵对数换算；可见点均值变化不直接当作材料平移，二维收缩也不当作已证明的三维聚集。

独立合成种子和不同源图块的18类已知视觉运动中，两遍PIV在216个相关时间窗对比里输出全六参数151次，其中147次有全局仿射真值，21次至少一项超过预设误差容差；12个局部反向运动病例仍误接受4个。CPD输出117次，其中111次有全局真值、12次超容差，并误接受6个局部反向病例。清晰平移、旋转、收缩和剪切可以恢复，但固定纹理、漏检和置信区间覆盖仍未通过通用科学验收。官方SEA-RAFT-S实际完成24组双帧对照；共同的18组合成输入中接受17组、6组超容差，而PIV接受14组、0组超容差，因此没有改用AI主线。没有训练模型，也没有把这些相关窗口称为独立实验。

本次原Windows接收器仍为0帧，LIVE_CAMERA_PASS=false。真实录像因果回放中，两遍PIV软件处理p95约148.31 ms，100 ms新鲜度门槛下有效六参数更新率为0；单遍候选增强p95约48.02 ms，120帧处理104帧，采集活动时段的新鲜六参数约17.51 Hz，含启动解码时段约13.16 Hz。它们是软件质量与新鲜度结果，不代表真实直播、绝对六参数准确或控制器端到端延迟。接口已支持逐字段null、valid与原因、有限队列和断流过期，并输出未来动作字段，但未执行任何PZT/AFG、LED或世界模型训练。

本轮还只读使用了F10S已有32泡、完整3 ms的真实保存轨迹，用于二维投影和数学接口验证，没有新跑声学仿真或把毫秒轨迹拼成长时间物理视频。五组英文PNG/PDF、两段按原时间播放并完整解码的MP4、全部失败、SHA和23项中文结论位于external_connect/vision_tracking_research/collective_motion_v1。原标注冻结数据和历史结果哈希保持不变，training_ready=false。下一步先解决固定背景错误接受、独立显微运动真值与不确定性校准，再恢复相机、校准像素/时间并验证真实接流；不能仅凭本轮软件完成就启动自动PZT控制。

'''
    marker='---\n\n### 尾声：故事里还没讲完的部分'
    assert orig.count(marker)==1
    updated=orig.replace(marker,addition+marker,1)
    assert updated.replace(addition,'',1)==orig
    assert re.findall(r'^#+ .*$',orig,re.M)==re.findall(r'^#+ .*$',updated,re.M)
    progress.write_bytes(updated.encode('utf-8'));(ev/'after.md').write_bytes(progress.read_bytes())
    (ev/'diff.patch').write_text(''.join(difflib.unified_diff(orig.splitlines(True),updated.splitlines(True),fromfile='before',tofile='after')))
    dump(ev/'verification.json',{'PASS':True,'backup':str(backup),'before_sha256':sha(backup),'after_sha256':sha(progress),'headings_unchanged':True,'old_text_bytewise_preserved_when_addition_removed':True})

visual={'figures':['Fig1_collective_distribution','Fig2_motion_field_six_parameters','Fig3_ground_truth_recovery','Fig4_window_tradeoff','Fig5_runtime_and_admission'],
        'inspection':'All five PNGs inspected with view_image; readable English labels, dark axes on white, aspect-preserving image panels. Video1 first/middle/last inspected; Video2 individual motion and invalid-global-model frames inspected. No hidden failed cases. Non-affine estimates displayed as algorithm failure evidence, not true global motion.',
        'video1_frames_inspected':[0,40,119],'video2_screenshots_inspected':['rotation_cw','contraction','fixed_background','local_opposite'],
        'full_decode_record':'evidence/video_decode_verification.json','decoded_visual_contact_sheet':'evidence/decoded_video_contact_sheet.png','decoded_samples':'evidence/decoded_visual_samples.json','data_altered_for_appearance':False}
for stem in visual['figures']:
    for ext in ['png','pdf']:assert (ROOT/'figures'/f'{stem}.{ext}').is_file()
dump(ROOT/'evidence/visual_inspection.json',visual)

inputs={'protected':before,'source_video':{'path':str(VIDEO),'sha256':sha(VIDEO)},'AI_weights':json.loads((ROOT/'results/ai_comparison.json').read_text())['model_files_sha256']}
outputs={}
for directory in ['src','scripts','configs','results','figures','videos','evidence','logs']:
    for p in sorted((ROOT/directory).rglob('*')):
        if p.is_file() and '__pycache__' not in str(p) and p.name not in ['delivery_git.json','finalize.log']:outputs[str(p.relative_to(ROOT))]=sha(p)
for p in sorted(ROOT.glob('*.md')):outputs[p.name]=sha(p)
for name in ['FINAL_GATE.json','requirements-lock.txt','.gitignore','.gitattributes']:outputs[name]=sha(ROOT/name)
dump(ROOT/'INPUT_OUTPUT_SHA256.json',{'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'execution':identity(),'inputs':inputs,'outputs':outputs,'excluded':'manifest itself, delivery Git stamp, finalize.log (written after manifest), virtual environments, downloaded repository and model binary files (weights separately identified above)'})
print('protected inputs unchanged:',len(before),'schema states:',count,'output files hashed:',len(outputs),'progress updated:',tag in progress.read_text())
