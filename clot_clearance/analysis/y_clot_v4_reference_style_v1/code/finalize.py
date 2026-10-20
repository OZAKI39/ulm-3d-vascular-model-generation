"""Validate encoded artifacts and exact input preservation; write local index."""
from render import *


def main():
    assert all(json.loads((O/f'{m}_COMPLETE.json').read_text())['status']=='PASS'
               for m in ['RESULTS','ORBIT'])
    source_hashes=json.loads((O/'SOURCE_HASHES.json').read_text())
    assert all(sha(R/p)==digest for p,digest in source_hashes.items())
    ref,field,paths=inputs()
    # Explicitly quantify the partly regular reference, rather than implying a full remesh.
    old=load_npz(R/'verification/y_clot_current_method_v2/20261019T131047Z/provenance/reference.npz')
    shift=np.linalg.norm(ref['X']-old['X'],axis=1)
    n=int(ref['clot_count'])
    moved=shift>1e-14
    assert moved[:n].sum()==412 and not moved[n:].any()
    geometry=dict(material_nodes=n,perturbed_interior_material_nodes=int(moved[:n].sum()),
                  unperturbed_boundary_interface_material_nodes=int((~moved[:n]).sum()),
                  fixed_regular_anchors=len(shift)-n,
                  maximum_reference_shift_m=float(shift.max()),
                  description='Interior-jittered point cloud with retained regular boundary and anchor layers; not a fully irregular remesh',
                  boundary_material_nodes_are_dynamic=True)
    dump(O/'POINT_CLOUD_REGULARITY.json',geometry)
    dyn=json.loads((O/'DYNAMICS_FRAMES.json').read_text())
    frames=dyn['frames'];assert len(frames)==len(paths)==149
    for row,(t,path) in zip(frames,paths):
        assert row['source']==str(path.relative_to(R)) and row['sha256']==sha(path)
        assert row['time_s']==t
        s=load_npz(path);lab,attached=topology(ref,s['d'])
        assert row['free_nodes']==int((~attached[lab[:n]]).sum())
        for typ in [1,2]:
            assert row['bonds'][str(typ)]['shown_surviving_edges']==int(((s['d']<1)&(ref['bond_types']==typ)).sum())
        x=s['x']*1e3;b=np.array(dyn['bounds_mm'])
        assert (x>=b[:,0]).all() and (x<=b[:,1]).all()
    assert frames[0]['time_s']==0 and abs(frames[-1]['time_s']-.0036)<1e-12
    expected_dynamic_frames=sum(q['repeat'] for q in frames)
    videos=[]
    for path in sorted(VID.glob('*.mp4')):
        probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-count_frames',
              '-select_streams','v:0','-show_entries',
              'stream=width,height,avg_frame_rate,nb_read_frames,duration,pix_fmt',
              '-of','json',str(path)]))['streams'][0]
        expected=540 if 'ROTATION' in path.name else expected_dynamic_frames
        assert int(probe['nb_read_frames'])==expected
        assert probe['avg_frame_rate']=='30/1' and probe['pix_fmt']=='yuv420p'
        assert [probe['width'],probe['height']]==([3840,2160] if '_4K' in path.name else [1920,1080])
        subprocess.run(['ffmpeg','-v','error','-xerror','-i',str(path),'-f','null','-'],check=True)
        videos.append(dict(file=str(path.relative_to(F)),sha256=sha(path),full_decode='PASS',**probe))
    assert len(videos)==4
    pngs=[]
    for path in sorted(FIG.glob('*.png')):
        with Image.open(path) as im:
            assert im.width>=3840 and im.height>=1600
            im.verify()
        pngs.append(dict(file=str(path.relative_to(F)),sha256=sha(path)))
    assert len(pngs)==16
    assert json.loads((O/'ORBIT_LAYOUT.json').read_text())['layout']['status']=='PASS'
    load=json.loads((O/'LOAD_FIELD_VISUALIZATION.json').read_text())
    assert load['force_readback_max_N']<1e-15
    validation=dict(status='PASS',input_file_count=len(source_hashes),
                    original_data_and_physics_code_unchanged=True,
                    native_dynamics_runs_started=0,displacement_scale=1.,
                    geometry=geometry,actual_saved_states=149,dynamics_frames=expected_dynamic_frames,
                    orbit_frames=540,orbit_physical_time_s=0.,
                    phase_coloring='frozen weak/strong material phase; not component identity',
                    marker_size='14 px material / 8 px anchors in 4K; visual markers, not Vi spheres',
                    surviving_edges_only=True,full_six_face_grid=True,world_Z_vertical=True,
                    videos=videos,figures=pngs)
    dump(O/'VALIDATION.json',validation)
    report=f"""# Y clot V4：按指定视频风格重绘

本轮只读取已完成的 `y_clot_irregular_current_method_v4/20261020T041100Z` 实际数据，未启动新物理求解。{len(source_hashes)}个输入文件的 SHA-256 前后保持一致。

## 风格与科学含义

指定参考视频为 `clot_discretization_current_method_rotation_v1/20261019T205403Z/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4`。沿用其白底、弱相浅红/强相深红、蓝灰锚点、六面网格、正交低角度相机、世界Z竖直、单一图例和旋转标签渐隐算法。

4K下材料点14px、锚点8px；这些是可视化标记，**不再表示Vi等效球半径**。全部实际点均保留。图中细线来自真实当前存活本体/界面键；参考距离≤50.5µm的存活短键加深，固定支架键不显示。红色深浅表示冻结材料相区，不表示碎块编号或破坏程度。

点云仍是**保留规则边界的内部扰动点云**：1927材料点中412个内部点发生初始扰动，1515个边界/界面点保留原始排列；1561固定锚点仍规则。边界材料点在动力学中可以运动，只有锚点固定。不能把本轮图像称为完整不规则重网格结果。坐标直接取自接受态，没有为视觉效果重新布点或放大位移。

原V4任务第4.1节指定从规则点阵出发扰动内部点，允许边界点切向扰动或固定；当前Y应用采用保留边界位置的分支。未扰动边界/界面占材料点约78.6%，这是外观仍规则的主要原因。本轮参考风格又加深了真实短键，进一步突出外层规则连接；所以外观由实际离散方式与短键显示权重共同决定。

## 输出

- 初态绕Z轴360°旋转：18s、30fps、540帧，4K和1080p。只改变相机，物理时间始终为0。
- 实际碎裂视频：149个保存态、0–3.6ms、{expected_dynamic_frames}视频帧、30fps，4K和1080p。固定全程相机范围，无节点轨迹插值；首尾停留、各保存态重复显示，视频时间不是物理时间，画面标注实际物理时间。
- 初态、较大游离块态、末态、Y管定位和实际载荷场静态图，以及8个旋转视角。单张状态图按各自范围取景；动态视频保持固定范围。
- 历史曲线、能量账本和粒径/质量分布使用原始数据重新绘制，统一红色/灰蓝色、字体和网格。

载荷图直接读取已核验 `ACTUAL_LOAD_FIELD.npz`，时刻{load['state']['time_s']*1e3:.5f}ms；显示55个原始面力箭头与137个实际受载面质心。箭头有统一比例，2µN图例使用相同标尺；最大牵引{load['maximum_traction_Pa']:.3f}Pa。受力读取误差{load['force_readback_max_N']:.3e}N。没有新求解CFD、声场或声流。

## 核验与复现

4个MP4均通过全帧解码、尺寸与帧数检查；旋转540帧的Z竖直、标签碰撞/裁切/过渡检查沿用参考算法。149态逐一核对源哈希、时刻、游离节点数、存活键数量和全程显示范围。

在服务器激活 `/workspace/y_clot/env.sh` 后，运行 `python analysis/{TASK}/code/render.py results`、`python analysis/{TASK}/code/render.py orbit`、`python analysis/{TASK}/code/finalize.py`。已完成视频拒绝覆盖；再次生成应使用新输出目录。完整中间轨迹保留于服务器，原项目依赖保持。

查看 [结果索引](../../results/{TASK}/{RUN}/OPEN_RESULTS.html) 和 [校验记录](../../verification/{TASK}/{RUN}/VALIDATION.json)。
"""
    (A/'REVIEW.md').write_text(report)
    title='Y clot V4 · Reference video style'
    parts=['<!doctype html><html lang="zh"><meta charset="utf-8"><title>'+title+'</title>',
           '<style>body{max-width:1300px;margin:40px auto;padding:0 24px;color:#293a45;font:17px/1.65 Arial,sans-serif}h1{color:#88152b}a{color:#88152b}img,video{width:100%;height:auto}section{margin:40px 0}figure{margin:24px 0}small{color:#697d89}</style>',
           f'<h1>{title}</h1>',
           '<p>最新 Y 血块 V4 接受态，按指定视频的红色小节点、细键线、蓝灰锚点和六面网格重绘。</p>',
           '<p><strong>离散方式：</strong>412个内部材料点扰动；1515个边界/界面材料点和1561个锚点保留初始规则排列。并非完整不规则重网格。点标记大小仅用于显示；深浅红表示材料相区。</p>']
    for stem,label,poster in [
        ('VIDEO_CLOT_DISCRETIZATION_Z_ROTATION','初态绕 Z 轴旋转 · 18秒','figures/ROTATION_VIEW_000.png'),
        ('VIDEO_01_Y_CLOT_V4_FRAGMENTATION','真实碎裂过程 · 物理时间0–3.6ms','visual_review/DYNAMICS_074.png')]:
        parts += [f'<section><h2>{label}</h2>',
                  f'<video controls preload="metadata" poster="{poster}" src="videos/{stem}_1080P.mp4"></video>',
                  f'<p><a href="videos/{stem}_4K.mp4">4K 视频</a> · <a href="videos/{stem}_1080P.mp4">1080p 视频</a></p></section>']
    labels={
        'FIG_01_Y_VESSEL_LOCATION.png':'原 Y 形血管中的血块位置',
        'FIG_02_INITIAL_CLOT.png':'初始离散点云',
        'FIG_02_LARGE_FRAGMENT_CLOT.png':'较大游离块保存态 · 2.15212ms',
        'FIG_02_FINAL_CLOT.png':'末态 · 3.6ms',
        'FIG_03_ACTUAL_SURFACE_LOAD.png':'实际载荷场 · 1.95187ms',
        'FIG_04_LOAD_AND_FRAGMENT_HISTORY.png':'实际载荷与碎块历史',
        'FIG_05_ENERGY_LEDGER.png':'能量账本',
        'FIG_06_FRAGMENT_SIZE_DISTRIBUTION.png':'末态碎块数量与质量分布'}
    for name,label in labels.items():
        parts += [f'<figure><figcaption>{label}</figcaption><a href="figures/{name}"><img loading="lazy" src="figures/{name}"></a></figure>']
    parts += ['<h2>旋转视角</h2><a href="visual_review/ORBIT_CONTACT_SHEET.jpg"><img src="visual_review/ORBIT_CONTACT_SHEET.jpg"></a>',
              '<p>'+ ' · '.join(f'<a href="{q.relative_to(F)}">{i+1}</a>' for i,q in enumerate(sorted(FIG.glob('ROTATION_VIEW*'))))+'</p>',
              f'<p><a href="../../../analysis/{TASK}/REVIEW.md">重绘说明</a> · <a href="../../../verification/{TASK}/{RUN}/VALIDATION.json">校验记录</a></p>',
              '<small>工程构造Y血管血块；探索性、未标定的指定剪切代理计算。旋转视频为初态展示，碎裂视频才是实际时间演化。</small></html>']
    (F/'OPEN_RESULTS.html').write_text('\n'.join(parts))
    files=list(F.rglob('*'))+[A/'REVIEW.md']+list((A/'code').glob('*'))
    files += [p for p in O.glob('*.json') if p.name!='ARTIFACT_MANIFEST.json']
    manifest={str(p.relative_to(R)):dict(bytes=p.stat().st_size,sha256=sha(p))
              for p in files if p.is_file() and p.suffix!='.pyc'}
    dump(O/'ARTIFACT_MANIFEST.json',manifest)
    dump(A/'FINAL_STATUS.json',dict(status='COMPLETE',validation='PASS',task=TASK,run=RUN,
          source_trajectory='y_clot_irregular_current_method_v4/20261020T041100Z',
          scope='Visualization only',new_physics_runs=0,
          result_index=str((F/'OPEN_RESULTS.html').relative_to(R)),
          regularity=geometry))
    print('DELIVERY_PASS',len(videos),'videos',len(pngs),'figures',len(source_hashes),'unchanged inputs',flush=True)


if __name__=='__main__':
    main()
