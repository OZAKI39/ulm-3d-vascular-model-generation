"""Analyze existing outputs and render an offline review. No launch/deploy imports."""
from common import *
import csv, math, html
import numpy as np
from rbc_checks import scan_frames, probe_checks, preparation, outcome
from rbc_h5 import read_xmf
from cloud_geometry import read_off
from rbc_geometry import membership, checked_mesh, fluid_summary, local_flow_summary


def handoff_audit(work, spec):
    initial=work/'native_shear/outer00000.xmf'
    if not initial.exists():return dict(status='NOT_REACHED',seamless_full_state=False)
    states=[read_xmf(work/f'native_shear/{side}00000.xmf') for side in ['outer','inner']]
    ids=np.concatenate([s['id'] for s in states]);order=np.argsort(ids)
    position=np.concatenate([s['position'] for s in states])[order];velocity=np.concatenate([s['velocity'] for s in states])[order]
    members=np.concatenate([np.full(len(s['id']),i) for i,s in enumerate(states)])[order]
    with np.load(work/'handoff_saved_arrays.npz',allow_pickle=False) as old:
        if len(ids)!=len(old['fluid_ids']) or len(np.unique(ids))!=len(ids):return dict(status='FAILED',reason='HANDOFF_PARTICLE_COUNT_OR_ID_ERROR')
        pe=float(np.max(np.abs(position-old['fluid_positions'])));ve=float(np.max(np.abs(velocity-old['fluid_velocities'])))
        changed=members!=old['fluid_membership']
        result=dict(status='PASS_COORDINATE_VELOCITY_MAPPING' if max(pe,ve)<=1e-5 else 'FAILED',old_ids_sha256=hashlib.sha256(old['fluid_ids'].tobytes()).hexdigest(),
                    new_ids_sha256=hashlib.sha256(ids[order].tobytes()).hexdigest(),particles=len(ids),position_max_error=pe,velocity_max_error=ve,
                    classification_changes=int(changed.sum()),changed_old_ids=old['fluid_ids'][changed].tolist(),seamless_full_state=False,
                    regenerated='IDs, native initial classification, oldPositions, forces/task graph and stochastic generator; not a full-state restart')
        _,faces=read_off(spec['mesh']);L=spec['config']['geometry']['periodic_length'];v,_=checked_mesh(old['rbc_ids'],old['rbc_positions'],faces,L);points=old['fluid_positions'][changed].copy()
        points[:,:2]+=L*np.rint((v.mean(0)[:2]-points[:,:2])/L)
        if len(points):
            before=membership(points,old['fluid_membership'][changed]==1,v,faces,spec['config']['repair']['containment_band'])
            after=membership(points,members[changed]==1,v,faces,spec['config']['repair']['containment_band'])
            result['changed_membership_audit']=dict(checked=len(points),confirmed_before=int(before['confirmed_mismatch'].sum()),confirmed_after=int(after['confirmed_mismatch'].sum()),
                                                    near_surface=int(before['uncertain'].sum()),ambiguous=int(before['ambiguous'].sum()),scope='all IDs changed by the stage-initial classifier; prepared geometry and positions before reset')
            if before['confirmed_mismatch'].any() or after['confirmed_mismatch'].any():result['status']='FAILED'
        else:result['changed_membership_audit']=dict(checked=0,confirmed_before=0,confirmed_after=0,near_surface=0,ambiguous=0)
    result['membrane']=read(work/'handoff_membrane_initialization.json')
    return result


def analyze(root,spec,execution,work=None):
    root=Path(root);work=Path(work) if work is not None else root/'simulation';start=time.monotonic();frames=list(scan_frames(work,spec,final=True).values())
    frames.sort(key=lambda x:(x['phase']=='shear',x['phase_step'],x.get('source','')))
    _,faces=read_off(spec['mesh']);write(root/'frames.json',dict(frames=frames,faces=faces.tolist(),animation='saved frames only; no interpolation'))
    completion=read(work/'completion.json') if (work/'completion.json').exists() else None
    result=outcome(spec,completion,execution.get('exit_code'),frames)
    result['actual_successful_returned_prep_steps']=0;result['actual_successful_returned_shear_steps']=0
    result['native_return_records']={}
    for phase,key in [('relaxation','actual_successful_returned_prep_steps'),('shear','actual_successful_returned_shear_steps')]:
        p=work/f'phase_{phase}.json'
        if p.exists():
            phase_record=read(p)
            result['native_return_records'][phase]=phase_record
            result[key]=phase_record.get('successful_returned_steps') if phase_record.get('native_returned') else None
    result['preparation']=read(work/'preparation_gate.json') if (work/'preparation_gate.json').exists() else dict(status='NOT_COMPLETED',passed=False)
    try:members=probe_checks(work,spec,frames)
    except Exception as ex:members=dict(strict_impermeability='NOT_VERIFIED',error=str(ex),checks=[],confirmed_mismatches=0,tested_point_frames=0)
    result['membership']=members
    try:result['handoff']=handoff_audit(work,spec)
    except Exception as ex:result['handoff']=dict(status='INCONCLUSIVE',error=str(ex))
    collision=[];overflow=[]
    for path in work.glob('bounce_counts_*.csv'):
        with path.open() as f:
            for row in csv.DictReader(f):
                if None in row.values():continue
                collision.append(row)
                if int(row['coarse_count'])>int(row['coarse_capacity']) or int(row['fine_count'])>int(row['fine_capacity']):overflow.append(row)
    errors=[]
    for path in sorted(work.glob('*.log')):
        for line in path.read_text(errors='replace').splitlines():
            if any(s in line.lower() for s in ['overflow','too many','illegal memory','cuda error','mpi_abort','fatal error','segmentation']):errors.append(dict(file=path.name,line=line[:2000]))
    result['collision']=dict(sampled_counter_rows=len(collision),overflow_rows=overflow,matching_error_lines=errors,
                             max_coarse=max([int(r['coarse_count']) for r in collision],default=None),max_fine=max([int(r['fine_count']) for r in collision],default=None),
                             capacity_status='existing fixed capacity; no resizing implementation or modification',continuous_collision_correctness='NOT_VERIFIED')
    valid=[f for f in frames if 'geometry' in f]
    result['geometry']=dict(saved_frames=len(valid),unreadable_frames=[f for f in frames if 'unreadable' in f],hard_failures=[f for f in frames if f.get('hard_failures')],
                            maximum_area_relative_drift=max([f['area_relative_drift'] for f in valid],default=None),maximum_volume_relative_drift=max([f['volume_relative_drift'] for f in valid],default=None),
                            confirmed_intersection_frames=sum(f['intersections']['confirmed_count']>0 for f in valid),uncertain_contact_frames=sum(f['intersections']['uncertain_count']>0 for f in valid),
                            scope='saved frames only; all pairs sharing any vertex excluded; X/Y periodic unwrapping with edge consistency; wall Z not periodic')
    with (root/'geometry.csv').open('w') as out:
        writer=csv.writer(out);writer.writerow(['phase','phase_step','time_star','strain','area','volume','area_relative_drift','volume_relative_drift','D','theta_deg','cx','cy','cz','max_membrane_speed','intersections'])
        for f in valid:
            g=f['geometry'];writer.writerow([f['phase'],f['phase_step'],f['time_star'],f['strain'],g['area'],g['volume'],f['area_relative_drift'],f['volume_relative_drift'],g['D'],g['theta_deg'],*g['center'],f['max_membrane_speed'],f['intersections']['confirmed_count']])
    statistics=[]
    for path in sorted(work.glob('native_stats_*.csv')):
        with path.open() as f:rows=[r for r in csv.DictReader(f) if None not in r.values()]
        statistics.append(dict(file=path.name,rows=rows,definition='native m*mean(|v|^2)/3; includes bulk flow, NOT a peculiar-velocity thermal estimator; recorded afterIntegration at native step label'))
    result['raw_native_statistics']=statistics
    result['actual_mpi_ranks']=[read(p) for p in sorted(work.glob('rank_*.json'))]
    result['timings']=read(work/'phase_timings.json') if (work/'phase_timings.json').exists() else {'status':'UNMEASURED'}
    result['timing_accounting']='setup includes wall preparation; native HDF5/trace output included in solver stage times; endpoint output and state capture/teardown separate; live observer runs concurrently and cannot be added to the critical-path wall clock'
    resource_file=work.parent/'resources.jsonl';resource_rows=[]
    if resource_file.exists():
        for line in resource_file.read_text().splitlines():
            try:resource_rows.append(json.loads(line))
            except ValueError:pass
    result['resources']=dict(samples=len(resource_rows),peak_registered_process_rss_bytes=max([sum(p['rss_bytes'] for p in r['processes']) for r in resource_rows],default=None),
                             minimum_disk_free_bytes=min([r['disk_free_bytes'] for r in resource_rows],default=None),
                             gpu_whole_device_samples=[dict(elapsed_s=r['elapsed_s'],value=r['gpu_whole_device']) for r in resource_rows],attribution='GPU measurements cover the whole device; CPU/RSS list registered process group only')
    moments=[]
    if (work/'moments.csv').exists():
        with (work/'moments.csv').open() as f:moments=list(csv.DictReader(f))
    result['endpoint_moments']=dict(rows=moments,temperature_definition='m*sum(|v-mean_velocity_in_24_z_bins|^2)/(3*(N-nonempty_bins)); target DPD kBT=1',all_steps_temperature='UNMEASURED by this peculiar-velocity estimator')
    def csv_array(path):return np.atleast_1d(np.genfromtxt(path,names=True,delimiter=',',dtype=None,encoding='utf-8'))
    result['fluid_response']=fluid_summary(work,spec['config'],csv_array)
    result['local_flow_response']=local_flow_summary(work,csv_array)
    fields=[]
    for path in sorted(work.glob('native_*/local_flow*.xmf')):
        try:
            raw=read_xmf(path);fields.append(dict(file=path.relative_to(work).as_posix(),datasets={k:dict(shape=list(v.shape),finite=bool(np.isfinite(v).all())) for k,v in raw.items()}))
        except Exception as ex:fields.append(dict(file=path.relative_to(work).as_posix(),error=str(ex)))
    result['native_local_flow_files']=fields
    hard=bool(result['geometry']['hard_failures'] or overflow or errors or members['confirmed_mismatches'] or execution.get('exit_code') not in (0,None) or result['handoff'].get('status')=='FAILED')
    soft=any(f['soft_failures'] for f in valid) or result['preparation'].get('passed') is False
    for r in moments:
        if int(r['N'])!=110592 or int(r['wall_crossings'])!=0:hard=True
        if abs(float(r['temperature'])-spec['config']['dpd']['kBT'])>spec['config']['screen']['temperature_relative_error']:soft=True
    for key in ['effective_shear_relative_error','near_wall_density_relative_error','temperature_relative_error','mass_relative_drift']:
        if result['fluid_response'] and result['fluid_response'].get(key,0)>spec['config']['screen'][key]:soft=True
    inconclusive=not valid or result['geometry']['unreadable_frames'] or not members['tested_point_frames'] or any(r['status']!='CHECKED' for r in members['checks'])
    inconclusive=inconclusive or any(r.get('second_method_ambiguous',0) for r in members['checks']) or any('error' in f or not all(v['finite'] for v in f['datasets'].values()) for f in fields)
    result['CLOUD_RBC_NUMERICAL_SCREEN']='FAILED' if hard or soft else 'INCONCLUSIVE' if inconclusive or result['CLOUD_RBC_RUN_COMPLETE']!='PASS' else 'PASS'
    result['screen_scope']='Only saved geometries and enumerated probes; not all-time/all-particle impermeability proof. Native bulk kinetic statistics do not replace thermal fluctuation measurements.'
    result['RESULTS_RETURN_VERIFIED']='PENDING_LOCAL_FETCH_VERIFY'
    result['analysis_cpu_wall_s']=time.monotonic()-start;result['browser_check']='NOT_RUN';result['human_review']='PENDING'
    write(root/'numerical_checks.json',result)
    if result['CLOUD_RBC_RUN_COMPLETE']!='PASS' or result['CLOUD_RBC_NUMERICAL_SCREEN']!='PASS':
        write(root/'failure_summary.json',dict(run_complete=result['CLOUD_RBC_RUN_COMPLETE'],numerical_screen=result['CLOUD_RBC_NUMERICAL_SCREEN'],exit_code=execution.get('exit_code'),stop_reason=execution.get('stop_reason'),
                                               returned_prep_steps=result['actual_successful_returned_prep_steps'],returned_shear_steps=result['actual_successful_returned_shear_steps'],last_saved_progress_lower_bound=result['last_saved_progress_lower_bound'],retry=False))
    lines=['# 云端单红细胞完整验证',f"是否完成 Γ=4：{result['CLOUD_RBC_RUN_COMPLETE']}",f"数值筛选：{result['CLOUD_RBC_NUMERICAL_SCREEN']}",
           f"云端库：{spec['native_library']}；SHA-256：{spec['native_library_sha256']}",f"准备：{result['preparation']['status']}；原判据未放宽。",
           f"原生正常返回步数：准备 {result['actual_successful_returned_prep_steps']}，剪切 {result['actual_successful_returned_shear_steps']}。",
           f"最后保存帧确认的应变下界：{result['strain_lower_bound']}；失败的连续调用不把计划步数记为完成。",
           f"几何完整帧 {len(valid)}；最大面积漂移 {result['geometry']['maximum_area_relative_drift']}，最大体积漂移 {result['geometry']['maximum_volume_relative_drift']}。",
           f"确证自交帧 {result['geometry']['confirmed_intersection_frames']}；几何硬错误 {len(result['geometry']['hard_failures'])}。近接触另列。",
           f"成员探针：{members['tested_point_frames']} 点次，确证不符 {members['confirmed_mismatches']}。每阶段初态及正常返回末态最多各 192 点；不证明所有粒子严格不渗透。",
           f"碰撞计数采样 {len(collision)} 条；容量溢出记录 {len(overflow)} 条；原始报错匹配 {len(errors)} 行。没有修改原生容量。",
           f"退出码 {execution.get('exit_code')}；停止原因 {execution.get('stop_reason')}；求解进程墙钟 {execution.get('solver_process_wall_s')} 秒。",
           '壁粒子准备 4000 步不计入正式 Γ；setup 已包含壁准备，连续求解耗时包含原生输出，不能重复相加。',
           'GPU 显存来自全设备采样；求解进程墙钟不是 Vast 租赁账单。',
           '静壁到动壁沿用既有两协调器交接：坐标/速度/参考网格保留；ID、分类、oldPositions、力和随机状态重新生成。不是完整状态无缝延续。',
           'DPD 温度目标 kBT=1；原生 Stats 含流动动能，末态 moments 使用扣除分箱平均流速的估计。',
           f"物理验证：{result['PHYSICAL_MODEL_VALIDATION']}；qualified_speedup=null。",'下载校验由本地 LOCAL_ARCHIVE_VERIFIED.json 确认；人工验收 PENDING。']
    (root/'report_zh.md').write_text('\n\n'.join(lines)+'\n')
    return result


def render_review(archive, destination):
    """Derived files live beside the immutable raw archive, so verify remains exact."""
    root=Path(archive);verified=validate_result(root);destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    import plotly.graph_objects as go
    import plotly.io as pio
    data_root=root
    if not all((root/name).exists() for name in ['numerical_checks.json','frames.json','report_zh.md']):
        spec=read(root/'full_spec.json');provenance=read(root/'provenance.json')
        mesh=[r for r in provenance['python_inputs'] if r['remote']==spec['mesh']]
        if len(mesh)!=1 or sha(mesh[0]['local'])!=mesh[0]['sha256']:raise ValueError('LOCAL_ANALYSIS_MESH_IDENTITY_MISMATCH')
        spec['mesh']=mesh[0]['local'];analyze(destination,spec,read(root/'execution.json'),root/'simulation');data_root=destination
    checks=read(data_root/'numerical_checks.json');raw=read(data_root/'frames.json');frames=[f for f in raw['frames'] if 'vertices' in f];faces=np.array(raw['faces'])
    report_text=(data_root/'report_zh.md').read_text()
    postmortem=None
    # Local interpretation only: the uploaded worker and raw results retain
    # their frozen hashes. This case report is never used by a solver.
    from rbc_postmortem import REPORT_JOB_ID, analyze_failures, report_text as failure_report
    if root.name==REPORT_JOB_ID:
        postmortem=analyze_failures(root,destination,checks)
        report_text=failure_report(root,checks,verified,postmortem)
        checks=dict(checks,RESULTS_RETURN_VERIFIED='PASS',postmortem=postmortem)
        checks['timings_interpretation']='relaxation_s=0 is an unclosed failed-call placeholder, not measured zero runtime'
        checks['browser_check']='SEE_LOCAL_BROWSER_EVIDENCE: browser_final/browser_check.json; raw NOT_RUN was recorded before local fetch'
        write(destination/'numerical_checks.json',checks)
        failure=dict(read(root/'failure_summary.json'),stop_reason='NATIVE_COARSE_COLLISION_OVERFLOW',
                     stop_reason_source='console.log and bounce_counts_relaxation.csv; original execution.stop_reason remains null',
                     failure_snapshots=postmortem['failure_snapshots'],raw_archive_unmodified=True)
        write(destination/'failure_summary.json',failure)
    plots=[]
    if frames:
        def mesh(f):
            v=np.array(f['vertices']);return go.Mesh3d(x=v[:,0],y=v[:,1],z=v[:,2],i=faces[:,0],j=faces[:,1],k=faces[:,2],color='#c8454c',flatshading=False,opacity=.95)
        fig=go.Figure(data=[mesh(frames[0])],frames=[go.Frame(name=str(i),data=[mesh(f)]) for i,f in enumerate(frames)])
        all_v=np.concatenate([np.array(f['vertices']) for f in frames]);ranges=[[min(0,float(all_v[:,k].min())-1),max(24,float(all_v[:,k].max())+1)] for k in range(3)]
        fig.update_layout(title='真实膜帧：准备阶段 Γ=0；未保存时刻不插值',scene=dict(xaxis=dict(range=ranges[0]),yaxis=dict(range=ranges[1]),zaxis=dict(range=ranges[2]),aspectmode='cube'),height=650,
                          updatemenus=[dict(type='buttons',buttons=[dict(label='播放保存帧',method='animate',args=[None,dict(frame=dict(duration=180,redraw=True),transition=dict(duration=0),fromcurrent=True)]),dict(label='暂停',method='animate',args=[[None],dict(mode='immediate',frame=dict(duration=0,redraw=False))])])],
                          sliders=[dict(steps=[dict(method='animate',args=[[str(i)],dict(mode='immediate',frame=dict(duration=0,redraw=True),transition=dict(duration=0))],label=f"{f['phase']} {f['phase_step']} Γ={f['strain']:.3g}") for i,f in enumerate(frames)])])
        plots.append(pio.to_html(fig,include_plotlyjs=True,full_html=False,auto_play=False))
        for phase in ['relaxation','shear']:
            fs=[f for f in frames if f['phase']==phase]
            if not fs:continue
            x=[f['time_star'] if phase=='relaxation' else f['strain'] for f in fs]
            for title,fields in [('A/V 相对变化',[('area_relative_drift','面积'),('volume_relative_drift','体积')]),('形变与倾角',[('D','形变 D'),('theta_deg','倾角 °')])]:
                plot=go.Figure()
                for key,label in fields:plot.add_scatter(x=x,y=[f[key] if key in f else f['geometry'][key] for f in fs],mode='markers',name=label,connectgaps=False,yaxis='y2' if key=='theta_deg' else 'y')
                plot.update_layout(title=phase+'：'+title,xaxis_title='准备 t*（Γ=0）' if phase=='relaxation' else 'Γ',height=330)
                if fields[0][0]=='D':plot.update_layout(yaxis_title='D',yaxis2=dict(title='倾角 °',overlaying='y',side='right'))
                plots.append(pio.to_html(plot,include_plotlyjs=False,full_html=False))
    else:plots=['<p>没有完整膜帧；动画未生成。</p>']
    details=dict(checks);details.pop('raw_native_statistics',None);details['geometry']=dict(checks['geometry']);details['geometry'].pop('hard_failures',None)
    content='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Mirheo 单红细胞核查</title><style>body{font:16px system-ui;margin:24px auto;max-width:1200px;padding:0 20px;background:#f7f8fa;color:#182332}pre{white-space:pre-wrap;background:white;padding:18px}h1{font-size:28px}</style><h1>Mirheo 单红细胞核查</h1>'
    content+=f"<p>Γ=4：<b>{checks['CLOUD_RBC_RUN_COMPLETE']}</b>；数值筛选：<b>{checks['CLOUD_RBC_NUMERICAL_SCREEN']}</b>；回传：<b>{verified['status']}</b>。人工验收 PENDING。</p>"
    if postmortem:
        content+='<p style="background:#fee9e8;border-left:5px solid #bd3037;padding:16px"><b>准备阶段严重失稳，Γ=0。</b> 粗候选 18,074 &gt; 6,400；故障快照有 7 个膜顶点越墙，周期展开失败。下方动画仅含此前 46 个完整定期帧，不能代表报错瞬间几何仍正常。故障快照的自交数、A/V 与成员分类无法可靠判定。</p>'
    content+='<p>只展示真实保存帧。成员探针范围有限；物理材料匹配与收敛未在本轮验证。</p>'+''.join(plots)
    content+='<h2>版本、进度、异常、耗时与资源</h2><pre>'+html.escape(report_text)+'</pre><pre>'+html.escape(json.dumps(dict(execution=read(root/'execution.json'),output=read(root/'disk_after.json'),checks=details),ensure_ascii=False,indent=2))+'</pre></html>'
    (destination/'rbc_full_review.html').write_text(content)
    write(destination/'review_provenance.json',dict(raw_archive=str(root),raw_manifest_sha256=verified['manifest_sha256'],html_sha256=sha(destination/'rbc_full_review.html'),generated_at=now(),browser_check='NOT_RUN',human_review='PENDING',
          local_exporter_sha256=sha(__file__),local_postmortem_sha256=sha(Path(__file__).with_name('rbc_postmortem.py')),new_solver_runs=0))
    shutil_path=destination/'report_zh.md';shutil_path.write_text(report_text+f"\n本地回传校验：{verified['status']}，{verified['files']} 文件，{verified['bytes']} 字节。\n")
    return str(destination/'rbc_full_review.html')
