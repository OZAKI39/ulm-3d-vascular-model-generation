"""Self-contained Chinese Plotly HTML; offline reads never execute a solver."""
import html
import json
import os
from pathlib import Path
import numpy as np
from plotly.offline import get_plotlyjs
from py_scripts.fluid_physics.common import PROJECT_ROOT,write_json,read_json,atomic_state,fingerprint,sha256_file,now
from py_scripts.fluid_physics.analysis import periodic_shape
from .analysis import analyze_all,historical_reference
from .workflow import paths,gpu_budget


def fmt(v,digits=6):
    if v is None:return '未测 / null'
    if isinstance(v,bool):return '通过' if v else '未通过'
    if isinstance(v,(int,float)):return f'{v:.{digits}g}'
    return html.escape(str(v))


def table(headers,rows,ident=''):
    return '<div class="table-scroll"><table'+(f' id="{ident}"' if ident else '')+'><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+fmt(v)+'</td>' for v in r)+'</tr>' for r in rows)+'</tbody></table></div>'


def build_html(data,c):
    p=data['physics'];rows=data['comparison']['results']
    main=[r for r in rows if r['role']=='main'];actual=[r for r in rows if r.get('actual_steps') is not None]
    mir=next((r for r in reversed(rows) if r['backend']=='Mirheo' and r.get('actual_steps')),None)
    lbm=[r for r in main if r['backend']=='HemoCell' and r.get('actual_steps')]
    qualified=[r for r in main if r['qualified']]
    allow=bool(mir and mir['qualified'] and any(r['qualified'] for r in lbm))
    lines=[f'实际运行：HemoCell 纯流体部署冒烟、16³ 的 1/2/4 rank 主试验与低 I/O 成本段、32³ 的 1 rank 分辨率检查；主试验独立重复两次，另补四次 CPU 内存审计。SDPD：'+(f'完成 {mir["actual_steps"]:,} 步连续试验。' if mir else '新流体步进未完成，启动失败与缺项见下。'),
           f'共同准确性初筛：{len([r for r in qualified if r["backend"]=="HemoCell"])} 条 LBM 主试验通过；SDPD '+('通过。' if mir and mir['qualified'] else '尚未通过（具体门槛及测量见下）。'),
           '合格解速度比：'+('可按固定窗口范围给出，包含准备与分辨率核查成本。' if allow else '不允许给出；qualified_speedup=null。速度、安装成功与科学验收分别判断。')]
    sections=[]
    sections.append('<section><h2>共同物理定义与设备</h2>'+table(['参数','共同目标或口径'],[
        ['ρ / ν / μ',f'{p["rho_si"]!r} kg/m³ / {p["nu_si"]!r} m²/s / {p["mu_si"]!r} Pa·s'],
        ['周期盒 / 宏观初速度',f'{p["box_si"]!r} m / [0,0,0]；x,y,z 均周期'],
        ['加速度 / 力密度',f'{p["acceleration_si"]!r} m/s² / {p["force_density_si"]!r} N/m³；下半 −x、上半 +x'],
        ['每粒子力',f'{p["particle_force_si"]!r} N；F*={p["particle_force_star"]!r}；F=m a，只转换一次'],
        ['驱动设计',c['physics']['force_rationale']],
        ['目标峰值速度 / 半盒平均速度',f'{p["expected_max_velocity_si"]*1000:.9g} / {p["expected_half_mean_si"]*1000:.9g} mm/s'],
        ['Re / SDPD 目标流动 Mach',f'{p["expected_Re_box"]:.9g} / {p["expected_SDPD_Mach"]:.9g}'],
        ['背景温度',f'{p["temperature_K"]} K；旧数值测试目标，并非已测 PDMS 实验液体。LBM 不显式模拟热噪声，热温度测量不适用。'],
        ['共同时间安排',f'0–1 µs 无驱动准备；1–3 µs 建立流场；3–8 µs 固定正式统计。黏性最慢模态时间 {p["relaxation_time_si"]*1e6:.9g} µs；这不保证 SDPD 热松弛完成。'],
        ['共同分箱 / 采样',f'{p["bins"]} 个体积 bin，宽 {p["bin_width_si"]*1e6:.9g} µm。请求间隔 {p["mapping"]["t0"]*c["sampling"]["interval_star"]*1e9:.9g} ns，各求解器取最近整数步并记录。'],
        ['WSL 实际硬件','Intel i7-13700HX；24 个可用逻辑 CPU；WSL 内存 7.62 GiB。RTX 4060 Laptop GPU，8,188 MiB 显存。OpenMPI 4.1.6；OMP_NUM_THREADS=1。'],
        ['公平性边界','CPU 双精度 LBM 对 GPU 单精度 SDPD，属于当前实际部署对照。粗格距 0.25 µm、细格距 0.125 µm；SDPD 平均间距 0.25 µm、核半径 0.5 µm，不能认为精度相同。']])+'</section>')
    design=[]
    for r in main:
        comp=r.get('completion',{});taskdt=r.get('dt_si');duration=r.get('actual_time_si')
        design.append([r['task_id'],r['backend'],r['precision'],r['rank_count'],comp.get('lattice_nodes',comp.get('fluid_particles')),
                       taskdt*1e9 if taskdt else None,r.get('actual_steps'),duration*1e6 if duration is not None else None,
                       comp.get('tau'),r.get('Mach_measured_max'),r.get('lattice_max_velocity'),'CPU' if r['backend']=='HemoCell' else 'GPU + CPU',r['qualified']])
    sections.append('<section><h2>A · 实验对照表</h2>'+table(['任务','后端','精度','rank','格点/粒子','dt (ns)','实际步数','实际物理时间 (µs)','τ','峰值速率/cs','LBM峰值格子速度','设备','共同初筛'],design,'experiments')+
        '<p>格子节点是 (i+½)dx，N 个唯一周期节点对应 Ndx 盒长。所有格点均为流体，无壁面粒子、细胞或细胞回调。Guo 接口外场按加速度 a·dt²/dx 输入，computeVelocity 已包含 +a·dt/2；原始动量速度另存用于核对。τ=½+3νdt/dx²，ρLBM=1，绝不使用 DPD 数密度充当 LBM 质量密度。</p></section>')
    sections.append('<section><h2>B · 平均速度剖面</h2><p>独立解析解：下半盒 u=−a·y(H−y)/(2ν)，上半盒符号相反，H=Ly/2；曲线不传给求解器。误差使用二次曲线的精确分箱平均。测量、拟合与目标分开显示，全部空间 bin 保留。</p><div id="profile" class="plot"></div><p>LBM 的确定性无波动不表示零总误差；粗细网格差异 '+fmt(data['comparison']['refinement_relative_l2']*100 if data['comparison']['refinement_relative_l2'] is not None else None)+'%。统计 CI 只在足够完整独立块且满足稳态条件时显示；缺测后端不画假曲线。</p></section>')
    accuracy=[]
    for r in main:
        a=r.get('accuracy') or {};accuracy.append([r['task_id'],(a.get('profile_relative_l2')*100 if a.get('profile_relative_l2') is not None else None),
            (a.get('half_flow_relative_error')*100 if a.get('half_flow_relative_error') is not None else None),
            (a.get('apparent_nu_relative_error')*100 if a.get('apparent_nu_relative_error') is not None else None),a.get('apparent_nu_si'),
            a.get('stationarity_relative_drift'),a.get('global_density_relative_error'),r['sampling_status'],r.get('temperature_status'),
            ', '.join(k for k,v in a.get('gates',{}).items() if not v) or ('无局部门槛失败' if a else '未测')])
    sections.append('<section><h2>C · 收敛、统计与状态</h2><p>PROPOSED：剖面 L2、两个半盒通量、拟合 ν 各≤5%；剖面漂移≤1%；全局密度误差≤0.1%；局部密度变化≤5%；粗细网格差≤5%。SDPD 还需有效剖面 CI≤5%，温度误差连同 CI≤2%。这些只是本轮建议筛选，不是实验最终标准，也不修改历史 SDPD 门槛。</p>'+table(['任务','剖面 L2 (%)','半盒通量误差 (%)','ν 误差 (%)','拟合 ν (m²/s)','剖面漂移','全局密度误差','统计状态','温度状态','未过门槛'],accuracy,'accuracy')+
        '<div id="flow" class="plot"></div><div id="error" class="plot"></div><div id="temperature" class="plot"></div><p>全周期净通量接近零是反向周期流对称性的结果；验收分别计算上下半盒通量，绝不以相互抵消的总通量代替流量。所有启动峰值、变号区域和不合格数据保留。温度先扣除每个瞬时局部 bin 的流速并做自由度修正，另以更细 bin 核对剪切残差。</p></section>')
    if mir:
        a=mir['accuracy'];t=mir['temperature_statistics'];s=a['block_statistics']
        sections.append('<section><h2>SDPD 统计缺口的具体数值</h2>'+table(['诊断','实测或分析结果'],[
            ['固定窗口温度均值 / 目标',f'{t["measured_mean_K"]:.9g} K / {p["temperature_K"]} K'],
            ['温度相对偏差 / 前后半窗漂移',f'{abs(t["measured_mean_K"]/p["temperature_K"]-1)*100:.6g}% / {t["drift_relative_to_target"]*100:.6g}%'],
            ['温度完整块 / 有效 CI',f'{t["block_count"]} / null；非稳态 ACF 仅描述，不证明独立性'],
            ['剖面完整块 / 使用样本',f'{s["block_count"]} / {s["estimator_sample_count"]}，至少需 {c["criteria"]["min_blocks"]} 块才允许 CI'],
            ['剖面均值实际估计时间段 (µs)',str([x*1e6 for x in a['profile_estimator_interval_si']])],
            ['全正式窗口描述性剖面 L2（另列）',f'{a["full_window_descriptive_metrics"]["profile_relative_l2"]*100:.6g}%'],
            ['局部瞬时分箱密度最大变化',f'{a["local_density_relative_variation"]*100:.6g}%；含粒子热涨落，仅为本次保守局部筛选，不证明两种方法的热统计应相同'],
            ['细化温度 bin 的差异',f'{t["bin_refinement_difference_relative"]*100:.6g}%'],
            ['表观黏度解释','由实际非稳态流场独立拟合的描述值，不是已验证的材料黏度。']])+
            '<p>剖面均值与块统计使用相同完整块；未使用的时间尾部单列。原始全窗样本和全窗描述性剖面也保留。空间点没有删减，不把相关样本数当成独立实验数，不用温度失败隐藏已经测到的流场和执行成本。</p></section>')
    costs=[]
    for r in rows:
        t=r.get('timing') or {};mem=r.get('memory') or {};host=mem.get('host',{})
        costs.append([r['task_id'],r['role'],r['rank_count'],t.get('total_wall_s'),t.get('compute_ms_per_step'),t.get('workflow_wall_s_per_us'),
            t.get('setup_s'),t.get('compute_s'),t.get('sampling_transfer_reduction_output_s'),
            host.get('sampled_peak_tree_rss_bytes')/2**20 if host.get('sampled_peak_tree_rss_bytes') else None,mem.get('device_sampled_peak_used_MiB'),r.get('time_to_qualified_solution_s')])
    sections.append('<section><h2>D · 执行成本与达到筛选要求的成本</h2><p>短段数字是执行成本，不是获得正确稳态解的耗时。表中保留全部 rank 和两次独立执行；缓存读取不算重跑。低 I/O 段仅作核心成本辅助，其采样间隔与正式工作流不同。</p>'+table(['任务','用途','rank','总墙钟 (s)','计算 ms/步','工作流 s/µs','初始化 (s)','计算 (s)','采样/输出 (s)','主机 RSS 采样峰值 (MiB)','GPU 全设备采样峰值 (MiB)','合格固定窗口+网格核查 (s)'],costs,'costs')+
        '<div id="cost" class="plot"></div><div id="memory" class="plot"></div><p>MPI 每段取所有 rank 的最大完成时间；SDPD 计时结束前执行 cudaDeviceSynchronize。下载和编译不计入 ms/步。初始化、时空采样、主机传输、MPI 归约和输出分列；细小计时日志与退出等未分类耗时保留在总墙钟中。补充 CPU 内存为每 0.1 s 对本任务完整进程树 RSS 求和，包含重复映射的共享页，不是 PSS；GPU 为设备整体采样值，含桌面/其他程序，绝不是 RSS。各 rank 样本和实际亲和性保存在 host_memory.json。</p><p>“合格固定窗口”不是最早停止时间：从全新静止状态完成预定 8 µs 窗口，再计一次另一网格核查。未通过则 null。实验整体安装/校核成本另列，不能用一次缓存核查冒充另一次实测。</p></section>')
    rankrows=[]
    for r in rows:
        mem=r.get('memory') or {};host=mem.get('host',{})
        if host.get('coverage')=='DESCENDANT_TREE':
            for rank,item in host['ranks'].items():rankrows.append([r['task_id'],rank,item['sampled_peak_rss_bytes']/2**20,item['cpus_allowed'],host['minimum_WSL_available_bytes']/2**20])
    sections.append('<section><h2>内存覆盖范围修正与各 rank 实测</h2><p>原进程组采样遗漏了 OpenMPI 自建进程组的计算 rank，只采到约 13 MiB 的 mpirun 启动器；原记录保留，但其求解器主机 RSS 已改标未测。补充 CPU 试验沿真实父子关系跨进程组采样，以下是独立重复的资源结果，不回填成原计时试验的峰值。SDPD 的 rank 主机 RSS 仍缺测；设备显存采样有效。没有因此启动第三次 GPU 试验。</p>'+table(['补充运行','MPI rank','RSS 采样峰值 (MiB)','实际 CPU 亲和性','WSL 最低可用内存 (MiB)'],rankrows,'rank-memory')+'</section>')
    h=data['historical_reference']
    sections.append('<section><h2>E · 预算、来源与限制</h2><p>HemoCell 数值预算上限 1,200 s；本轮累计实收 '+fmt(data['cpu_charged_s'])+' s。新 GPU 授权 600 s；含失败启动累计实收 '+fmt(data['gpu_charged_s'])+' s。所有历史用量、原授权和用途限制保留，共享账本按新授权追加记录。失败/中断不重置额度，不自动重试。</p>'+table(['历史参考（单独列出，不参与新对照）','数值'],[
        ['旧 SDPD sdpd_flow 真实步数 / 物理时长',f'{h["actual_steps"]} / {h["actual_time_si"]*1e6:.9g} µs'],
        ['旧工作流墙钟 / 计算成本',f'{h["total_wall_s"]:.6g} s / {h["compute_ms_per_step"]:.6g} ms/步'],
        ['旧工作流 s/µs / ν 误差',f'{h["workflow_wall_s_per_us"]:.6g} / {h["nu_relative_error"]*100:.6g}%'],
        ['旧温度均值 / CI',f'{h["temperature_mean_K"]:.9g} K / 无有效 CI'],['历史限制',h['limitations']]])+
        '<p>首次新 SDPD 启动失败发生在导入 mpi4py 阶段（未导入 Mirheo、未推进流体），实收 2.458043783 s；修正版直接使用现有 OpenMPI C ABI，通过真实双 rank CPU 检查。原失败目录和版本均保留。是否显式重试、实际成功轨迹与预算详见嵌入结果。</p>'+table(['来源/构建','记录'],[
            ['HemoCell 固定提交',c['hemocell_commit']],['Palabos 固定提交',c['palabos_commit']],
            ['编译器/构建模式','GCC/G++ 13.3.0，系统 OpenMPI wrapper；Release -O3 -march=native；HemoCell 原库保留上游调试符号；2 个编译进程。'],
            ['MPI / HDF5','OpenMPI 4.1.6；并行 HDF5 1.10.10，C 与 HL 组件已找到并链接。无系统包安装/升级。'],
            ['独立可执行文件',str(Path(c['hemocell_root'])/'build/benchmark/pure_fluid_benchmark')],
            ['HemoCell 核心库编译时间 (s)',data['installation'].get('library_build_s')],
            ['HemoCell 基准编译含首次 API 修正 (s)',data['installation'].get('case_build_s')],
            ['许可','HemoCell、Palabos 和派生 C++ 案例：AGPL-3.0-or-later；Mirheo：MIT。未修改任何上游核心求解代码。'],
            ['配置 SHA256',c['_config_sha256']],['保护核查',data.get('protection_status','最终核查待完成')],
            ['自动浏览器检查',data.get('browser_status','独立 browser_check.json 为准')],['人工验收','PENDING']])+
        '<p>官方来源：<a href="https://github.com/UvaCsl/HemoCell/tree/'+c['hemocell_commit']+'">固定 HemoCell 源码</a>；<a href="https://hemocell.eu/user_guide/QuickStart.html">部署说明</a>；<a href="https://hemocell.eu/user_guide/advanced_cases/pure_flow_simulations.html">纯流体路径</a>；<a href="https://hemocell.eu/user_guide/concepts/units_and_scaling.html">单位说明</a>；<a href="https://hemocell.eu/user_guide/cases/pipeflow.html">官方圆管案例</a>。在线文档的 v2.3.0 标签说明未覆盖固定 setup.sh 锁定的具体 Palabos 提交，实际使用后者及配套补丁。</p><p>PIPE_BENCHMARK_DEFERRED：SDPD 近壁核截断、冻结壁面密度和滑移尚未验证。本次没有 RBC、微泡、黏附、真实血管、PDMS 实验迁移或完整压力出口验收。小盒内存/速度不能线性承诺完整稀疏血管。建议先完成 SDPD 平均流场、温度、时间步和壁面验证，再考虑 RBC 对照；后续还需匹配细胞尺寸、膜力学、浓度、内外黏度及耦合精度。</p></section>')
    encoded=json.dumps(data,ensure_ascii=False,allow_nan=False,separators=(',',':')).replace('</',r'<\/')
    return '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Mirheo 与 HemoCell 纯流体核查</title><style>
body{margin:0;background:#f3f5f8;color:#1c2937;font:16px/1.65 "Segoe UI","Microsoft YaHei",sans-serif}main{max-width:1500px;margin:auto;padding:26px}h1{font-size:29px;margin:0 0 16px}h2{font-size:23px;margin:0 0 12px}section,.lead{background:white;margin:22px 0;padding:24px;border:1px solid #dbe1e7;border-radius:10px}.lead{border-left:6px solid #b07016;background:#fffaf0}.lead p{margin:8px 0}.table-scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px 11px;border:1px solid #dfe5ec;vertical-align:top;text-align:left}th{background:#eaf0f6;white-space:nowrap}tr:nth-child(even){background:#f7f9fb}.plot{height:420px;min-width:450px}a{color:#165aa7}code{word-break:break-all}button{border:1px solid #c2cbd5;border-radius:5px;background:white;padding:9px 15px;cursor:pointer}.muted{color:#536779;font-size:14px}details{padding:15px;background:#edf2f7}pre{white-space:pre-wrap;word-break:break-all}
</style></head><body><main><h1>Mirheo vs HemoCell/Palabos · 纯流体核查</h1><p class="muted">当前 WSL 电脑的实测部署对照 · 单文件离线 HTML · 非全血性能排名</p><div class="lead">'''+''.join('<p><strong>'+str(i+1)+'. </strong>'+html.escape(line)+'</p>' for i,line in enumerate(lines))+'''</div><button id="reset-plots">重置图表缩放</button>'''+''.join(sections)+'''<details><summary>完整结果、各 rank 内存样本、门槛与 provenance（离线嵌入）</summary><pre id="json-view"></pre></details><p class="muted">人工验收始终 PENDING。图形仅使用本地嵌入数据与 JavaScript，不需要 Web 服务。</p></main><script id="benchmark-audit-data" type="application/json">'''+encoded+'''</script><script>'''+get_plotlyjs()+'''</script><script>
const audit=JSON.parse(document.getElementById('benchmark-audit-data').textContent);const R=audit.comparison.results;const P=audit.physics;
const selected=R.filter(r=>r.role==='main'&&r.accuracy&&r.accuracy.measured_profile_si&&r.repetition===0&&(r.rank_count===1||r.backend==='Mirheo'));
const latestMir=R.filter(r=>r.backend==='Mirheo'&&r.accuracy&&r.accuracy.measured_profile_si).at(-1);if(latestMir&&!selected.includes(latestMir))selected.push(latestMir);
const colors=['#236aaf','#b55223','#7b3da2'];const name=r=>r.backend==='Mirheo'?'Mirheo SDPD':`${r.candidate_id} / ${r.rank_count} rank`;
const common={responsive:true,scrollZoom:true,displaylogo:false,toImageButtonOptions:{format:'png',scale:2}};
const make=(id,traces,title,x,y,extra={})=>Plotly.newPlot(id,traces,{title,xaxis:{title:x},yaxis:{title:y},margin:{t:60,l:80,r:25,b:65},legend:{orientation:'h',y:-.25},...extra},common);
const profile=[];if(selected.length){const a=selected[0].accuracy;profile.push({x:a.y_si.map(v=>v*1e6),y:a.target_profile_si.map(v=>v*1000),name:'目标解析分箱平均',mode:'lines+markers',line:{color:'#222',dash:'dot'}});}
selected.forEach((r,i)=>{const a=r.accuracy,col=colors[i%colors.length],x=a.y_si.map(v=>v*1e6),y=a.measured_profile_si.map(v=>v*1000);if(a.profile_ci95_si){const z=a.profile_ci95_si.map(v=>v*1000);profile.push({x:[...x,...x.slice().reverse()],y:[...y.map((v,j)=>v+z[j]),...y.map((v,j)=>v-z[j]).reverse()],fill:'toself',fillcolor:'rgba(100,80,160,.15)',line:{width:0},name:name(r)+' 有效95% CI',hoverinfo:'skip'});}profile.push({x,y,name:name(r)+' 实测',mode:'lines+markers',line:{color:col}});profile.push({x,y:a.fitted_profile_si.map(v=>v*1000),name:name(r)+' 数据拟合',mode:'lines',line:{color:col,dash:'dash'}});});
const plots=[make('profile',profile,'共同分箱的平均速度：目标 / 测量 / 拟合','y (µm)','uₓ (mm/s)')];
plots.push(make('flow',selected.map((r,i)=>({x:r.series.time_si.map(t=>t*1e6),y:r.series.upper_mean_si.map(v=>v*1000),name:name(r),line:{color:colors[i%3]}})),'上半盒平均速度随物理时间','t (µs)','上半盒平均 uₓ (mm/s)',{shapes:[{type:'line',x0:0,x1:8,y0:P.expected_half_mean_si*1000,y1:P.expected_half_mean_si*1000,line:{dash:'dot',color:'#444'}}]}));
plots.push(make('error',selected.map((r,i)=>({x:r.series.time_si.map(t=>t*1e6),y:r.series.profile_relative_l2.map(v=>v*100),name:name(r),line:{color:colors[i%3]}})),'瞬时剖面误差（正式验收使用固定窗口平均）','t (µs)','相对 L2 (%)'));
const therm=latestMir? [{x:latestMir.series.temperature_time_si.map(t=>t*1e6),y:latestMir.series.temperature_K,name:'SDPD 局部流速扣除后温度',line:{color:'#7b3da2'}}]:[];
plots.push(make('temperature',therm,therm.length?'SDPD 温度独立诊断（保留启动全过程）':'SDPD 新温度未测；LBM 热统计不适用','t (µs)','T (K)',{shapes:[{type:'line',x0:0,x1:8,y0:P.temperature_K,y1:P.temperature_K,line:{dash:'dot',color:'#444'}}]}));
const measured=R.filter(r=>r.timing&&r.timing.compute_ms_per_step!=null&&r.role!=='smoke');
plots.push(make('cost',[{x:measured.map(r=>r.task_id),y:measured.map(r=>r.timing.workflow_wall_s_per_us),name:'完整工作流',type:'bar',marker:{color:'#416a91'}},{x:measured.map(r=>r.task_id),y:measured.map(r=>r.timing.core_wall_s_per_us),name:'计算段',type:'bar',marker:{color:'#b8844f'}}],'相同物理时间的成本（包含全部重复）','任务','墙钟 s / 物理 µs',{barmode:'group',margin:{t:60,l:80,r:25,b:150}}));
plots.push(make('memory',[{x:measured.map(r=>r.task_id),y:measured.map(r=>r.memory?.host?.sampled_peak_tree_rss_bytes==null?null:r.memory.host.sampled_peak_tree_rss_bytes/1048576),name:'主机 RSS 采样峰值（补充审计）',type:'bar',marker:{color:'#416a91'}},{x:measured.filter(r=>r.backend==='Mirheo').map(r=>r.task_id),y:measured.filter(r=>r.backend==='Mirheo').map(r=>r.memory.device_sampled_peak_used_MiB),name:'GPU 全设备显存采样峰值',type:'bar',marker:{color:'#b8844f'}}],'主机内存与显存分开记录','任务','MiB',{barmode:'group',margin:{t:60,l:80,r:25,b:150}}));
document.getElementById('json-view').textContent=JSON.stringify(audit,null,2);document.getElementById('reset-plots').onclick=()=>document.querySelectorAll('.js-plotly-plot').forEach(e=>Plotly.relayout(e,{'xaxis.autorange':true,'yaxis.autorange':true}));
Promise.all(plots).then(()=>{window.benchmarkReady=true;});
</script></body></html>'''


def export(c,frozen):
    base,runs=paths(c);comparison=analyze_all(c,frozen);root=Path(c['hemocell_root'])
    def charges(backend):
        path=runs/backend/'budget_ledger.json'
        return sum(x.get('charged_s',x['reserved_s']) for x in read_json(path)['attempts']) if path.exists() else 0
    installation={'library_build_s':read_json(root/'metadata/build_library.json')['elapsed_monotonic_s'],
                  'case_build_s':sum(read_json(x)['elapsed_monotonic_s'] for x in (root/'metadata').glob('build_benchmark*.json')),
                  'sources':read_json(root/'metadata/sources.json')}
    data={'physics':frozen['physics'],'comparison':comparison,'criteria':c['criteria'],'config_sha256':c['_config_sha256'],
          'frozen_plan_sha256':sha256_file(base/'frozen_plan.json'),'historical_reference':historical_reference(),
          'installation':installation,'cpu_charged_s':charges('cpu'),'gpu_charged_s':charges('gpu'),
          'budget':gpu_budget(c,frozen['gpu_plan']),'human_review':'PENDING'}
    data['data_categories']={'MEASURED':'Actual solver CSV, steps, process timing and resource samples',
        'DERIVED':'SI mappings, analytic profile, fits, errors and costs computed from measured time',
        'HISTORICAL_REFERENCE':'Old SDPD evidence shown separately, never pooled with new runs',
        'ESTIMATED':'Preflight duration/cost forecasts only, not measurements',
        'NOT_RUN':'Pipe/RBC/pressure outlet and missing measurements remain null'}
    prot=root/'metadata/protection_after.json'
    if prot.exists():data['protection_status']=read_json(prot)['status']
    key=fingerprint({'data':data,'reporting_code':sha256_file(Path(__file__)),'analysis_code':sha256_file(Path(__file__).with_name('analysis.py'))})[:16]
    run_id='review_'+key;dest=PROJECT_ROOT/'test_code/outputs/solver_benchmark'/run_id;dest.mkdir(parents=True,exist_ok=True)
    result=dest/'results.json';page=dest/'benchmark_review.html'
    if not result.exists():write_json(result,data)
    elif read_json(result)!=data:raise ValueError('RESULT_COLLISION')
    if not page.exists():page.write_text(build_html(data,c),encoding='utf-8')
    text=page.read_text();encoded=text.split('<script id="benchmark-audit-data" type="application/json">',1)[1].split('</script>',1)[0]
    if json.loads(encoded)!=read_json(result):raise ValueError('HTML_JSON_READBACK_MISMATCH')
    if not (dest/'artifact_sha256.json').exists():write_json(dest/'artifact_sha256.json',{'results.json':sha256_file(result),'benchmark_review.html':sha256_file(page)})
    atomic_state(PROJECT_ROOT/'data/solver_benchmark/LATEST.json',{'run_id':run_id,'html':str(page),'results':str(result),'html_sha256':sha256_file(page)})
    return page
