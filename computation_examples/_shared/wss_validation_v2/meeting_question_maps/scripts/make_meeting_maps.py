"""Read-only WSS location maps. No solver import, job control, or CFD execution.

All writes are confined to this script's grandparent meeting_question_maps/.
Run with the project .venv, Python -B, and MPLCONFIGDIR inside this folder.
"""
from pathlib import Path
import os, sys, csv, json, hashlib, argparse, xml.etree.ElementTree as ET
os.environ.setdefault('LP_NUM_THREADS', '2')
os.environ.setdefault('VTK_SMP_MAX_THREADS', '2')
sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parents[1]
V = OUT.parent
A = V.parent / 'wss_audit'
os.environ.setdefault('MPLCONFIGDIR', str(OUT / '.mplconfig'))
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Ellipse
from matplotlib.colors import Normalize

FONT = Path('/mnt/c/Windows/Fonts/msyh.ttc')
font_manager.fontManager.addfont(str(FONT))
font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/msyhbd.ttc')
plt.rcParams.update({'font.family': font_manager.FontProperties(fname=FONT).get_name(),
                     'font.size': 12, 'axes.unicode_minus': False,
                     'pdf.fonttype': 42, 'axes.spines.top': False,
                     'axes.spines.right': False})
COL = {'A': '#b74425', 'B': '#246a9b', 'C': '#815aa0', 'neutral': '#555c61'}
BG = '#ffffff'
CMAP = 'viridis'
J1 = np.array([92., 49., 111.])
J2 = np.array([130.04, 82.04, 87.18])
FRONT = np.array([1., -2.2, 1.1])
CASES = ['vessel_baseline', 'vessel_medium']
CASE_ZH = dict(zip(CASES, ['原网格', '中档网格']))
DISCLAIMER = '会议原图与当前数据版本尚未核实对应；标记用于疑问位置对照，不是已确认错误位置。'
SOURCES = {}
VIEWS = {}
ANNOTATIONS = []

def source(path):
    path = Path(path).resolve()
    if str(path) not in SOURCES:
        SOURCES[str(path)] = {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size}
    return path

def rows(path):
    return list(csv.DictReader(source(path).open()))

def dump(path, value):
    path = OUT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item()) + '\n')

def csvout(name, rr):
    path = OUT / 'data' / name
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rr[0]))
        w.writeheader(); w.writerows(rr)

def path_project(points, path):
    xyz = path['xyz']; radii = path['radii']
    arc = path['arc']; best = np.full(len(points), np.inf)
    s = np.zeros(len(points)); r = np.zeros(len(points))
    for i, (p, v, L) in enumerate(zip(xyz[:-1], np.diff(xyz, axis=0), np.diff(arc))):
        f = np.clip((points-p) @ v / L**2, 0, 1)
        d = np.linalg.norm(points-p-f[:, None]*v, axis=1)
        ok = d < best
        best[ok] = d[ok]; s[ok] = arc[i]+f[ok]*L
        r[ok] = (1-f[ok])*radii[i] + f[ok]*radii[i+1]
    return s, best, r

def path_at(path, s):
    s = np.clip(s, 0, path['arc'][-1])
    return np.array([np.interp(s, path['arc'], path['xyz'][:, k]) for k in range(3)]).T

def get_paths(ports):
    graph = np.load(source(A/'inputs/network/analysis_A_H0_graph_si.npz'))
    ids, xyz, radii = graph['ids'], graph['xyz_m']*1e6, graph['radius_m']*1e6
    edges = graph['edges'][graph['roi_internal_edge_mask']]
    adj = {}
    for i, j in edges:
        adj.setdefault(int(i), []).append(int(j)); adj.setdefault(int(j), []).append(int(i))
    def chain(a, b):
        start, stop = [int(np.flatnonzero(ids == x)[0]) for x in [a, b]]
        q = [start]; parent = {start: None}
        for i in q:
            for j in adj[i]:
                if j not in parent: parent[j] = i; q.append(j)
        path = [stop]
        while path[-1] != start: path.append(parent[path[-1]])
        return path[::-1]
    specs = [('INLET_to_J1', ports['INLET']['network_node_id'], 3238, 'INLET'),
             ('J1_to_O2', 3238, ports['O2']['network_node_id'], 'O2'),
             ('J1_to_J2', 3238, 3274, None),
             ('J2_to_O1', 3274, ports['O1']['network_node_id'], 'O1'),
             ('J2_to_O3', 3274, ports['O3']['network_node_id'], 'O3')]
    result = {}; rr = []
    for name, a, b, port in specs:
        ix = chain(a, b); xx = xyz[ix]; rad = radii[ix]; nodeids = ids[ix].tolist()
        if port:
            p = ports[port]; cap = np.array(p['fem_cap']['centroid_um'])
            caprad = p['fem_cap']['equivalent_radius_um']
            if port == 'INLET':
                xx = np.vstack([cap, xx]); rad = np.r_[caprad, rad]; nodeids = ['INLET_cap'] + nodeids
            else:
                xx = np.vstack([xx, cap]); rad = np.r_[rad, caprad]; nodeids += [port+'_cap']
        arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xx, axis=0), axis=1))]
        result[name] = dict(xyz=xx, radii=rad, arc=arc, node_ids=nodeids)
        for s, x, r, node in zip(arc, xx, rad, nodeids):
            rr.append(dict(path=name, node_id=node, s_um=s, x_um=x[0], y_um=x[1], z_um=x[2], radius_um=r))
    assert np.allclose(result['J1_to_J2']['xyz'][[0, -1]], [J1, J2], atol=1e-9)
    csvout('reference_paths_um.csv', rr)
    return result

def stats(mesh, mask, case, region):
    idx = np.flatnonzero(mask); w = mesh.cell_data['WSS_raw_Pa'][idx]
    area = mesh.cell_data['Area_m2'][idx]*1e12
    cent = mesh.cell_centers().points[idx]
    order = np.argsort(w)
    q = np.interp([.05, .5, .95], (np.cumsum(area[order])-.5*area[order])/area.sum(), w[order])
    mn, mx = np.argmin(w), np.argmax(w)
    out = dict(case=case, region=region, facets=len(idx), area_um2=float(area.sum()),
               mean_Pa=float(np.average(w, weights=area)), p05_Pa=q[0], p50_Pa=q[1], p95_Pa=q[2],
               min_Pa=float(w[mn]), max_Pa=float(w[mx]),
               min_cell_index=int(idx[mn]), max_cell_index=int(idx[mx]),
               min_boundary_facet=int(mesh.cell_data['Global_boundary_facet_zero_based'][idx[mn]]),
               max_boundary_facet=int(mesh.cell_data['Global_boundary_facet_zero_based'][idx[mx]]))
    for label, i in [('min', mn), ('max', mx)]:
        for axis, val in zip('xyz', cent[i]): out[label+'_'+axis+'_um'] = val
    for threshold in [1, 2]:
        out[f'below{threshold}_count'] = int((w < threshold).sum())
        out[f'below{threshold}_area_pct'] = float(100*area[w < threshold].sum()/area.sum())
    out['zero_count'] = int(np.sum(w == 0))
    return out

def xyz_of(stat, which):
    return np.array([stat[f'{which}_{k}_um'] for k in 'xyz'])

def prepare():
    source(V/'WSS_VALIDATION_REPORT_V2.md')
    source(V/'scripts/make_figures.py'); source(V/'scripts/vessel_regions.py')
    ports = {p['name']: p for p in json.loads(source(A/'inputs/network/roi_ports_in_a.json').read_text())['ports']}
    identities = rows(V/'data/port_identity.csv')
    portpos = {}; checks = []
    xml=ET.parse(source(V/'stage3/vessel_medium/run/solver.xml')).getroot()
    for short, long in [('INLET', 'INLET'), ('O1', 'OUTLET_01'), ('O2', 'OUTLET_02'), ('O3', 'OUTLET_03')]:
        rec = next(r for r in identities if r['case'] == 'vessel_medium' and r['boundary'] == long)
        assert Path(xml.find(f'.//Add_face[@name="{long}"]/Face_file_path').text).name == long+'.vtp'
        assert float(xml.find(f'.//Add_BC[@name="{long}"]/Value').text)==float(rec['XML_value'])
        cap = pv.read(source(V/'stage3/vessel_medium/SV_MESH/mesh-surfaces'/f'{long}.vtp'))
        tri = cap.points[cap.faces.reshape(-1, 4)[:, 1:]]
        area = .5*np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)
        center = np.average(tri.mean(axis=1), weights=area, axis=0)*1e6
        ref = np.array([float(rec[f'center_{k}_um']) for k in 'xyz'])
        assert np.linalg.norm(center-ref) < 1e-8
        portpos[short] = center
        checks.append(dict(label=short, named_surface=long+'.vtp', facet_tag=int(rec['facet_tag']),
                           x_um=center[0], y_um=center[1], z_um=center[2],
                           identity_csv_difference_um=np.linalg.norm(center-ref)))
    csvout('verified_port_coordinates.csv', checks)
    paths = get_paths(ports)
    meshes = {}; summaries = []; profiles = []; assignments = {}; branchmeta = {}
    existing = rows(V/'data/vessel_region_wss.csv')
    validation = []
    trunk = paths['J1_to_J2']; L = trunk['arc'][-1]
    # Bins fixed in physical path coordinates, chosen without consulting WSS.
    bin_edges = np.r_[np.arange(5., L-5., 2.), L-5.]
    for case in CASES:
        receipt=json.loads(source(V/'stage3'/case/'wss/COMPUTE_VALIDATION.json').read_text())
        wallfile=source(V/'stage3'/case/'wss/data/wall_wss_si.vtp')
        assert SOURCES[str(wallfile)]['sha256']==receipt['outputs_sha256']['data/wall_wss_si.vtp']
        m = pv.read(wallfile)
        assert m.bounds[1] < .001  # original coordinates are SI metres
        assert np.isfinite(m.cell_data['WSS_raw_Pa']).all()
        m.points = m.points*1e6
        cent = m.cell_centers().points
        tri = m.points[m.faces.reshape(-1, 4)[:, 1:]]
        area = .5*np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)
        assert np.allclose(area, m.cell_data['Area_m2']*1e12, rtol=1e-10, atol=1e-13)
        masks = {'all_wall': np.ones(m.n_cells, bool),
                 'J1_r5um': np.linalg.norm(cent-J1, axis=1) < 5,
                 'J2_r5um': np.linalg.norm(cent-J2, axis=1) < 5}
        for name in ['O1', 'O2', 'O3']:
            p = ports[name]; origin = np.array(p['real_cut_xyz_um']); n = np.array(p['outward_normal'])
            s = (cent-origin)@n
            r = np.linalg.norm(cent-origin-s[:, None]*n, axis=1)
            masks[name+'_extension'] = (s > 0) & (s < p['fem_cap']['real_cut_to_cap_axial_um']) & (r < 3*p['radius_um'])
        for reg, mask in masks.items():
            out = stats(m, mask, case, reg); summaries.append(out)
            old = next(r for r in existing if r['case'] == case and r['region'] == reg)
            err = max(abs(out[k]-float(old[k])) for k in ['mean_Pa','p05_Pa','p50_Pa','p95_Pa','min_Pa','max_Pa','below1_area_pct','below2_area_pct'])
            assert err < 1e-9
            validation.append(dict(case=case, region=reg, maximum_difference_from_existing_summary=err))
        projs = {name: path_project(cent, path) for name, path in paths.items()}
        names = list(paths)
        dd = np.array([projs[name][1] for name in names])
        assign = np.argmin(dd, axis=0)
        assignments[case] = dict(index=assign, names=names, projections=projs)
        s, d, r = projs['J1_to_J2']
        # Require nearest among all five branches; omit both 5 um junction balls.
        trunk_mask = (assign == names.index('J1_to_J2')) & (d < 2.5*r) & ~masks['J1_r5um'] & ~masks['J2_r5um']
        branchmeta[case] = dict(trunk_length_um=L, path_type='existing ROI SWC graph unique connected path 3238→3274',
            bin_edges_um=bin_edges, nearest_branch_competition=True, radius_gate='d < 2.5 interpolated SWC radius',
            excluded_junction_ball_radius_um=5., selected_trunk_facets=int(trunk_mask.sum()),
            selected_distance_to_path_max_um=float(d[trunk_mask].max()), selected_distance_over_radius_p95=float(np.quantile(d[trunk_mask]/r[trunk_mask], .95)))
        for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
            mask = trunk_mask & (s >= lo) & (s < hi)
            assert mask.sum() > 0
            rec = stats(m, mask, case, 'J1_to_J2_bin')
            rec.update(s_start_um=lo, s_end_um=hi, s_mid_um=.5*(lo+hi))
            profiles.append(rec)
        meshes[case] = m
    csvout('region_summary.csv', summaries)
    csvout('J1_J2_wall_profile.csv', profiles)
    csvout('existing_statistics_crosscheck.csv', validation)
    dump('data/path_statistics_definition.json', branchmeta)
    flow = rows(V/'data/vessel_flow_split_comparison.csv')
    o1 = next(r for r in flow if r['boundary'] == 'OUTLET_01')
    receipt=json.loads((V/'stage3/vessel_medium/wss/COMPUTE_VALIDATION.json').read_text())
    f=receipt['measurements']['outward_flows_by_tag_m3_s']
    assert abs(100*f['3']/(-f['4'])-float(o1['medium_fraction_pct']))<1e-10
    return meshes, summaries, profiles, portpos, ports, paths, assignments, o1

def render(name, mesh, center=None, direction=FRONT, up=(0,0,1), scale=None, size=(1400,1200), edges=False, limits=(0,55)):
    direction = np.array(direction, float); direction /= np.linalg.norm(direction)
    up = np.array(up, float); right = np.cross(up, direction); right /= np.linalg.norm(right)
    vup = np.cross(direction, right)
    if center is None:
        coords = mesh.points @ np.array([right, vup, direction]).T
        center = .5*(coords.min(axis=0)+coords.max(axis=0)) @ np.array([right, vup, direction])
    center = np.asarray(center)
    if scale is None:
        pts = (mesh.points-center) @ np.array([right,vup]).T
        scale = max(np.max(abs(pts[:,1])), np.max(abs(pts[:,0]))*size[1]/size[0])*1.18
    p = pv.Plotter(off_screen=True, window_size=size)
    p.set_background('white')
    p.add_mesh(mesh, scalars='WSS_raw_Pa', preference='cell', cmap=CMAP, clim=limits,
               lighting=False, smooth_shading=False, show_edges=edges, edge_color='#333b40',
               line_width=.5, show_scalar_bar=False, interpolate_before_map=False)
    p.camera.position = center+direction*400
    p.camera.focal_point = center; p.camera.up = vup
    p.enable_parallel_projection(); p.camera.parallel_scale = scale
    im = p.screenshot(OUT/'assets'/f'{name}.png'); p.close()
    VIEWS[name] = dict(center_um=center, direction=direction, up=vup, right=right,
                        parallel_scale_um=scale, pixel_size=size, limits_Pa=limits,
                        field='WSS_raw_Pa', association='cell', smoothing=False, lighting=False,
                        interpolation='none; constant value per original wall triangle')
    def project(xyz):
        pts = np.atleast_2d(xyz)-center
        xy = np.column_stack([size[0]/2+(pts@right)*size[1]/(2*scale),
                              size[1]/2-(pts@vup)*size[1]/(2*scale)])
        return xy[0] if np.ndim(xyz)==1 else xy
    return im, project

def panel(fig, rect, image):
    ax = fig.add_axes(rect); ax.imshow(image); ax.set_axis_off()
    return ax

def mark(ax, project, xyz, text, pos, color='#37434d', fontsize=12, align='left', marker='o', view='', occluded=False):
    xy = project(np.asarray(xyz))
    ax.plot(*xy, marker=marker, ms=6, mfc='none', mec=color, mew=1.5, zorder=8)
    ax.annotate(text, xy, xycoords='data', xytext=pos, textcoords='axes fraction',
                ha=align, va='center', fontsize=fontsize, color=color, linespacing=1.45,
                bbox=dict(boxstyle='round,pad=.3', fc='white', ec='none', alpha=.93),
                arrowprops=dict(arrowstyle='-', color=color, lw=1.2, linestyle='--' if occluded else '-'), zorder=9)
    ANNOTATIONS.append(dict(view=view, text=text, target_xyz_um=np.asarray(xyz), target_projected_pixel=xy,
                            text_axes_fraction=pos, occluded_in_this_view=occluded,
                            note='projected 3D coordinate; only local extrema visibility ray-tested'))

def header(fig, num, title, subtitle):
    fig.text(.04,.963, f'{num}  {title}', fontsize=23, weight='bold', color='#172b3b')
    fig.text(.04,.93, subtitle, fontsize=12, color='#4b5b65')
    fig.text(.04,.033, DISCLAIMER, fontsize=11, color='#5a6267')

def colorbar(fig, rect, limits=(0,55), label='原始壁面面片 WSS（Pa）｜线性 0–55'):
    ca = fig.add_axes(rect)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(*limits), cmap=CMAP), cax=ca, orientation='horizontal')
    cb.set_label(label, fontsize=11)
    if limits == (0,55): cb.set_ticks([0,15,30,45,55])
    cb.ax.tick_params(labelsize=10)

def save(fig, name):
    for ext in ['png','pdf']:
        fig.savefig(OUT/f'{name}.{ext}', dpi=220, facecolor='white')
    plt.close(fig)

def figure1(mesh, summaries, portpos, paths):
    fig = plt.figure(figsize=(18,12), facecolor='white')
    header(fig, '01', '全血管疑问位置对照', '当前中档真实壁面 · A / B / C 对应三个会议疑问 · 不判定为错误')
    im, p = render('map1_global', mesh, size=(1600,1400))
    ax = panel(fig, [.02,.13,.65,.77], im)
    # Region outlines lie outside the vessel silhouette; do not alter WSS colors.
    o2 = paths['J1_to_O2']['xyz']
    xy = p(np.vstack([J1, o2])); center = .5*(xy.min(axis=0)+xy.max(axis=0))
    width, height = np.ptp(xy,axis=0)+np.array([120,180])
    ax.add_patch(Ellipse(center, width, height, fill=False, ec=COL['A'], lw=2.1))
    bb = p(paths['J1_to_J2']['xyz'])
    # Twin offset dashed guides bracket the candidate trunk in the projection.
    delta = np.gradient(bb, axis=0); perp = np.column_stack([-delta[:,1],delta[:,0]])
    perp /= np.linalg.norm(perp,axis=1)[:,None]
    for sign in [-1,1]: ax.plot(*(bb+sign*20*perp).T, '--', color=COL['B'], lw=1.5, zorder=5)
    cc = p(paths['J2_to_O1']['xyz'])
    ax.plot(cc[:,0]+24, cc[:,1], ':', color=COL['C'], lw=2, zorder=5)
    for name, pos in [('INLET',(.32,.96)),('O1',(.82,.86)),('O2',(.06,.54)),('O3',(.94,.23))]:
        mark(ax,p,portpos[name],name,pos,fontsize=14,align='center',view='map1_global')
    mark(ax,p,J1,'J1',(.24,.49),color=COL['A'],fontsize=14,view='map1_global')
    mark(ax,p,J2,'J2',(.56,.13),color=COL['B'],fontsize=14,view='map1_global')
    for label, pos in [('A',(.12,.31)), ('B',(.46,.35)), ('C',(.55,.59))]:
        ax.text(*pos,label,transform=ax.transAxes,fontsize=20,weight='bold',color=COL[label],
                bbox=dict(boxstyle='round,pad=.18',fc='white',ec=COL[label],lw=1.2))
    low = next(r for r in summaries if r['case']=='vessel_medium' and r['region']=='all_wall')
    mark(ax,p,xyz_of(low,'min'),f"◇ 当前最低值\n{low['min_Pa']:.3f} Pa",(.87,.08),
         marker='D',fontsize=11,align='center',view='map1_global')
    def card(y, label, title, body):
        fig.text(.695,y,label, fontsize=25, weight='bold',color=COL[label])
        fig.text(.735,y+.003,title,fontsize=16,weight='bold',color=COL[label])
        fig.text(.695,y-.040,body,fontsize=13,color='#26343e',va='top',linespacing=1.65)
    card(.805,'A','第一处分叉 J1 + O2 连接区',
         '会议疑问：相邻壁面高低反差为何明显？\n所述近零点尚未定位。\n“颜色像被分开”与“连接处近零”\n是同一处疑问，不拆成两个位置。')
    card(.565,'B','J1—J2 主干：候选区域',
         '会议疑问：为何中间升高、随后降低？\n具体指向待确认。\n“经过 O2 后”指经过其分叉位置，\n继续沿主干看，不是从 O2 流出再返回。')
    card(.325,'C','O1 支路',
         '会议疑问：为何这条支路整体较低？\n当前全局图中的低值支路是合理对应。')
    fig.text(.695,.166,'◇ 当前 J2 最低值单独标记',fontsize=13,weight='bold',color=COL['neutral'])
    fig.text(.695,.143,'当前数据最低值；\n并非已确认的会议指示位置。',fontsize=12,color=COL['neutral'],linespacing=1.6,va='top')
    colorbar(fig,[.10,.105,.45,.017])
    save(fig,'Figure_01_meeting_question_locations')

def figure2(meshes, summaries, paths):
    fig = plt.figure(figsize=(18,16), facecolor='white')
    header(fig,'02','J1 局部：高低反差在哪里？','两档相同视角、相同物理窗口与 0–55 Pa 色标；反向视角用于查看背面，不平滑面片值')
    directions = [FRONT, -FRONT]
    branch_targets = {'入口侧':path_at(paths['INLET_to_J1'], paths['INLET_to_J1']['arc'][-1]-6.8),
                      '通往 O2':path_at(paths['J1_to_O2'],6.8), '通往 J2':path_at(paths['J1_to_J2'],6.8)}
    for col, case in enumerate(CASES):
        m = meshes[case]; c = m.cell_centers().points
        local = m.extract_cells(np.linalg.norm(c-J1,axis=1)<8.5).extract_surface(algorithm='dataset_surface')
        rec = next(r for r in summaries if r['case']==case and r['region']=='J1_r5um')
        for row, d in enumerate(directions):
            name = f'map2_{case}_{row}'
            im,p = render(name,local,center=J1,direction=d,scale=9.5,size=(1500,1200),edges=True)
            ax = panel(fig,[.035+col*.475,.51-row*.32,.445,.31],im)
            ax.set_title(CASE_ZH[case]+('｜正向视角' if row==0 else '｜反向视角'),fontsize=15,pad=6)
            cp=p(J1); rad=5*1200/(2*9.5)
            ax.add_patch(Circle(cp,rad,fill=False,ec='#7a8188',lw=1,ls='--',alpha=.85))
            for label, xyz in branch_targets.items():
                loc=p(xyz)/[1500,1200]; loc=np.array([loc[0],1-loc[1]])
                target=.5+(loc-.5)*1.4
                target=np.clip(target,[.07,.10],[.93,.90])
                if label=='入口侧': target=np.array([.50,.88])
                mark(ax,p,xyz,label,target,fontsize=11,align='center',view=name)
            for which,pos,color,align in [('max',(.02,.96),'#ad402c','left'),('min',(.98,.04),'#345a85','right')]:
                xyz=xyz_of(rec,which); dr=np.asarray(d)/np.linalg.norm(d)
                hit,_=local.ray_trace(xyz+dr*100,xyz-dr*.1,first_point=True)
                hidden=bool(len(hit) and np.linalg.norm(hit-xyz)>.001)
                label=('统计区高值' if which=='max' else '统计区低值')+f" {rec[which+'_Pa']:.2f} Pa"
                if hidden: label+='（背面投影）'
                mark(ax,p,xyz,label,pos,color=color,fontsize=11,align=align,view=name,occluded=hidden,
                     marker='x' if hidden else 'o')
        fig.text(.055+col*.475,.859,
                 f"5 μm 球区：均值 {rec['mean_Pa']:.2f} Pa；P5–P95 {rec['p05_Pa']:.2f}–{rec['p95_Pa']:.2f} Pa\n"
                 f"<1 Pa：{rec['below1_count']} 片；<2 Pa：{rec['below2_count']} 片；统计面片数 {rec['facets']}",
                 fontsize=12,color='#243844',linespacing=1.55)
    fig.text(.055,.141,'当前 J1 统计区未复现相应阈值的近零区（<1 Pa 或 <2 Pa）。',fontsize=15,weight='bold',color='#253e52')
    fig.text(.055,.114,'显示：三角形中心距 J1 <8.5 μm；统计：距 J1 <5 μm 的球区。灰虚线为该球的投影轮廓。',fontsize=11)
    fig.text(.055,.093,'高/低值引线指向统计区极值面片中心；虚线 × 为背面投影，反向图可见该低值位置。',fontsize=11)
    colorbar(fig,[.66,.113,.27,.012])
    save(fig,'Figure_02_J1_two_mesh_views')

def figure3(mesh, summaries, profiles, portpos, ports, paths, assignment, o1flow):
    fig = plt.figure(figsize=(20,14), facecolor='white')
    header(fig,'03','主干变化与 O1 低值：把两个疑问分开看','B 是 J1—J2 主干候选区；C 是 O1 支路。曲线读取壁面 WSS，不采样中心线速度。')
    c = mesh.cell_centers().points; names=assignment['names']; assign=assignment['index']
    trunk_mask = (assign==names.index('J1_to_J2')) | (np.linalg.norm(c-J1,axis=1)<5) | (np.linalg.norm(c-J2,axis=1)<5)
    trunkmesh=mesh.extract_cells(trunk_mask).extract_surface(algorithm='dataset_surface')
    path=paths['J1_to_J2']; L=path['arc'][-1]
    # Orient the trunk horizontally in its best-fit plane; this is a display rotation only.
    right=(J2-J1)/np.linalg.norm(J2-J1)
    up=np.array([0.,0.,1.]); up=up-right*(up@right); up/=np.linalg.norm(up)
    direction=np.cross(right,up)
    im,p=render('map3_trunk',trunkmesh,direction=direction,up=up,size=(1700,780))
    ax=panel(fig,[.035,.49,.47,.39],im)
    ax.set_title('B｜J1 → J2 主干真实表面（中档）',fontsize=16,loc='left',color=COL['B'])
    xx=p(path['xyz']); guide=xx+np.array([0,42])
    ax.plot(*guide.T,'--',color=COL['B'],lw=1.2)
    for frac,label in [(0,'J1 / 0'),(.25,f'{L*.25:.1f}'),(.5,f'{L*.5:.1f}'),(.75,f'{L*.75:.1f}'),(1,f'J2 / {L:.1f}')]:
        xy=p(path_at(path,L*frac))
        ax.annotate(label+' μm',xy,xytext=(0,-31),textcoords='offset points',ha='center',fontsize=10,
                    bbox=dict(fc='white',ec='none',alpha=.9,pad=1.5),arrowprops=dict(arrowstyle='-',lw=.7,color=COL['B']))
    # Indicative arrows along the wall-adjacent guide; not trajectories.
    for frac in [.2,.5,.8]:
        a=p(path_at(path,L*frac))+[0,42]; b=p(path_at(path,L*frac+2))+[0,42]
        ax.annotate('',b,a,arrowprops=dict(arrowstyle='->',color=COL['B'],lw=1.5))
    ax.text(.01,.02,'虚线为参考路径的旁置投影，仅指示沿程方向；不是壁面流线。',transform=ax.transAxes,fontsize=10,color='#4b5b65')
    fig.text(.055,.805,'当前分箱均值：先升高、后降低，末端再回升。\n这支持 B 作为候选区；会议具体指向仍待确认。',fontsize=13,color='#345971',va='top',linespacing=1.6)
    plot=fig.add_axes([.075,.25,.395,.225])
    for case,color,label in [('vessel_baseline','#8e969d','原网格均值'),('vessel_medium','#216594','中档均值')]:
        rr=[r for r in profiles if r['case']==case]
        ss=[r['s_mid_um'] for r in rr]
        if case=='vessel_medium':
            plot.fill_between(ss,[r['p05_Pa'] for r in rr],[r['p95_Pa'] for r in rr],color=color,alpha=.18,label='中档 P5–P95（面积加权）')
        plot.plot(ss,[r['mean_Pa'] for r in rr],'-o',color=color,lw=1.5,ms=3,label=label)
    plot.set(xlim=(0,L),ylim=(0,55),xlabel='从 J1 沿已有参考中心线的弧长（μm）',ylabel='壁面 WSS（Pa）')
    plot.axvspan(0,5,color='#e9ecef',alpha=.8);plot.axvspan(L-5,L,color='#e9ecef',alpha=.8)
    plot.legend(fontsize=9,loc='upper right');plot.grid(alpha=.2)
    current=[r for r in profiles if r['case']=='vessel_medium']
    imax=int(np.argmax([r['mean_Pa'] for r in current])); imin=imax+int(np.argmin([r['mean_Pa'] for r in current[imax:]]))
    for r,offset in [(current[imax],(-40,38)),(current[imin],(-25,-36))]:
        plot.annotate(f"{r['mean_Pa']:.2f} Pa\ns ≈ {r['s_mid_um']:.0f} μm",(r['s_mid_um'],r['mean_Pa']),
                      xytext=offset,textcoords='offset points',fontsize=9,color='#244b6b',ha='center',
                      arrowprops=dict(arrowstyle='-',color='#244b6b',lw=.8))
    fig.text(.055,.164,'2 μm 分箱（最后一箱不足 2 μm）；面片按最近参考路径段归箱，面积加权。\n与其他支路竞争归属，并排除 J1/J2 半径 5 μm 球区；端部分叉区不作沿程定量解释。',fontsize=11,linespacing=1.5,va='top')
    o1mask=(assign==names.index('J2_to_O1')) | (np.linalg.norm(c-J2,axis=1)<5)
    branch=mesh.extract_cells(o1mask).extract_surface(algorithm='dataset_surface')
    # Identical camera and exact same cells for the two O1 color scales.
    for i,lim in enumerate([(0,55),(0,5)]):
        name=f'map3_O1_scale{lim[1]}'
        im,p=render(name,branch,direction=FRONT,size=(800,1300),limits=lim)
        ax=panel(fig,[.545+i*.22,.36,.205,.50],im)
        ax.set_title('C｜'+('全局色标 0–55 Pa' if i==0 else '局部色标 0–5 Pa'),fontsize=13,color=COL['C'])
        mark(ax,p,portpos['O1'],'O1',(.72,.94),align='center',fontsize=12,view=name)
        mark(ax,p,J2,'J2',(.30,.05),align='center',fontsize=12,view=name)
        cut=np.array(ports['O1']['real_cut_xyz_um'])
        mark(ax,p,cut,'人工延伸段起点',(.05,.66),fontsize=9,view=name)
        colorbar(fig,[.565+i*.22,.32,.16,.012],lim,'线性 WSS（Pa）'+('；>5 饱和为黄色' if i else ''))
    rec=next(r for r in summaries if r['case']=='vessel_medium' and r['region']=='O1_extension')
    low=next(r for r in summaries if r['case']=='vessel_medium' and r['region']=='J2_r5um')
    fig.text(.545,.245,f"O1 流量占入口 {float(o1flow['medium_fraction_pct']):.2f}%（原网格 {float(o1flow['baseline_fraction_pct']):.2f}%）",fontsize=14,weight='bold',color='#303b46')
    fig.text(.545,.208,f"O1 人工延伸段：均值 {rec['mean_Pa']:.3f} Pa\nP5–P95：{rec['p05_Pa']:.3f}–{rec['p95_Pa']:.3f} Pa；不是整条支路均值。",fontsize=12,linespacing=1.6)
    fig.text(.545,.145,f"全局色标下的深蓝紫色 ≠ 恰好为零。\n中档所有壁面无 0 Pa 面片；J2 最低值 {low['min_Pa']:.3f} Pa 另属当前位置。",fontsize=12,linespacing=1.6)
    colorbar(fig,[.10,.091,.35,.013])
    save(fig,'Figure_03_trunk_and_O1')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--analyze-only',action='store_true');args=ap.parse_args()
    for sub in ['assets','data','logs']: (OUT/sub).mkdir(exist_ok=True)
    meshes,summaries,profiles,portpos,ports,paths,assignments,o1 = prepare()
    if not args.analyze_only:
        figure1(meshes['vessel_medium'],summaries,portpos,paths)
        figure2(meshes,summaries,paths)
        figure3(meshes['vessel_medium'],summaries,profiles,portpos,ports,paths,assignments['vessel_medium'],o1)
    dump('data/view_metadata.json',VIEWS);dump('data/annotation_targets.json',ANNOTATIONS)
    # Hash evidence before/after; all referenced source inputs must remain unchanged.
    for filename,rec in SOURCES.items():
        assert hashlib.sha256(Path(filename).read_bytes()).hexdigest()==rec['sha256'],filename
    dump('data/source_manifest.json',dict(inputs=SOURCES,all_referenced_inputs_unchanged=True,
        pyvista=pv.__version__,vtk=pv.vtk_version_info,numpy=np.__version__,matplotlib=matplotlib.__version__,
        CFD_started=False,originals_modified=False,all_writes_under=str(OUT)))
    print('Completed read-only data checks and meeting maps:',OUT)
    for case in CASES:
        print(next(r for r in summaries if r['case']==case and r['region']=='J1_r5um'))
    rr=[r for r in profiles if r['case']=='vessel_medium']
    print('TRUNK PROFILE',[(round(r['s_mid_um'],2),round(r['mean_Pa'],3)) for r in rr])

if __name__=='__main__': main()
