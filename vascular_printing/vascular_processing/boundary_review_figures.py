"""Static, publication-size review figures; no vessel surface reconstruction."""
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import proj3d
import numpy as np

COLORS={'M1':'#d77d00','MeVO':'#008b91','UNKNOWN':'#8b2cb5','Current boundary':'#db2345','Context':'#aeb7c2'}
COMPONENT_COLORS=['#008b91','#d77d00','#5565bd','#ba457d','#55942c']
SIZE=(20,13.34)
DPI=120


def short(identity):return identity.split('_')[-1]


def points(cache,key):return np.array([cache.coordinates[n][:3] for n in cache.branches[key]['node_ids']])


def bounds(ax,coords):
    coords=np.asarray(coords);center=(coords.min(axis=0)+coords.max(axis=0))/2
    half=max(float(np.ptp(coords,axis=0).max())*.63,1.)
    ax.set_xlim(center[0]-half,center[0]+half);ax.set_ylim(center[1]-half,center[1]+half);ax.set_zlim(center[2]-half,center[2]+half)
    ax.set_box_aspect((1,1,1));ax.set_proj_type('ortho')
    ax.set_xlabel('X (mm)',fontsize=11,labelpad=6);ax.set_ylabel('Y (mm)',fontsize=11,labelpad=6);ax.set_zlabel('Z (mm)',fontsize=11,labelpad=6)
    ax.tick_params(labelsize=9,pad=1)
    for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.set_major_locator(plt.MaxNLocator(3))


def lines(ax,cache,keys,*,colors=None,alpha=1.,width=2.8):
    for key in keys:
        xyz=points(cache,key);color=(colors or {}).get(key,COLORS[cache.branches[key]['label']])
        ax.add_collection3d(Line3DCollection([xyz],colors=color,linewidths=width,alpha=alpha))


def direction_arrow(ax,xyz,color):
    # Use an actual directed edge near the branch midpoint, never a guessed axis.
    i=max(0,min(len(xyz)//2,len(xyz)-2));vector=xyz[i+1]-xyz[i]
    norm=np.linalg.norm(vector)
    if norm>0:
        # A short glyph cannot be mistaken for an additional straight branch.
        # This is a display length only; direction is the original directed edge.
        length=min(max(np.linalg.norm(np.diff(xyz,axis=0),axis=1).sum()*.06,norm),3.)
        vector=vector/norm*length
        ax.quiver(*xyz[i],*vector,color=color,linewidth=2,arrow_length_ratio=.35)


def label(ax,xyz,text,*,color='#111827',size=12):
    ax.text(*xyz,text,fontsize=size,color=color,zorder=20,
        bbox=dict(boxstyle='round,pad=.20',fc='white',ec=color,alpha=.92,lw=.6))


def projected_labels(ax,entries,*,size=12):
    """Readable static callouts: sorted outside the geometry, with leader lines."""
    if not entries:return
    projection=ax.get_proj()
    projected=[proj3d.proj_transform(*xyz,projection)[:2] for xyz,text in entries]
    order=sorted(range(len(entries)),key=lambda i:projected[i][0])
    middle=max(1,len(order)//2)
    for column,indices in enumerate([order[:middle],order[middle:]]):
        indices=sorted(indices,key=lambda i:projected[i][1],reverse=True)
        ys=np.linspace(.88,.15,len(indices)) if len(indices)>1 else [.5]
        for index,y in zip(indices,ys):
            ax.annotate(entries[index][1],xy=projected[index],xycoords='data',
                xytext=(.01 if column==0 else .99,float(y)),textcoords='axes fraction',
                ha='left' if column==0 else 'right',va='center',fontsize=size,color='#111827',
                bbox=dict(boxstyle='round,pad=.2',fc='white',ec='#64748b',alpha=.96,lw=.7),
                arrowprops=dict(arrowstyle='-',color='#64748b',lw=.8),zorder=40)


def legend(fig,component=False):
    handles=[Line2D([0],[0],color=color,lw=4,label=name) if name!='Current boundary' else
        Line2D([0],[0],marker='o',color='white',markerfacecolor=color,markeredgecolor='#711020',markersize=10,label=name)
        for name,color in COLORS.items()]
    fig.legend(handles=handles,loc='lower center',ncol=5,fontsize=13,frameon=False,bbox_to_anchor=(.5,.012))


def new_figure(title):
    fig=plt.figure(figsize=SIZE,dpi=DPI,facecolor='white')
    fig.suptitle(title,fontsize=19,fontweight='bold',x=.025,ha='left',y=.976)
    fig.text(.025,.944,'BG001 RMCA | Original SWC world coordinates (mm) | Candidate evidence only | UNREVIEWED',fontsize=12,color='#475569')
    return fig


def finish(fig,path):
    fig.savefig(path,dpi=DPI,facecolor='white',pil_kwargs={'compress_level':9})
    plt.close(fig)


def boundary_figure(cache,event,detail,local,path):
    fig=new_figure(f'{event["boundary_id"]}  |  {event["boundary_type"]}')
    grid=fig.add_gridspec(2,2,left=.035,right=.965,bottom=.095,top=.90,wspace=.14,hspace=.22)
    a=fig.add_subplot(grid[0,0],projection='3d');b=fig.add_subplot(grid[0,1],projection='3d')
    c=fig.add_subplot(grid[1,0],projection='3d');d=fig.add_subplot(grid[1,1]);d.axis('off')
    a.set_title('A  GLOBAL CONTEXT',loc='left',fontsize=16,fontweight='bold',pad=4)
    b.set_title('B  LOCAL 3D VIEW',loc='left',fontsize=16,fontweight='bold',pad=4)
    c.set_title('C  SECOND VIEW (+90 degrees)',loc='left',fontsize=16,fontweight='bold',pad=4)
    d.set_title('D  CACHED NUMERICAL EVIDENCE',loc='left',fontsize=16,fontweight='bold',pad=8)
    all_keys=list(cache.branches);component=event['roi_component']
    selected=[k for k,v in cache.branch_component.items() if v==component]
    lines(a,cache,all_keys,colors={k:COLORS['Context'] for k in all_keys},alpha=.25,width=1.3)
    lines(a,cache,selected,width=3.2)
    xyz=cache.coordinates[event['boundary_node_id']][:3]
    a.scatter(*xyz,s=120,c=COLORS['Current boundary'],edgecolors='black',depthshade=False,zorder=30)
    label(a,xyz,f'part{component:02d} / node {event["boundary_node_id"]}',size=11)
    bounds(a,np.array(list(cache.coordinates.values()))[:,:3]);a.view_init(elev=24,azim=-65)
    keys=[branch['branch_id'] for branch in local['branches']]
    coords=np.concatenate([points(cache,k) for k in keys])
    for ax,azim in [(b,-65),(c,25)]:
        lines(ax,cache,keys,width=3.4)
        for key in keys:direction_arrow(ax,points(cache,key),COLORS[cache.branches[key]['label']])
        ax.scatter(*xyz,s=130,c=COLORS['Current boundary'],edgecolors='black',depthshade=False,zorder=30)
        callouts=[]
        for side,word in [('upstream','PROXIMAL'),('downstream','DISTAL')]:
            key=event[side+'_branch_id']
            if key in cache.branches:
                branch_points=points(cache,key)
                anchor=branch_points[max(0,len(branch_points)//2-1)]
                callouts.append((anchor,f'{word}\n{key} [{cache.branches[key]["label"]}]'))
        if event['downstream_branch_id'] is None:callouts.append((xyz,'DISTAL\nSOURCE TERMINAL'))
        bounds(ax,coords);ax.view_init(elev=25,azim=azim)
        projected_labels(ax,callouts,size=12)
    fmt=lambda value: 'n/a' if value is None else f'{value:.3f}' if isinstance(value,(int,float)) else str(value)
    lines_text=[f'Component part{component:02d}   |   boundary node {event["boundary_node_id"]}',
        f'UPSTREAM   {event["upstream_branch_id"] or "none"} [{event["upstream_label"] or "none"}]',
        f'DOWNSTREAM {event["downstream_branch_id"] or "terminal"} [{event["downstream_label"] or "no branch"}]',
        '                              UPSTREAM       DOWNSTREAM',
        f'p(M1) / p(MeVO)      {fmt(event["upstream_p_M1"])} / {fmt(event["upstream_p_MeVO"])}      {fmt(event["downstream_p_M1"])} / {fmt(event["downstream_p_MeVO"])}',
        f'Known length fraction       {fmt(event["upstream_known_fraction"]):<12}   {fmt(event["downstream_known_fraction"])}',
        f'Adjacent-point agreement    {fmt(detail["donor_vote_breakdown"]["upstream"].get("agreement")):<12}   {fmt(detail["donor_vote_breakdown"]["downstream"].get("agreement"))}',
        f'Normalized support distance {fmt(event["upstream_median_support_distance"]):<12}   {fmt(event["downstream_median_support_distance"])}',
        f'Median radius (mm)          {fmt(event["upstream_radius_median"]):<12}   {fmt(event["downstream_radius_median"])}',
        f'Branch length (mm)          {fmt(event["upstream_length_mm"]):<12}   {fmt(event["downstream_length_mm"])}',
        f'Root arc distance: {event["distance_from_RMCA_root_mm"]:.2f} mm   |   angle: {fmt(detail["boundary_angle_degrees"])} deg',
        f'VascularMD: {event["vascularmd_component_status"]}',
        'Donor       UP vote / dist(mm)      DOWN vote / dist(mm)']
    up=detail['donor_vote_breakdown']['upstream']['donors'];down=detail['donor_vote_breakdown']['downstream']['donors']
    donor_ids=sorted({r['donor_case'] for r in up+down})
    for donor in donor_ids:
        values=[]
        for data in [up,down]:
            entry=next((v for v in data if v['donor_case']==donor),None)
            values.append(f'{entry["vote"]} / {entry["nearest_distance"]:.2f}' if entry and entry['nearest_distance'] is not None else 'n/a')
        lines_text.append(f'MRA{donor}      {values[0]:<20}    {values[1]}')
    if not donor_ids:lines_text.append('DONOR_DETAIL_UNAVAILABLE / no sample on absent side')
    lines_text+=['Directions show native parent -> child; not measured flow.',
        'Model entry failure does NOT imply semantic boundary failure.' if component==3 else 'Candidate status is separate from model/manufacturing status.']
    d.text(0,.99,'\n'.join(lines_text),va='top',fontsize=11.2,fontfamily='DejaVu Sans Mono',linespacing=1.35)
    legend(fig);finish(fig,path)


def overview(cache,events,path,kind,*,azim=-65,component=None):
    fig=new_figure(f'BG001 RMCA | {kind.upper()}'+(f' | part{component:02d}' if component else ''))
    ax=fig.add_axes([.035,.12,.70,.79],projection='3d')
    keys=list(cache.branches)
    if component:
        selected={k for k,v in cache.branch_component.items() if v==component}
        context={cache.branches[k]['parent_branch'] for k in selected if cache.branches[k]['parent_branch'] in cache.branches}
        neighbors={v for k in selected for v in cache.branch_graph.successors(k) if cache.branches[v]['label']=='UNKNOWN'}
        keys=sorted(selected|context|neighbors)
    color_map={}
    if kind=='components':color_map={k:COMPONENT_COLORS[(cache.branch_component[k]-1)%5] if k in cache.branch_component else COLORS['Context'] for k in keys}
    lines(ax,cache,keys,colors=color_map,width=3.3)
    bounds(ax,np.concatenate([points(cache,k) for k in keys]));ax.view_init(elev=26,azim=azim)
    sidebar=[]
    if kind=='components':
        for c in sorted(cache.components):
            points_all=np.concatenate([points(cache,k) for k,v in cache.branch_component.items() if v==c])
            label(ax,points_all.mean(axis=0),f'part{c:02d}',color=COMPONENT_COLORS[c-1],size=17)
            sidebar.append(f'part{c:02d}: {cache.components[c]["strict"]["original_branch_count"]} branches')
    if kind=='boundaries' or component:
        relevant=[e for e in events if component is None or e['roi_component']==component]
        groups={}
        for e in relevant:groups.setdefault(e['boundary_node_id'],[]).append(e)
        callouts=[]
        for node,entries in groups.items():
            xyz=cache.coordinates[node][:3]
            ax.scatter(*xyz,s=90,c=COLORS['Current boundary'],edgecolors='black',depthshade=False)
            # Events at a shared physical junction retain their individual IDs.
            # A compact location key connects markers to the fully listed IDs.
            marker=short(entries[0]['boundary_id'])
            callouts.append((xyz,marker+(' +' if len(entries)>1 else '')))
            if len(entries)>1:
                sidebar.append(marker+' +\n'+ '\n'.join('  '+short(e['boundary_id']) for e in entries[1:]))
        projected_labels(ax,callouts,size=11)
        sidebar.insert(0,'+ = multiple events at this node\nAdditional event IDs at shared nodes:')
    if component:
        sidebar=[f'part{component:02d}: '+cache.components[component]['VascularMD_status'],
                 'Semantic candidate: UNREVIEWED',*sidebar]
        for key in keys:direction_arrow(ax,points(cache,key),COLORS[cache.branches[key]['label']])
    if kind=='labels':
        count={label:sum(b['label']==label for b in cache.branches.values()) for label in ['M1','MeVO','UNKNOWN']}
        sidebar=[f'{label}: {n} branches' for label,n in count.items()]+[
            'MeVO = M2 OR M3','TopBrain-informed candidates','No manual ground truth',
            'TYPE 4 = RMCA','BG001-specific','DERIVED_SOURCE_SUPPORTED',
            'Coordinates unchanged','No classification recomputed']
    fig.text(.765,.87,'\n\n'.join(textwrap.fill(s,36,replace_whitespace=False) if '\n' not in s else s for s in sidebar),
        va='top',fontsize=11.5,linespacing=1.4)
    legend(fig);finish(fig,path)


def write_local_vtp(cache,event,local,path):
    """Compact VTP + JSON manifest: lines, boundary vertex and directed arrows."""
    import pyvista as pv
    vertices=[cache.coordinates[event['boundary_node_id']][:3]];node_ids=[event['boundary_node_id']]
    connectivity=[];roles=['boundary_marker'];branch_names=[''];label_codes=[-1]
    for b in local['branches']:
        xyz=points(cache,b['branch_id']);start=len(vertices)
        vertices.extend(xyz);node_ids.extend(b['node_ids']);connectivity.extend([len(xyz),*range(start,start+len(xyz))])
        roles.append('native_branch');branch_names.append(b['branch_id']);label_codes.append({'UNKNOWN':0,'M1':1,'MeVO':2}[b['label']])
        i=max(0,min(len(xyz)//2,len(xyz)-2));tail,head=xyz[i],xyz[i+1];direction=head-tail
        length=np.linalg.norm(direction)
        if length>0:
            direction/=length;axis=np.eye(3)[np.argmin(abs(direction))];normal=np.cross(direction,axis);normal/=np.linalg.norm(normal)
            tip=head;wing=max(length*.3,.15)
            for end in [tip-direction*wing+normal*wing*.5,tip-direction*wing-normal*wing*.5]:
                start=len(vertices);vertices.extend([tip,end]);node_ids.extend([-1,-1]);connectivity.extend([2,start,start+1])
                roles.append('direction_arrow');branch_names.append(b['branch_id']);label_codes.append(-1)
    mesh=pv.PolyData(np.asarray(vertices),verts=np.array([1,0]),lines=np.asarray(connectivity))
    mesh.point_data['original_swc_id']=node_ids;mesh.cell_data['role']=roles;mesh.cell_data['branch_id']=branch_names
    mesh.cell_data['label_code']=label_codes;mesh.field_data['boundary_id']=[event['boundary_id']]
    mesh.field_data['units']=['mm'];mesh.save(path)
    write_manifest=dict(file=Path(path).name,format='VTP with vertex and directed line cells',
        role_array='boundary_marker / native_branch / direction_arrow',boundary_id=event['boundary_id'],
        labels={'0':'UNKNOWN','1':'M1','2':'MeVO','-1':'review marker/arrow'},
        direction='Native parent-to-child; arrowheads are display glyphs only',colors=COLORS)
    from .topbrain_qc import write_json
    write_json(Path(path).with_suffix('.json'),write_manifest)
