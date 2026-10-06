"""Scientific figures from saved arrays only; no smoothing or hidden invalid rows."""
from common import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import csv

def table(name,header,rows):
    with (HERE/name).open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n');writer.writerow(header);writer.writerows(rows)

def main():
    plt.rcParams.update({'font.size':11,'axes.titlesize':12,'axes.spines.top':False,'axes.spines.right':False})
    b=read(HERE/'P2_REPRESENTATION_FIXTURES.json');a=np.load(HERE/'P2_REPRESENTATION_FIXTURES.npz')
    pipe=read(HERE/'P2_POISEUILLE_REPRESENTATION.json');geom=read(HERE/'P2_CURVED_GEOMETRY_AUDIT.json')
    fig=plt.figure(figsize=(13,10),layout='constrained');gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[0,0]);keys=['u','G','curl','strain','s'];xx=np.arange(len(keys))
    for offset,key,color in [(-.18,'P1','#d57225'),(.18,'P2','#2366aa')]:
        vals=[b['B1'][key][k]['scaled_max_error'] for k in keys]
        ax.bar(xx+offset,vals,width=.36,color=color,label=key)
    ax.axhline(1e-11,ls='--',color='k',label='Declared roundoff tolerance')
    ax.set(yscale='log',xticks=xx,xticklabels=keys,ylabel='Scaled maximum absolute error',title='A. Full quadratic field: native basis')
    ax.legend(fontsize=9)
    ax=fig.add_subplot(gs[0,1]);order=np.argsort(a['B2_points'][10:,2]);z=a['B2_points'][10:,2][order]
    ax.plot(z,a['B2_exact_normal'][10:][order],color='k',lw=2,label='Exact: z²')
    ax.plot(z,a['B2_P1_normal'][10:][order],color='#d57225',ls='--',label='P1')
    ax.scatter(z[::65],a['B2_P2_normal'][10:][order][::65],s=14,color='#2366aa',label='P2 native samples')
    ax.set(xlabel='Distance z from planar wall (dimensionless)',ylabel='Signed normal velocity (dimensionless)',title='B. One wall face: representational rank freedom')
    ax.text(.03,.7,'4096 interior points\nP1 rank(G) = 1\nP2 rank(G) = 3',transform=ax.transAxes,fontsize=10)
    ax.legend(loc='lower right',fontsize=9)
    ax=fig.add_subplot(gs[1,0],projection='3d');nodes=a['real_138584_nodes'];origin=nodes[:4].mean(0);local=(nodes-origin)*1e6
    faces=[local[f['local_vertices']] for f in b['real_138584']['wall_faces']]
    ax.add_collection3d(Poly3DCollection(faces,alpha=.18,facecolor='#68798b',edgecolor='none'))
    for edge in b['real_138584']['edges']:
        i,j=edge['local_edge'];line=local[[i,j]];free=not edge['fixed_velocity']
        ax.plot(*line.T,color='#168349' if free else '#7b7b7b',lw=3 if free else 1)
    fixed=a['real_138584_fixed_nodes'];ax.scatter(*local[fixed].T,color='#bd4743',s=23,label='Wall-constrained nodes')
    ax.scatter(*local[~fixed].T,color='#168349',s=60,label='One free edge node')
    for i in [0,1,2,3,8]:ax.text(*local[i],str(i),fontsize=10)
    ax.set(xlabel='x offset (µm)',ylabel='y offset (µm)',zlabel='z offset (µm)',title='C. Actual original tetra 138584 topology')
    ax.set_box_aspect(np.ptp(local,axis=0),zoom=.8);ax.legend(loc='upper right',fontsize=8)
    ax=fig.add_subplot(gs[1,1]);ax.axis('off')
    rows=[[str(e['local_edge']),str(e['wall_triangle_count']),'Fixed' if e['fixed_velocity'] else 'FREE'] for e in b['real_138584']['edges']]
    t=ax.table(cellText=rows,colLabels=['Local edge','WALL triangles','Midpoint'],loc='upper center',cellLoc='center',bbox=[0,.42,1,.5]);t.auto_set_font_size(False);t.set_fontsize(11)
    ax.text(0,.31,'4 constrained vertices + 5 constrained midpoints\n= 1 free velocity node (3 scalar DOFs).',fontsize=12)
    ax.text(0,.12,'P2 need not be identically zero here.\nOne free node still implies rank(G) ≤ 1.\nNo actual P2 CFD field was computed.',fontsize=11)
    ax.set_title('D. Topology conclusion, not a CFD fix')
    fig.suptitle('P1 vs native P2 representation — validation only',fontsize=16)
    fig.savefig(HERE/'01_P1_VS_P2_REPRESENTATION.png',dpi=300);plt.close(fig)
    table('01_REPRESENTATION_METRICS.csv',['fixture','space','quantity','scaled_max_error'],
          [[f,s,k,v['scaled_max_error']] for f in ['B1','B2'] for s in ['P1','P2'] for k,v in b[f][s].items()])
    table('01_REAL_138584_EDGES.csv',['midpoint_local_id','endpoint0_GlobalNodeID','endpoint1_GlobalNodeID','wall_triangle_count','fixed_velocity'],
          [[e['local_midpoint_node'],*e['GlobalNodeIDs_1based'],e['wall_triangle_count'],e['fixed_velocity']] for e in b['real_138584']['edges']])
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    groups=pipe['groups'];xx=np.arange(len(groups));labels=[f"{r['wall_distance_m']*1e6:.3f}" for r in groups]
    for ax,k in zip(axes.flat,['u','G','s','curl']):
        for space,color,marker in [('P1','#d57225','o'),('P2','#2366aa','s')]:
            vals=np.array([r[space][k]['relative_L2']*100 for r in groups])
            positive=vals>0;ax.plot(xx[positive],vals[positive],marker=marker,color=color,label=f'Exact-node {space}')
            for i in np.flatnonzero(~positive):ax.annotate('exact zero',(i,1e-14),fontsize=8)
        ax.set(yscale='log',xticks=xx,xticklabels=labels,xlabel='Analytic wall distance (µm)',ylabel='Relative L2 error (%)',title=f'{k}: exact-node interpolation vs analytic field')
        ax.tick_params(axis='x',rotation=45);ax.grid(axis='y',alpha=.2);ax.legend(fontsize=9)
    fig.suptitle('Frozen V2: affine P2 representation (no Navier–Stokes solve)\n960 unchanged probes + 4096 random interior points; all valid',fontsize=15)
    fig.savefig(HERE/'02_P1_VS_P2_POISEUILLE_REPRESENTATION.png',dpi=300);plt.close(fig)
    table('02_POISEUILLE_REPRESENTATION.csv',['wall_distance_um','role','samples','quantity','P1_relative_L2','P2_relative_L2'],
          [[r['wall_distance_m']*1e6,r['role'],r['samples'],k,r['P1'][k]['relative_L2'],r['P2'][k]['relative_L2']] for r in groups for k in ['u','G','s','curl']])
    fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained');colors=['#2366aa','#a23c69']
    for f,color in zip(geom['fixtures'],colors):
        q=f['variants']['curved'];name=f['label'].replace('_',' ')
        axes[0].plot(np.arange(4),q['continuity_used_over_velocity_Jacobian'],marker='o',color=color,label=name)
    axes[0].axhline(1,color='k',ls='--',label='Common curved geometry')
    axes[0].set(xlabel='Native continuity quadrature index',ylabel='Used Jacobian / velocity Jacobian',xticks=range(4),title='A. Only the first point uses affine Jacobian');axes[0].legend(fontsize=8)
    for i,f in enumerate(geom['fixtures']):
        for offset,var,color in [(-.17,'affine','#7d8793'),(.17,'curved',colors[i])]:
            axes[1].bar(i+offset,100*f['variants'][var]['continuity_geometry_weight_relative_L2_mismatch'],.34,color=color)
    axes[1].set(yscale='log',xticks=[0,1],xticklabels=['Central wall','Inlet rim'],ylabel='Native Rp vs common-geometry Rp (%)',title='B. Direct native construct_fluid residual')
    axes[1].text(.04,.9,'Gray: affine control\nColor: projected midpoints',transform=axes[1].transAxes,fontsize=9)
    axes[2].axis('off');axes[2].set_title('C. Geometry gate blocks GPU CFD')
    axes[2].text(0,.92,'PASS\n• Native P2 exact representation\n• Positive curved Jacobians (certified bound)\n• Wall midpoint and cap-plane constraints\n\nBLOCKED\n• Mixed pressure / velocity mappings differ\n• Native continuity uses inconsistent weights\n\nP2 GPU pipe / BraVa: NOT_RUN\nP2 trajectory reliability: NOT_EVALUABLE',va='top',fontsize=11)
    fig.suptitle('Native curved mixed geometry audit — two actual V2 cells, no CFD time steps',fontsize=15)
    fig.savefig(HERE/'03_CURVED_GEOMETRY_NATIVE_ASSEMBLY.png',dpi=300);plt.close(fig)
    table('03_CURVED_GEOMETRY_AUDIT.csv',['fixture','geometry','Rp_relative_L2_mismatch','pressure_gradient_relative_L2_mismatch','max_mapping_distance_m'],
          [[f['label'],v,q['continuity_geometry_weight_relative_L2_mismatch'],q['pressure_gradient_common_mapping_relative_L2'],q['max_pressure_coordinate_mismatch_m']] for f in geom['fixtures'] for v,q in f['variants'].items()])
    put('FIGURE_MANIFEST.json',dict(dpi=300,smoothing=False,invalid_states_dropped=0,
        logarithmic_axis_reason='Representation and residual errors span roundoff to percent levels; linear quadrature ratio and geometry axes preserve interpretation.',
        figures=[dict(file=p.name,sha256=sha(p)) for p in sorted(HERE.glob('0*.png'))],
        data=['P2_REPRESENTATION_FIXTURES.npz','P2_POISEUILLE_SAMPLES.npz','P2_CURVED_GEOMETRY_SAMPLES.npz'],
        D1_warning='0.1 µm is diagnostic only. Exact-node affine wall-midpoint values are not a no-slip CFD solution; curved-geometry semantic failure is a separate gate.'))
    print('3 figures saved, 300 dpi, CSV/NPZ source retained')

if __name__=='__main__':main()
