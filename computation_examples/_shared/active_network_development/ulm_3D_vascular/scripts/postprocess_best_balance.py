"""Read-only design/CFD comparison, using the archived validated WSS core."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from network_1d0d.audit import sha256,write_json,write_csv


def area_quantile(values,areas,q):
    order=np.argsort(values,kind='stable')
    cumulative=np.cumsum(areas[order])/areas.sum()
    return float(values[order[min(np.searchsorted(cumulative,q),len(order)-1)]])


def summarize(label,region,values,areas):
    return dict(case=label,region=region,face_count=len(values),area_um2=float(areas.sum()*1e12),
        area_mean_Pa=float(np.average(values,weights=areas)),
        P5_area_Pa=area_quantile(values,areas,.05),P50_area_Pa=area_quantile(values,areas,.5),
        P95_area_Pa=area_quantile(values,areas,.95),min_Pa=float(values.min()),max_Pa=float(values.max()),
        area_below_1Pa_um2=float(areas[values<1].sum()*1e12),
        area_above_30Pa_um2=float(areas[values>30].sum()*1e12))


def run(case,source,flow_root,design,report):
    case,source,flow_root,design,report=map(Path,(case,source,flow_root,design,report))
    acceptance=json.loads((report/'independent_validation.json').read_text())
    assert acceptance['status']=='PASS'
    core=report/'evidence/wss_core'
    provenance=json.loads((core/'provenance.json').read_text())
    for name,digest in provenance['source_hashes'].items():
        assert sha256(core/'flow_solver_support'/name)==digest
    sys.path[:0]=[str(core),str(flow_root/'src')]
    from flow_solver_support.wss_case import recover,material,coordinate_identity
    from sv_validation.postprocess import SolutionMeasurements
    import pyvista as pv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    policy=json.loads((case/'policy.json').read_text())
    meas=SolutionMeasurements(case/'SV_MESH/mesh_arrays.npz',policy['Q_target_m3_s'],policy['Umean_m_s'])
    paths={'old_H0_3D':source/'frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu',
           'best_feasible_3D':case/'frozen_flow/steady_flow_mean_2p0_mmps_A_best_feasible_balance.vtu'}
    data=report/'data';data.mkdir(exist_ok=True)
    figures=report/'figures';figures.mkdir(exist_ok=True)
    stats=[];meshes={};metadata={}
    for label,path in paths.items():
        actual_case=source if label=='old_H0_3D' else case
        mat=material(actual_case)
        assert mat['mu_Pa_s']==.00345312
        grid=pv.read(path);identity=coordinate_identity(grid.points,meas.points)
        np.testing.assert_array_equal(np.sort(grid.cells.reshape(-1,5)[:,1:],axis=1),np.sort(meas.tetra,axis=1))
        u,p=meas.read(path)
        surface,detail=recover(meas.points,meas.tetra,meas.boundary,meas.tags,u,p,mat['mu_Pa_s'])
        output=data/(label+'_wall_wss.vtp');surface.save(output)
        values,areas,centers=detail['magnitude'],detail['area'],detail['centers']
        np.savez_compressed(data/(label+'_wss_raw_si.npz'),WSS_Pa=values,area_m2=areas,centers_m=centers,
            wall_boundary_facet_zero_based=detail['wall_ids'],parent_tetra_zero_based=detail['owners'])
        masks={'whole_wall_including_extensions':np.ones(len(values),bool)}
        for region,center in [('J1_ball_5um',(92,49,111)),('J2_ball_5um',(130.04,82.04,87.18))]:
            masks[region]=np.linalg.norm(centers-np.asarray(center)*1e-6,axis=1)<=5e-6
        for region,mask in masks.items():
            assert mask.any(),region
            stats.append(summarize(label,region,values[mask],areas[mask]))
        metadata[label]=dict(source=str(path.resolve()),source_sha256=sha256(path),material=mat,
            coordinate_identity=identity,wall_surface_sha256=sha256(output),wall_face_count=len(values))
        meshes[label]=surface
    write_csv(report/'wss_comparison.csv',stats)
    # P1 raw wall-facet magnitudes only: no nodal interpolation or smoothing in the figure.
    panel=[]
    for label,surface in meshes.items():
        plotter=pv.Plotter(off_screen=True,window_size=(1400,1300))
        plotter.set_background('white')
        plotter.add_mesh(surface,scalars='WSS_raw_Pa',preference='cell',clim=(0,55),cmap='turbo',
            smooth_shading=False,lighting=False,show_scalar_bar=False)
        plotter.view_vector((1,-1,1),viewup=(0,0,1))
        plotter.enable_parallel_projection();plotter.reset_camera()
        panel.append(plotter.screenshot(return_img=True));plotter.close()
    fig,axs=plt.subplots(1,2,figsize=(12,6.4),layout='constrained')
    for ax,img,label in zip(axs,panel,('Old H0: same mesh','Best feasible design: one CFD')):
        ax.imshow(img);ax.set_title(label);ax.axis('off')
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,55),cmap='turbo'),ax=axs,
        label='Raw wall-facet WSS [Pa]; linear 0–55',shrink=.7)
    fig.suptitle('Same camera, geometry and scale; no smoothing / no CFD retuning')
    for suffix in ('png','pdf'):fig.savefig(figures/('wss_raw_comparison.'+suffix),dpi=240)
    plt.close(fig)
    comp=acceptance['comparison'];summary=json.loads((design/'design_summary.json').read_text())
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    grid=np.genfromtxt(design/'coarse_grid.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    x=np.unique(grid['s_O1']);y=np.unique(grid['s_O2'])
    z=np.empty((len(y),len(x)))
    for row in grid:z[np.searchsorted(y,row['s_O2']),np.searchsorted(x,row['s_O1'])]=row['J_balance']
    heat=axs[0].pcolormesh(x,y,z,shading='nearest',cmap='viridis')
    opt=summary['optimized']['parameters']
    axs[0].plot(opt['s_O1'],opt['s_O2'],'*',ms=14,color='red',clip_on=False,label='Frozen design')
    axs[0].plot(1,1,'wo',mec='black',label='Baseline')
    axs[0].set(xscale='log',yscale='log',xlabel='s_O1',ylabel='s_O2',title='Pure 0D: 25 × 25 log grid')
    axs[0].legend(fontsize=8);fig.colorbar(heat,ax=axs[0],label='J balance')
    positions=np.arange(3)
    for offset,(label,key) in zip((-.26,0,.26),[('Old H0 3D','old_H0_3D_fractions'),('Frozen 0D design','design_0D_fractions'),('New 3D forward','new_3D_fractions')]):
        vals=100*np.asarray(comp[key]);bars=axs[1].bar(positions+offset,vals,width=.25,label=label)
        axs[1].bar_label(bars,fmt='%.2f',fontsize=7)
    axs[1].axhline(100/3,color='black',ls='--',lw=.8,label='Equal split: not a gate')
    axs[1].set(xticks=positions,xticklabels=['O1','O2','O3'],ylabel='Signed Qout / Qin [%]',ylim=(0,65),title='Actual forward prediction')
    axs[1].legend(fontsize=7)
    for suffix in ('png','pdf'):fig.savefig(figures/('balance_design_and_forward.'+suffix),dpi=240)
    plt.close(fig)
    write_json(report/'wss_postprocess.json',dict(status='PASS',method='Unchanged archived production P1 tensor traction; only wall tag 1',
        formula='||(I-nn^T) mu (grad u + grad u^T) n||',case_sources=metadata,core=provenance,
        statistics='Area-weighted empirical inverse-CDF quantiles, raw face magnitudes; 5 um balls selected by face centroids',
        global_wall_includes_artificial_extensions=True,render='Same orthographic camera, raw cell values, 0-55 Pa, no smoothing',
        region_centers_um={'J1':[92,49,111],'J2':[130.04,82.04,87.18]},
        no_new_solver_calls=True,mesh_independence_claim=False))
    print(json.dumps({'status':'PASS','statistics':stats},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('case','source-case','flow-root','design','report'):parser.add_argument('--'+key,type=Path,required=True)
    a=parser.parse_args();run(a.case,a.source_case,a.flow_root,a.design,a.report)
