#!/usr/bin/env python3
"""WSL-only review figures from final volume artifacts; display units are um."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fem3d.audit import sha256, timestamp, write_json

BG="#f7f9fc"
INK="#172e40"
ORANGE="#e67e36"
TEAL="#09828b"


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--coarse-only",action="store_true")
    args=parser.parse_args()
    profiles=["coarse"] if args.coarse_only else ["coarse","medium"]
    out=ROOT/"reports/stage01"/("coarse_review" if args.coarse_only else "")
    out.mkdir(parents=True,exist_ok=True)
    selected=profiles[-1]
    data=np.load(ROOT/f"outputs/stage01/{selected}/mesh/volume_mesh.npz")
    qc=json.loads((ROOT/f"outputs/stage01/{selected}/qc/geometry_qc.json").read_text())
    points=data["points_m"]*1e6  # display copy only; meshing and QC stay in SI
    triangles=data["boundary_triangles"]
    surface=pv.PolyData(points,np.c_[np.full(len(triangles),3),triangles].ravel())
    surface.cell_data["facet_tags"]=data["facet_tags"]
    tetra=data["tetra"]
    grid=pv.UnstructuredGrid(np.c_[np.full(len(tetra),4),tetra].ravel(),np.full(len(tetra),pv.CellType.TETRA,dtype=np.uint8),points)
    grid.cell_data["quality"]=data["min_sicn"]
    wall=surface.extract_cells(data["facet_tags"]==1).extract_surface()
    ports=list(qc["ports"].items())

    def plotter(shape,size):
        plot=pv.Plotter(off_screen=True,shape=shape,window_size=size,border=False)
        for renderer in plot.renderers:
            renderer.set_background(BG)
        return plot

    def save(plot,name):
        plot.screenshot(out/name)
        plot.close()

    plot=plotter((1,2),(2200,1200))
    for i in range(2):
        plot.subplot(0,i)
        plot.add_mesh(wall,color="#a4adb5",smooth_shading=False)
        labels=[]
        positions=[]
        for name,port in ports:
            color=ORANGE if name=="inlet" else TEAL
            plot.add_mesh(surface.extract_cells(data["facet_tags"]==port["entity_id"]),color=color,lighting=False)
            center=np.array(port["centroid_m"])*1e6
            normal=np.array(port["outward_normal"])
            plot.add_mesh(pv.Arrow(start=center,direction=normal,scale=9),color=color)
            positions.append(center+normal*13)
            labels.append(name)
        plot.add_point_labels(np.array(positions),labels,font_size=21,text_color=INK,
            shape_color="white",shape_opacity=.93,show_points=False,always_visible=True)
        plot.add_text(f"Final {selected} boundary | "+("isometric" if i==0 else "XY view"),font_size=18,color=INK)
        plot.add_text("Gray: WALL   Orange: INLET   Teal: OUTLETS\nArrows: normals from adjacent tetrahedra | display: um",position="lower_left",font_size=12,color=INK)
        plot.enable_parallel_projection()
        plot.view_isometric() if i==0 else plot.view_xy()
        plot.reset_camera()
        plot.camera.zoom(.85)
        plot.add_axes()
    save(plot,"boundary_tags.png")

    plot=plotter((2,2),(1800,1500))
    for i,(name,port) in enumerate(ports):
        plot.subplot(i//2,i%2)
        center=np.array(port["centroid_m"])*1e6
        normal=np.array(port["outward_normal"])
        color=ORANGE if name=="inlet" else TEAL
        cap=surface.extract_cells(data["facet_tags"]==port["entity_id"]).extract_surface()
        nearby=wall.extract_cells(np.linalg.norm(wall.cell_centers().points-center,axis=1)<5)
        plot.add_mesh(nearby,color="#b5bcc3",opacity=.65)
        plot.add_mesh(cap,color=color,show_edges=True,edge_color="#35505b",line_width=1)
        plot.add_mesh(pv.Arrow(start=center,direction=normal,scale=3.5),color=color)
        plot.add_text(f"{name} | entity {port['entity_id']} | {port['facet_count']} frozen triangles",font_size=15,color=INK)
        plot.add_text(f"Outward arrow | area {port['area_m2']*1e12:.4f} um^2\nTriangle edges shown; cap geometry unchanged",position="lower_left",font_size=12,color=INK)
        side=np.cross(normal,[0,0,1])
        if np.linalg.norm(side)<.1:
            side=np.cross(normal,[0,1,0])
        side/=np.linalg.norm(side)
        plot.camera_position=[center+12*(normal+.65*side),center,[0,0,1]]
        plot.enable_parallel_projection()
        plot.reset_camera()
        plot.camera.zoom(.9)
    save(plot,"port_closeups.png")

    inlet=qc["ports"]["inlet"]
    normal=np.array(inlet["outward_normal"])
    center=np.array(inlet["centroid_m"])*1e6-4*normal
    centers=grid.cell_centers().points
    local=grid.extract_cells(np.linalg.norm(centers-center,axis=1)<5)
    clipped=local.clip(normal=normal,origin=center,invert=True)
    sliced=local.slice(normal=normal,origin=center)
    if not sliced.n_cells or not clipped.n_cells:
        raise RuntimeError("Empty volume cutaway")
    plot=plotter((1,2),(2000,1100))
    for i,piece in enumerate((clipped,sliced)):
        plot.subplot(0,i)
        plot.add_mesh(piece,scalars="quality",cmap="viridis",clim=(0,1),show_edges=True,edge_color="#627480",line_width=.55,
            scalar_bar_args={"title":"Tetra quality (higher is better)","vertical":False,"position_y":.04,"height":.06,"color":INK})
        plot.add_text(f"{selected}: "+("cut volume, inlet branch" if i==0 else "actual tetra intersections"),font_size=17,color=INK)
        plot.add_text("Plane 4 um inside inlet; local radius 5 um\nColored polygons/cells come from volume tetrahedra",position=(20,140),font_size=12,color=INK)
        side=np.cross(normal,[0,0,1]); side/=np.linalg.norm(side)
        direction=normal+.65*side if i==0 else normal
        plot.camera_position=[center+15*direction,center,[0,0,1]]
        plot.enable_parallel_projection()
        plot.reset_camera()
        plot.camera.zoom(.82)
    save(plot,"tetrahedral_cutaway.png")

    fig,axes=plt.subplots(1,len(profiles),figsize=(7*len(profiles),5),squeeze=False)
    for ax,profile in zip(axes[0],profiles):
        values=np.load(ROOT/f"outputs/stage01/{profile}/mesh/volume_mesh.npz")["min_sicn"]
        v=np.quantile(values,[0,.05,.5])
        ax.hist(values,bins=np.linspace(0,1,61),color=TEAL,alpha=.85)
        for number,label,color in zip(v,["min","P5","median"],["#ca4d32","#d99b24",INK]):
            ax.axvline(number,color=color,linewidth=1.5,label=f"{label} = {number:.4f}")
        ax.set(title=f"{profile}: tetra quality (higher is better)",xlabel="Gmsh minSICN: 0 = flat, 1 = regular tetra",ylabel="Number of tetrahedra",xlim=(0,1))
        ax.legend(frameon=False)
        ax.grid(axis="y",alpha=.2)
    fig.tight_layout()
    fig.savefig(out/"mesh_quality_distribution.png",dpi=160)
    plt.close(fig)

    worst=np.argsort(data["min_sicn"])[:20]
    lowest=centers[worst[0]]
    plot=plotter((1,2),(2100,1200))
    plot.subplot(0,0)
    plot.add_mesh(wall,color="#b5bec6",opacity=.3)
    plot.add_mesh(grid.extract_cells(worst),color="#df462c",show_edges=True,edge_color="#852d20")
    plot.add_points(centers[worst],color="#df462c",point_size=9,render_points_as_spheres=True)
    plot.add_point_labels(np.array([lowest]),["worst"],font_size=19,always_visible=True,text_color=INK)
    plot.add_text(f"{selected}: 20 lowest-quality cells",font_size=18,color=INK)
    plot.add_text("Red markers show locations (enlarged for visibility)\nGray context is the frozen wall",position="lower_left",font_size=12,color=INK)
    plot.view_isometric(); plot.enable_parallel_projection(); plot.reset_camera(); plot.camera.zoom(.88)
    plot.subplot(0,1)
    localwall=wall.extract_cells(np.linalg.norm(wall.cell_centers().points-lowest,axis=1)<2.5)
    plot.add_mesh(localwall,color="#b5bec6",opacity=.2)
    nearworst=worst[np.linalg.norm(centers[worst]-lowest,axis=1)<2.5]
    plot.add_mesh(grid.extract_cells(nearworst),color="#df462c",show_edges=True,edge_color="#842816",line_width=2)
    plot.add_text(f"Worst cell #{worst[0]} | quality {data['min_sicn'][worst[0]]:.5f}",font_size=17,color=INK)
    plot.add_text("Actual tetrahedra, no enlargement in this panel\nPositive volume does not imply good element shape",position="lower_left",font_size=12,color=INK)
    plot.view_isometric(); plot.enable_parallel_projection(); plot.reset_camera(); plot.camera.zoom(.85)
    save(plot,"worst_elements.png")
    write_json(out/"visualization_manifest.json",{"timestamp":timestamp(),"profiles":profiles,"displayed_volume":selected,
        "display_units":"um; only temporary visualization copies multiplied by 1e6", "flow_fields_created":False,
        "cutaway":{"origin_m":(center*1e-6).tolist(),"normal":normal.tolist(),"cut_cells":clipped.n_cells,"slice_polygons":sliced.n_cells},
        "mesh_sha256":sha256(ROOT/f"outputs/stage01/{selected}/mesh/volume_mesh.npz"),
        "images":{p.name:{"sha256":sha256(p),"bytes":p.stat().st_size} for p in out.glob("*.png")}})
    print(f"Rendered 5 review figures in {out}")


if __name__=="__main__":
    main()
