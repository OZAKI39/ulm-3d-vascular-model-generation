"""One static renderer check on saved V2 coordinates; no new dynamics or movie."""
from render_videos import *

def main():
    case=V2_O/'cases/biphasic';ref=load_npz(V2_O/'provenance/reference.npz');surf=CurrentSurface(ref,'uniform')
    old=state(dict(path=case/'keys/final.npz'));lookup,reg,_,_=registry(case);top=topology(ref,old,lookup)
    ids=np.array(reg['3769']['atom_ids'])-1;label=int(top['labels'][ids[0]]);selected=top['labels']==label
    z=surf.evaluate(old['x'],old['velocity'],old['d'],float(old['time']));new=old.copy()
    for name,key in [('surface_active','active'),('surface_centroid','centroid'),('surface_force','face_force'),('surface_area','area'),('surface_weight','weight'),('surface_normal','normal'),('surface_velocity','face_velocity'),('surface_outer_velocity','outer_velocity'),('surface_unresolved','unresolved')]:new[name]=z[key]
    xyz=old['x'][ids]*1e6;center=(xyz.min(0)+xyz.max(0))/2;width=max(np.ptp(xyz,axis=0))*1.2
    bounds=np.ravel(np.column_stack([center-width/2,center+width/2]));pl=pv.Plotter(shape=(1,2),off_screen=True,window_size=(WIDTH,HEIGHT),border=False)
    details=[];widgets=[]
    for j,s in enumerate([old,new]):
        axes=base_view(pl,j,bounds)
        widgets.append(axes)
        _,info=draw_scene(pl,ref,s,top,surf,z,j==1,f'diag{j}_',selected,True);details.append(info)
        set_camera(pl,bounds);update_axis(axes,bounds)
    raw=pl.screenshot(return_img=True);raw=overlay(raw,1,dict(old=old,new=new),dict(old=top,new=top),details,dict(old='OFFLINE',new='OFFLINE'),'3769',annotations=axis_annotations(pl,widgets));pl.close()
    im=Image.fromarray(raw);draw=ImageDraw.Draw(im);font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',30)
    draw.rectangle((0,0,1920,198),fill='#fcfcfc');draw.text((28,12),'OFFLINE: identical saved V2 H fragment ID 3769 (148 nodes); no new dynamics',font=font,fill='#243c4b')
    draw.text((28,68),'V2 recorded surface force',font=font,fill='#243c4b');draw.text((990,68),'V3 quadrature replay on the same V2 state',font=font,fill='#243c4b')
    draw.text((28,125),'Renderer check: actual geometry, surviving original bonds, retained faces and partial dual patches',font=font,fill='#243c4b')
    draw.rectangle((0,895,1920,997),fill='#fcfcfc');small=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',22)
    for j,q in enumerate(details):
        x=j*960+30
        draw.text((x,904),f'Selected: {q["nodes"]} nodes, {q["internal_live_bonds"]} live bonds; L1/L2 = {q["L1_facets"]}/{q["L2_patches"]}',font=small,fill='#405563')
        draw.text((x,938),f'Area {q["surface_area_m2"]*1e6:.5g} mm²; |F| {np.linalg.norm(q["force_N"])*1e6:.5g} uN',font=small,fill='#405563')
        draw.text((x,972),'Local view of saved coordinates; no artificial surface shell.',font=small,fill='#405563')
    path=O/'checks/P0_OFFLINE_RENDER_PREVIEW.png';im.save(path)
    dump(O/'checks/P0_OFFLINE_RENDER_PREVIEW.json',dict(status='RENDERED_REQUIRES_VISUAL_REVIEW',file=str(path.relative_to(R)),sha256=sha(path),physical_runs_started=0,actual_coordinates_source=str((case/'keys/final.npz').relative_to(R))))
    print('OFFLINE_RENDER_PREVIEW_WRITTEN',path,flush=True)

if __name__=='__main__':main()
