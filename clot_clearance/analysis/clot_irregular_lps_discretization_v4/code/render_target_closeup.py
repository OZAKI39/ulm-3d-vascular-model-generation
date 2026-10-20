"""Supplementary final target views; original coordinates and complete live graph."""
from render_videos import *

def main():
 refs=[load_npz(V3_O/'provenance/reference.npz'),load_npz(O/'provenance/reference.npz')]
 end=json.loads((O/'NATIVE_V4_RESULT.json').read_text())['accepted_time_s']
 sources=[choose(state_paths(c),end)[1] for c in [OLD,CASE]]
 states=[state(p) for p in sources];records=[];selections=[]
 for ref,s,p in zip(refs,states,sources):
  lab,attached=topology(ref,s['d']);counts=np.bincount(lab)
  labels=np.flatnonzero(~attached&(counts>=20)&(counts<=200))
  assert len(labels)==1,'This supplementary comparison expects the observed one final target per case.'
  ids=np.flatnonzero(lab==labels[0]);b=ref['bonds'];edges=b[(s['d']<1)&(lab[b[:,0]]==labels[0])&(lab[b[:,1]]==labels[0])]
  selections.append((ids,edges))
  records.append(dict(source=str(p.relative_to(R)),source_sha256=sha(p),time_s=float(s['time']),original_atom_ids=(ids+1).tolist(),nodes=len(ids),internal_live_bonds=len(edges),mass_kg=float(RHO*ref['V'][ids].sum())))
 span=max(np.ptp(s['x'][ids]*1e3,axis=0).max() for s,(ids,_) in zip(states,selections))+.16
 images=[]
 for name,ref,s,(ids,edges),record in zip(['V3 regular','V4 irregular'],refs,states,selections,records):
  xyz=s['x'][ids]*1e3;center=(xyz.min(0)+xyz.max(0))/2;bounds=np.c_[center-span/2,center+span/2].ravel();pl=setup(bounds,1)
  points=pv.PolyData(xyz);points['radius']=np.cbrt(3*ref['V'][ids]/(4*np.pi))*1e3
  glyph=points.glyph(geom=pv.Sphere(radius=1.,theta_resolution=16,phi_resolution=16),scale='radius',orient=False)
  pl.add_mesh(glyph,color='#a51c3c',opacity=.6,smooth_shading=True,ambient=.25,show_scalar_bar=False)
  pl.add_mesh(line_mesh(s['x']*1e3,edges),color='#500e24',line_width=1.5,lighting=False)
  pl.render();im=annotate(pl,pl.screenshot(return_img=True),bounds,[dict(time_s=float(s['time']),N20_200=1)],[name+' — final target close-up'],1)
  draw=ImageDraw.Draw(im);draw.rectangle((0,925,WIDTH,HEIGHT),fill='white');font=lambda size:ImageFont.truetype(FONT,size)
  draw.text((WIDTH/2,940),f'{len(ids)} original nodes | {len(edges)} live internal bonds | {record["mass_kg"]*1e9:.3f} µg',font=font(26),fill='#642136',anchor='mt')
  draw.text((WIDTH/2,992),'Selected component only; all internal surviving original bonds shown.',font=font(24),fill='#303b44',anchor='mt')
  draw.text((WIDTH/2,1036),'Actual coordinates ×1; equal zoom; world Z vertical; no continuous shell.',font=font(22),fill='#59636c',anchor='mt')
  images.append(im);pl.close();record['camera_bounds_mm']=bounds
 combined=Image.new('RGB',(2*WIDTH,HEIGHT),'white')
 for j,im in enumerate(images):combined.paste(im,(j*WIDTH,0))
 path=F/'figures/FIG_SUPP_FINAL_TARGET_CLOSEUP.png';combined.save(path)
 dump(O/'TARGET_CLOSEUP_PROVENANCE.json',dict(status='PASS',source_sha256=sha(Path(__file__)),image_sha256=sha(path),common_view_span_mm=float(span),displacement_scale=1.,display_only_subset=True,all_internal_surviving_bonds=True,fragment_membership='SURVIVING_ORIGINAL_MATERIAL_GRAPH',panels=records))
 print('TARGET_CLOSEUP_COMPLETE',flush=True)

if __name__=='__main__':main()
