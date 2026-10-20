"""Tighter camera bounds for three static figures; no change to physical data/video."""
from render_results import *
import shutil
def bounds_for(s):
 x=s['x']*1e3;lo=x.min(0);hi=x.max(0);pad=.08*(hi-lo)+.08
 return np.c_[lo-pad,hi+pad].ravel()
def main():
 archive=O/'provenance/before_static_detail_refinement'
 assert not archive.exists(),'Detail refinement already recorded'
 archive.mkdir(parents=True)
 for name in ['ARTIFACT_MANIFEST.json','DELIVERY_VALIDATION.json','LOAD_VISUALIZATION_AUDIT.json']:
  shutil.copy2(O/name,archive/name)
 ref=load_npz(O/'provenance/reference.npz');pvlog=json.loads((O/'DYNAMICS_VIDEO_PROVENANCE.json').read_text())
 tasks=[('FIG_02_INITIAL_CLOT.png',CASE/'keys/initial.npz','Y clot — V4 irregular LPS'),
        ('FIG_02_LARGE_FRAGMENT_CLOT.png',R/pvlog['large_fragment_picture_source'],'Y clot — saved large-fragment state')]
 records=[]
 for name,path,title in tasks:
  image_path=F/'figures'/name;shutil.copy2(image_path,archive/name);s=state(path);bounds=bounds_for(s);pl=setup(bounds,1);st=scene(pl,0,ref,s,True);pl.render()
  annotate(pl,pl.screenshot(return_img=True),bounds,[st],[title],1).save(image_path);pl.close()
  records.append(dict(figure=name,source=str(path.relative_to(R)),source_sha256=sha(path),time_s=float(s['time']),camera_bounds_mm=bounds,coordinates_unchanged=True))
 name='FIG_03_ACTUAL_SURFACE_LOAD.png';shutil.copy2(F/'figures'/name,archive/name)
 old=load_npz(O/'ACTUAL_LOAD_FIELD.npz');old_hash=sha(O/'ACTUAL_LOAD_FIELD.npz')
 path=R/pvlog['load_picture_source'];s=load_npz(path);bounds=bounds_for(s)
 force_picture(ref,s,bounds,F/'figures'/name)
 new=load_npz(O/'ACTUAL_LOAD_FIELD.npz');assert old.keys()==new.keys()
 assert all(np.array_equal(old[k],new[k]) for k in old),'Physical load field must not change'
 records.append(dict(figure=name,source=str(path.relative_to(R)),source_sha256=sha(path),time_s=float(s['time']),camera_bounds_mm=bounds,coordinates_unchanged=True,all_load_field_arrays_exact=True,field_container_before_sha256=old_hash,field_container_after_sha256=sha(O/'ACTUAL_LOAD_FIELD.npz')))
 dump(O/'STATIC_DETAIL_RENDERING.json',dict(status='PASS',records=records,scope='Static figures use their own state bounds for legibility. Dynamics video retains one fixed global camera box across all times; rotation remains unchanged.',source_sha256=sha(A/'render_detail_stills.py')))
 print('STATIC_DETAIL_REFINEMENT_PASS',flush=True)
if __name__=='__main__':main()
