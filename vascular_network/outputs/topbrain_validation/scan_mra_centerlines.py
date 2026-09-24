from pathlib import Path
import json,time,sys
root=Path('/home/lzy/projects/ulm_3D_vascular');sys.path.insert(0,str(root))
from vascular_processing.topbrain_dataset import discover,load_case
from vascular_processing.topbrain_centerline import mca_centerline,graph_summary
from vascular_processing.topbrain_qc import components,TopBrainError,write_json
results=[]
for paths in discover(root/'data/TopBrain','2025'):
 if paths.modality!='mr':continue
 case=load_case(root/'data/TopBrain',paths.case_id,'mr',version='2025')
 for side in ['R','L']:
  m2,m3=case.mask(side+'-M2'),case.mask(side+'-M3')
  mask=m2|m3
  _,qc=components(mask,case.spacing_mm)
  entry={'case':paths.case_id,'side':side,**qc,'m2_voxels':int(m2.sum()),'m3_voxels':int(m3.sum())}
  if not m2.any() or not m3.any():entry['status']='MISSING_SEGMENT'
  elif qc['component_count']!=1:entry['status']='MULTIPLE_COMPONENTS'
  else:
   try:
    graph,_=mca_centerline(case,mask,side)
    entry.update(status='PASS',**graph_summary(graph),root_detection=graph.graph['root_qc'])
   except TopBrainError as exc:entry.update(status=exc.code,details=exc.details)
  results.append(entry)
  print(paths.case_id,side,entry['status'],flush=True)
 write_json(root/'outputs/topbrain_validation/mra_centerline_scan.json',results)
 if any(row['status']=='PASS' for row in results):break
