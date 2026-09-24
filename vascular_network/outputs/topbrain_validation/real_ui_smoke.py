from pathlib import Path
import sys
root=Path('/home/lzy/projects/ulm_3D_vascular');sys.path.insert(0,str(root))
from vascular_processing.topbrain_dataset import load_case
from vascular_processing.topbrain_pipeline import process_case,TopBrainOptions
from vascular_processing.topbrain_display_adapter import TopBrainViewer
from vascular_processing.topbrain_qc import write_json
results=[]
for number in ['001','002']:
 case=load_case(root/'data/TopBrain',number,'mr',version='2025')
 run,rois,manifest=process_case(case,root/'outputs/topbrain_mevo',options=TopBrainOptions(export_swc=True,allow_derived_aca=True),cache_root=root/'outputs/topbrain_mevo_cache')
 viewer=TopBrainViewer(case,rois,run,manifest,show=True)
 if number=='001':
  def next_case(step):
   if step==5:
    viewer.plotter.screenshot(run/'topbrain_case_switch.png')
    viewer.request_case(1)
  viewer.plotter.add_timer_event(max_steps=7,duration=100,callback=next_case)
 report=viewer.run_window(show=True,smoke_seconds=3)
 results.append({'run':str(run),'case':number,'gui':report})
 if number=='001':assert report['case_delta']==1
write_json(root/'outputs/topbrain_validation/real_ui_case_switch.json',results)
