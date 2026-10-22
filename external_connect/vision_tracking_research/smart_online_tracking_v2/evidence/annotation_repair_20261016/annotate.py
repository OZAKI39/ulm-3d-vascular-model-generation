"""Official micro-sam Annotator2D + minimal class/save controls on napari."""
from pathlib import Path
import argparse,sys,json,os
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('QT_API','pyqt6')
os.environ.setdefault('QT_QPA_PLATFORM','xcb')
os.environ.setdefault('PYOPENGL_PLATFORM','glx')
os.environ.setdefault('MICROSAM_CACHEDIR',str(ROOT/'models/microsam_cache'))
from annotation_store import CLASSES,load_annotation,save_annotation
import numpy as np
import imageio.v3 as iio
import napari
from magicgui.widgets import Container,Label,ComboBox,LineEdit,CheckBox,PushButton,SpinBox

def main():
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,default=0)
    p.add_argument('--manual-only',action='store_true');p.add_argument('--smoke',action='store_true')
    a=p.parse_args();records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']
    index=a.index;reviewer='';directory=ROOT/'data/annotations'
    if a.smoke:directory=ROOT/'evidence/gui_smoke_sandbox'
    while 0<=index<len(records):
        record=records[index];image=iio.imread(ROOT/record['image'])
        meta,masks=load_annotation(directory,record)
        if a.manual_only:
            viewer=napari.Viewer();viewer.add_image(image,name='image')
        else:
            from micro_sam.sam_annotator import annotator_2d
            viewer=annotator_2d(image,model_type='vit_b',checkpoint_path=str(ROOT/'models/sam_vit_b_01ec64.pth'),
                                device='cuda',return_viewer=True,precompute_amg_state=False,
                                embedding_path=str(ROOT/'models/embeddings'/f'{record["id"]}.zarr'))
            # micro-sam 1.8.14's sync helper skips dropdown index zero. Explicitly
            # display the ORIGINAL SAM family actually loaded by this adapter.
            from micro_sam.sam_annotator._state import AnnotatorState
            AnnotatorState().widgets['embeddings'].model_family_dropdown.setCurrentText('Natural Images (SAM)')
        viewer.title=f'F10T-Vision-2 | {index+1}/{len(records)} | {record["id"]} | UNREVIEWED until confirmed'
        shape=image.shape[:2]
        layers={k:viewer.add_labels(masks[k] if masks else np.zeros(shape,np.uint32),name=k,opacity=.45) for k in CLASSES}
        track_ids=dict(meta.get('track_ids',{})) if meta else {}
        navigation={'next':None}
        status=Label(value=f'{record["id"]} | {record["split"]} | '+(meta['status'] if meta else 'UNREVIEWED'))
        category=ComboBox(label='Object class',choices=list(CLASSES),value=CLASSES[0])
        person=LineEdit(label='Human reviewer',value=reviewer)
        checked=CheckBox(label='I reviewed the full frame, including missed objects',value=False)
        notes=LineEdit(label='Notes / focus / scene',value=meta.get('notes','') if meta else '')
        track=LineEdit(label='Optional short track ID',value='')
        label_id=SpinBox(label='Instance ID for short track',min=1,max=100000,value=1)

        def select_class(*_):
            layer=layers[category.value];viewer.layers.selection.active=layer
            layer.selected_label=int(layer.data.max())+1
        category.changed.connect(select_class)

        def assign():
            if 'current_object' not in viewer.layers:status.value='Manual mode: draw with napari Labels polygon/paint tools';return
            current=viewer.layers['current_object'].data
            if not np.any(current):status.value='First segment with the official SAM Segment Object button';return
            layer=layers[category.value];data=layer.data.copy()
            for value in np.unique(current):
                if value:data[current==value]=int(data.max())+1
            layer.data=data;viewer.layers['current_object'].data=np.zeros_like(current)
            checked.value=False;status.value='Mask assigned as DRAFT. Check boundary and class; clear SAM prompts before next object.'
            select_class()

        def save(confirm=False):
            nonlocal reviewer
            try:
                if confirm and not checked.value:raise ValueError('Check full-frame review box before confirming')
                if confirm and a.smoke:raise ValueError('Software smoke test cannot create expert labels')
                reviewer=person.value
                m=save_annotation(directory,record,{k:v.data for k,v in layers.items()},reviewer,confirm,track_ids,notes.value)
                status.value=m['status']+' saved; revision '+m['revision']
                return True
            except Exception as e:status.value='NOT SAVED: '+str(e);return False

        def set_track():
            key=f'{category.value}:{label_id.value}'
            if not np.any(layers[category.value].data==label_id.value):status.value='Instance ID does not exist';return
            if track.value.strip():track_ids[key]=track.value.strip()
            else:track_ids.pop(key,None)
            checked.value=False;status.value=f'Short track metadata updated for {key}; save draft or reconfirm'

        def go(delta):
            # Navigation saves a draft only when data differ from the latest save.
            old_meta,old_masks=load_annotation(directory,record)
            changed=old_masks is None or any(not np.array_equal(old_masks[k],layers[k].data) for k in CLASSES)
            changed=changed or (old_meta is not None and (old_meta.get('track_ids',{})!=track_ids or old_meta.get('notes','')!=notes.value))
            if changed and not save(False):return
            navigation['next']=index+delta;viewer.close()

        buttons=[]
        for title,callback in [('Assign SAM current mask to class',assign),('Save draft',lambda:save(False)),
                               ('Confirm full-frame expert review',lambda:save(True)),('Set short track ID',set_track),
                               ('Previous frame',lambda:go(-1)),('Next frame',lambda:go(1))]:
            button=PushButton(text=title);button.clicked.connect(callback);buttons.append(button)
        controls=Container(widgets=[status,category,person,checked,notes,track,label_id,*buttons])
        viewer.window.add_dock_widget(controls,name='Bubble classes and versioned save',area='left')
        viewer.window._qt_window.resize(1800,1040)
        viewer.reset_view()
        viewer.bind_key('Control-S',lambda _:save(False),overwrite=True)
        viewer.bind_key('Alt-Right',lambda _:go(1),overwrite=True)
        viewer.bind_key('Alt-Left',lambda _:go(-1),overwrite=True)
        select_class()
        if a.smoke:
            from qtpy.QtCore import QTimer
            def smoke():
                try:
                    if not a.manual_only:
                        from micro_sam import util
                        state=AnnotatorState()
                        util.set_precomputed(state.predictor,state.image_embeddings)
                        sm,score,_=state.predictor.predict(point_coords=np.array([[700,500]]),point_labels=np.array([1]),multimask_output=False)
                        np.savez_compressed(ROOT/'evidence/sam_prompt_unreviewed.npz',mask=sm,score=score)
                    viewer.screenshot(path=str(ROOT/'evidence/annotation_window.png'),canvas_only=False)
                    save(False)
                    loaded,loaded_masks=load_annotation(directory,record)
                    passed=all(np.array_equal(loaded_masks[k],layers[k].data) for k in CLASSES)
                    (ROOT/'evidence/annotation_gui_smoke.json').write_text(json.dumps({'visible':viewer.window._qt_window.isVisible(),
                        'platform':os.environ.get('DISPLAY'),'manual_only':a.manual_only,'save_reopen_equal':passed,
                        'expert_confirmed':loaded['expert_confirmed'],'annotation_directory':str(directory)},indent=2))
                finally:viewer.close()
            QTimer.singleShot(5000,smoke)
        napari.run()
        if navigation['next'] is None:break
        index=navigation['next']

if __name__=='__main__':main()
