import json,os
import sv
root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
items={}
for name,cls,methods in [('PolyData',sv.modeling.PolyData,['read','set_surface','get_polydata','get_face_ids','identify_caps']),('TetGen',sv.meshing.TetGen,['load_model','set_model','set_walls','get_model_face_ids','generate_mesh','get_mesh','get_surface','get_face_polydata','write_mesh'])]:
 for method in methods:
  obj=getattr(cls,method,None);items[name+'.'+method]=getattr(obj,'__doc__','NOT_AVAILABLE')
with open(os.path.join(root,'outputs/sv1/api_documentation.json'),'w') as f:json.dump(items,f,indent=2)
