"""One Bambu plate with six independent ABS objects and object-level support.

3MF object configuration follows Bambu's versioned importer:
https://github.com/bambulab/BambuStudio/blob/v02.07.01.57/src/libslic3r/Format/bbs_3mf.cpp
No G-code is written or spliced: the installed Bambu engine slices the plate.
"""
from __future__ import annotations
import json,shutil,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from .support_inputs import save,digest,verify,discover
from .compact_mold_revision import default_output as model_output
from .all_abs_pipeline import default_output as base_output
from .all_abs_slicing import profiles
from .support_slicer_adapter import slice_local
from .bambu_casting_mold import parse_paths,sample_segments

NS='http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
PROD='http://schemas.microsoft.com/3dmanufacturing/production/2015/06'
FILENAME='BG001_RMCA_ALL_ABS_X1C_02mm_single_plate.3mf'


def default_output():return model_output()/'print_plate_all_components'


def layout():
    source=model_output();summary=json.loads((source/'mold_revision_summary.json').read_text())
    if summary['status']!='READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED':raise ValueError('ACCEPTED_2MM_MODEL_REQUIRED')
    moldpath=source/'outer_mold/ABS_sacrificial_outer_mold_2mm.stl';corepath=base_output()/'bambu/vascular_reference/print_input.stl';cappath=source/'bambu/retainers/C_retainers_print_layout.stl';pkpath=source/'bambu/retainers/packing.json';records=[];meshes=[]
    def place(name,kind,m,center,support,path,part=None):
        m=m.copy();before=m.bounds.copy();delta=np.r_[center,0.]-np.r_[m.bounds.mean(0)[:2],m.bounds[0,2]];m.apply_translation(delta)
        records.append(dict(name=name,kind=kind,support_enabled=support,source=str(path),source_sha256=digest(path),source_part=part,source_to_plate_translation_mm=delta,source_bounds=before,plate_bounds=m.bounds,triangle_count=len(m.faces),volume_mm3=m.volume,scale_factor=1.));meshes.append(m)
    place('01_ABS_mold_2mm','mold',trimesh.load_mesh(moldpath),(128.,70.),False,moldpath)
    place('02_ABS_vascular_core','core',trimesh.load_mesh(corepath),(128.,181.),True,corepath)
    parts=sorted(trimesh.load_mesh(cappath).split(),key=lambda m:m.bounds[0,0]);packing=json.loads(pkpath.read_text())
    if len(parts)!=4 or len(packing)!=4:raise ValueError('EXPECTED_FOUR_PRINT_ORIENTED_C_CAPS')
    for i,(m,pk) in enumerate(zip(parts,packing)):
        place(f"{i+3:02}_C_retainer_{pk['port_id']}",'retainer',m,(219.,56.+31*i),False,cappath,pk['port_id'])
    return meshes,records,{str(p):digest(p) for p in [moldpath,corepath,cappath,pkpath,source/'mold_revision_summary.json',source/'vascular_core/BG001_RMCA_frozen_ABS_core.stl']}


def settings_dict(profiles):
    cfg={}
    for key in ['machine','filament','process']:cfg.update(profiles[key])
    cfg.update(curr_bed_type='Textured PEI Plate',printer_settings_id=profiles['machine']['name'],print_settings_id=profiles['process']['name'],filament_settings_id=[profiles['filament']['name']],print_sequence='by layer',enable_support='0')
    return single_filament_project_metadata(cfg)


def single_filament_project_metadata(settings):
    """Add project-level spool metadata absent from flattened preset files.

    The GUI counts spools via filament_colour, not filament_settings_id.
    Standard/high-flow entries are variants of ONE spool, not extra spools.
    See PresetBundle::load_config_file_config in Bambu Studio 2.7.1.57.
    """
    cfg=dict(settings)
    if len(cfg.get('filament_settings_id',[]))!=1:
        raise ValueError('EXPECTED_SINGLE_FILAMENT_PROJECT')
    variants=cfg.get('filament_extruder_variant') or ['Direct Drive Standard']
    if not isinstance(variants,list):raise ValueError('INVALID_FILAMENT_VARIANTS')
    color=cfg.get('filament_colour')
    if not color:
        color=['#91ACB6' if cfg.get('filament_type')==['PLA'] else '#C45142']
    cfg.update(filament_colour=color,filament_self_index=['1']*len(variants),
               filament_ids=[cfg.get('filament_id','')],filament_map=['1'])
    cfg.setdefault('inherits_group',['','',''])
    return cfg


def layout_audit(meshes,records):
    """Conservative XY envelopes include brim and an extra core support reserve."""
    envelopes=[];objects=[];pairs=[]
    for mesh,r in zip(meshes,records):
        margin=16. if r['support_enabled'] else 10.
        b=mesh.bounds[:,:2]+np.array([[-margin,-margin],[margin,margin]])
        envelopes.append(b)
        inside=bool((b[0]>=.11).all() and (b[1]<=255.89).all())
        excluded=bool(b[0,0]<18.11 and b[0,1]<28.11)
        objects.append(dict(name=r['name'],reserved_margin_mm=margin,reserved_bounds_xy_mm=b,passed=inside and not excluded and abs(mesh.bounds[0,2])<2e-4))
    for i,a in enumerate(envelopes):
        for j in range(i+1,len(envelopes)):
            b=envelopes[j];gap=float(np.max(np.maximum(b[0]-a[1],a[0]-b[1])))
            pairs.append(dict(first=records[i]['name'],second=records[j]['name'],separating_axis_gap_mm=gap,passed=gap>0))
    return dict(passed=len(meshes)==len(records)==6 and all(r['passed'] for r in objects+pairs),objects=objects,pairs=pairs,minimum_reserved_gap_mm=min(p['separating_axis_gap_mm'] for p in pairs),method='10 mm brim envelopes; core reserves an additional 6 mm; actual support is audited separately')


def write_project(path,meshes,records,profiles):
    root=ET.Element('model',{'unit':'millimeter','xmlns':NS});ET.SubElement(root,'metadata',{'name':'Application'}).text='BambuStudio-02.07.01.57';ET.SubElement(root,'metadata',{'name':'BambuStudio:3mfVersion'}).text='1';resources=ET.SubElement(root,'resources');build=ET.SubElement(root,'build');config=ET.Element('config')
    for j,(mesh,r) in enumerate(zip(meshes,records)):
        mid=str(j*2+1);oid=str(j*2+2);obj=ET.SubElement(resources,'object',{'id':mid,'type':'model','name':r['name']});m=ET.SubElement(obj,'mesh');vs=ET.SubElement(m,'vertices');ts=ET.SubElement(m,'triangles')
        for v in mesh.vertices:ET.SubElement(vs,'vertex',{k:format(float(x),'.12g') for k,x in zip('xyz',v)})
        for f in mesh.faces:ET.SubElement(ts,'triangle',{k:str(int(x)) for k,x in zip(['v1','v2','v3'],f)})
        parent=ET.SubElement(resources,'object',{'id':oid,'type':'model','name':r['name']});ET.SubElement(ET.SubElement(parent,'components'),'component',{'objectid':mid});ET.SubElement(build,'item',{'objectid':oid,'printable':'1'})
        c=ET.SubElement(config,'object',{'id':oid})
        for key,value in dict(name=r['name'],extruder='1',enable_support='1' if r['support_enabled'] else '0',support_type='tree(auto)',support_style='tree_organic').items():ET.SubElement(c,'metadata',{'key':key,'value':value})
        part=ET.SubElement(c,'part',{'id':mid,'subtype':'normal_part'});ET.SubElement(part,'metadata',{'key':'name','value':r['name']});ET.SubElement(part,'metadata',{'key':'matrix','value':'1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1'});ET.SubElement(part,'mesh_stat',{'face_count':str(len(mesh.faces)),'edges_fixed':'0','degenerate_facets':'0','facets_removed':'0','facets_reversed':'0','backwards_edges':'0'})
    plate=ET.SubElement(config,'plate')
    for key,value in dict(plater_id='1',plater_name='BG001 ALL ABS - 2mm mold + core + 4 C caps',locked='false',filament_map_mode='Auto For Flush',filament_maps='1',filament_volume_maps='1').items():ET.SubElement(plate,'metadata',{'key':key,'value':value})
    for j in range(len(records)):
        instance=ET.SubElement(plate,'model_instance')
        for key,val in dict(object_id=str(j*2+2),instance_id='0',identify_id=str(j+1)).items():ET.SubElement(instance,'metadata',{'key':key,'value':val})
    ET.SubElement(config,'assemble')
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml','<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>')
        z.writestr('_rels/.rels','<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>')
        z.writestr('3D/3dmodel.model',ET.tostring(root,encoding='utf-8',xml_declaration=True));z.writestr('Metadata/model_settings.config',ET.tostring(config,encoding='utf-8',xml_declaration=True));z.writestr('Metadata/project_settings.config',json.dumps(settings_dict(profiles)))


def read_objects(path):
    def matrix(s):
        m=np.eye(4)
        if s:
            a=np.fromstring(s,sep=' ').reshape(4,3);m[:3,:3]=a[:3].T;m[:3,3]=a[3]
        return m
    with zipfile.ZipFile(path) as z:
        cfg=ET.fromstring(z.read('Metadata/model_settings.config'));meta={o.attrib['id']:{r.attrib['key']:r.attrib.get('value') for r in o.findall('metadata') if 'key' in r.attrib} for o in cfg.findall('object')};cache={}
        def visit(file,oid,M):
            if file not in cache:cache[file]={o.attrib['id']:o for o in ET.fromstring(z.read(file)).findall(f'{{{NS}}}resources/{{{NS}}}object')}
            o=cache[file][oid];m=o.find(f'{{{NS}}}mesh');parts=[]
            if m is not None:
                v=np.array([[float(e.attrib[k]) for k in 'xyz'] for e in m.findall(f'{{{NS}}}vertices/{{{NS}}}vertex')]);f=np.array([[int(e.attrib[k]) for k in ['v1','v2','v3']] for e in m.findall(f'{{{NS}}}triangles/{{{NS}}}triangle')]);parts.append(trimesh.Trimesh(v@M[:3,:3].T+M[:3,3],f,process=False))
            for c in o.findall(f'{{{NS}}}components/{{{NS}}}component'):parts.extend(visit(c.attrib.get(f'{{{PROD}}}path',file).lstrip('/'),c.attrib['objectid'],M@matrix(c.attrib.get('transform'))))
            return parts
        root=ET.fromstring(z.read('3D/3dmodel.model'));items=[]
        for b in root.findall(f'{{{NS}}}build/{{{NS}}}item'):
            oid=b.attrib['objectid'];parts=visit('3D/3dmodel.model',oid,matrix(b.attrib.get('transform')));items.append(dict(object_id=oid,metadata=meta.get(oid,{}),mesh=trimesh.util.concatenate(parts)))
        return items,json.loads(z.read('Metadata/project_settings.config')),len(cfg.findall('plate')),z.testzip() is None


def validate_archive(path,meshes,records):
    objects,settings,nplates,crc=read_objects(path);rows=[]
    for mesh,r in zip(meshes,records):
        found=[o for o in objects if o['metadata'].get('name')==r['name']]
        if len(found)!=1:rows.append(dict(name=r['name'],passed=False,reason='MISSING_OR_DUPLICATE_OBJECT'));continue
        o=found[0];d=float(max(cKDTree(mesh.vertices).query(o['mesh'].vertices)[0].max(),cKDTree(o['mesh'].vertices).query(mesh.vertices)[0].max()));enabled=o['metadata'].get('enable_support',settings.get('enable_support'));expected='1' if r['support_enabled'] else '0'
        rows.append(dict(name=r['name'],object_id=o['object_id'],maximum_vertex_error_mm=d,actual_support_enabled=enabled,expected_support_enabled=expected,positive_bed_contact=abs(o['mesh'].bounds[0,2])<2e-4,passed=d<2e-4 and enabled==expected and abs(o['mesh'].bounds[0,2])<2e-4))
    result=dict(passed=len(objects)==6 and nplates==1 and crc and all(r['passed'] for r in rows) and settings.get('nozzle_diameter')==['0.2'] and settings.get('filament_type')==['ABS'] and settings.get('print_sequence')=='by layer',object_count=len(objects),plate_count=nplates,zip_crc_pass=crc,objects=rows,settings={k:settings.get(k) for k in ['printer_model','nozzle_diameter','filament_type','layer_height','enable_support','print_sequence','brim_width','brim_type']})
    return result


def slice_audit(record,meshes,records,out):
    if not record.get('slicer_success'):return dict(passed=False,reason='SLICE_FAILED')
    a=parse_paths(record['archive']);s=a['segments'];model=s[a['kind']==0];support=s[a['kind']==1];rows=[]
    for m,r in zip(meshes,records):
        center=model.mean(1);b=m.bounds;inside=np.all(center>=b[0]-.15,axis=1)&np.all(center<=b[1]+.15,axis=1);points=sample_segments(model[inside],.1);tree=cKDTree(points) if len(points) else None
        # Compare local model extrusion to its previously accepted independent
        # slice, translated to the new plate, on a 0.25 mm sampling grid.
        if r['kind']=='mold':ref=model_output()/'bambu/outer_mold/sliced.3mf'
        elif r['kind']=='core':ref=base_output()/'vascular_core/ABS_core_reference.3mf'
        else:ref=model_output()/'bambu/retainers/sliced.3mf'
        ra=parse_paths(ref);rs=ra['segments'][ra['kind']==0];refshift=np.asarray(r['source_to_plate_translation_mm']).copy()
        if r['kind']=='mold':
            M=np.asarray(json.loads((model_output()/'bambu/outer_mold/print_transform.json').read_text())['transform_4x4']);refshift-=M[:3,3]
        if r['kind']=='retainer':
            M=np.asarray(json.loads((model_output()/'bambu/retainers/print_transform.json').read_text())['transform_4x4']);bb=np.asarray(r['source_bounds'])+M[:3,3];middle=rs.mean(1);rs=rs[np.all(middle>=bb[0]-.15,axis=1)&np.all(middle<=bb[1]+.15,axis=1)];refshift-=M[:3,3]
        q=sample_segments(rs,.25)+refshift;dist=tree.query(q)[0] if tree is not None else np.array([1e9]);fraction=float(np.mean(dist<=.2));supportcenter=support.mean(1);near=np.all(supportcenter[:,:2]>=b[0,:2]-12,axis=1)&np.all(supportcenter[:,:2]<=b[1,:2]+12,axis=1)
        rows.append(dict(name=r['name'],model_segment_count=int(inside.sum()),support_segments_near_object=int(near.sum()),sampled_reference_path_coverage=fraction,reference_query_count=len(q),reference_path_tolerance_mm=.2,passed=bool(len(points) and fraction>=.99 and (int(near.sum())>0 if r['support_enabled'] else int(near.sum())==0))))
    # Audit every actual model and support endpoint against the bed/exclusion;
    # skip custom purge-line extrusion, which deliberately lies at the front.
    with zipfile.ZipFile(record['archive']) as z:
        single_plate=len([n for n in z.namelist() if n.endswith('.gcode')])==1
    deposited=s[(a['kind']==0)|(a['kind']==1)].reshape(-1,3);bed=bool((deposited[:,:2]>=.11).all() and (deposited[:,:2]<=255.89).all());excluded=bool(((deposited[:,0]<18.11)&(deposited[:,1]<28.11)).any());warnings=record['warning_audit'];r=dict(passed=all(x['passed'] for x in rows) and bed and not excluded and not warnings['critical'] and single_plate,objects=rows,model_segments=len(model),support_segments=len(support),model_and_support_inside_bed=bed,no_excluded_zone_extrusion=not excluded,warning_audit=warnings,single_gcode_plate=single_plate,method='Real native slice, object-level support retained, reference toolpath coverage within 0.2 mm; not mechanical validation')
    save(out/'actual_toolpath_audit.json',r);return r


def preview(out,meshes,records,archive=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig,ax=plt.subplots(figsize=(10,10));ax.add_patch(Rectangle((0,0),256,256,facecolor='#f4f6f7',edgecolor='#5f717f',lw=2));ax.add_patch(Rectangle((0,0),18,28,color='#c8cdd3'));colors=['#487f9a','#c95949']+['#29946d']*4
    if archive:
        p=parse_paths(archive);a=p['segments'];layer=(np.abs(a[:,:,2].mean(1)-.12)<.02);from matplotlib.collections import LineCollection
        for kind,col in [(1,'#bccbd0'),(2,'#e2c87e')]:
            seg=a[layer&(p['kind']==kind)];ax.add_collection(LineCollection(seg[:,:,:2],colors=col,linewidths=.35))
    for i,(m,r,c) in enumerate(zip(meshes,records,colors)):
        b=m.bounds;ax.add_patch(Rectangle(b[0,:2],*(b[1,:2]-b[0,:2]),fill=False,edgecolor=c,linestyle='--',lw=.7));pts=m.vertices;ax.scatter(pts[:,0],pts[:,1],s=.25,c=c);ax.text(b[:,0].mean(),b[1,1]+3,r['name'].split('_',1)[1],ha='center',fontsize=9,color=c)
    ax.set_xlim(-4,260);ax.set_ylim(-4,260);ax.set_aspect('equal');ax.set_xlabel('Build plate X (mm)');ax.set_ylabel('Build plate Y (mm)');ax.set_title('ONE PLATE | X1 Carbon / 0.2 mm / ABS\n6 independent objects | by layer | tree support for vascular core only');fig.tight_layout();fig.savefig(out/'single_plate_layout.png',dpi=170);plt.close(fig)


def run(out=None):
    out=Path(out or default_output()).resolve();out.mkdir(parents=True,exist_ok=True);meshes,records,h=layout();save(out/'source_protection_start.json',h);save(out/'plate_layout.json',dict(records=records,plate_size_mm=[256,256],layout_policy='Mold front, core rear, four individually named caps at right; inherited printing orientations; by-layer printing'))
    packing=layout_audit(meshes,records);save(out/'layout_audit.json',packing)
    if not packing['passed']:raise ValueError('PRINT_PLATE_ENVELOPES_OVERLAP_OR_OUTSIDE_BED')
    inputs=discover();values,actual,ref,engine=profiles(inputs,out/'profiles');actual['process']['print_sequence']='by layer';project=out/'input_project.3mf';write_project(project,meshes,records,actual);inputaudit=validate_archive(project,meshes,records);save(out/'input_project_audit.json',inputaudit)
    if not inputaudit['passed']:raise ValueError('INPUT_PROJECT_INVALID')
    trimesh.util.concatenate(meshes).export(out/'plate_layout_review.stl');preview(out,meshes,records);print('s1-10 plate: actual six-object native Bambu slice',flush=True);rec=slice_local(engine,actual,project,out/'bambu_all_abs_plate');print('native result',rec.get('slicer_success'),rec.get('warning_audit'),flush=True)
    archive=validate_archive(rec['archive'],meshes,records) if rec.get('archive') else dict(passed=False);paths=slice_audit(rec,meshes,records,out);protection=verify(h);protection.update(start_hashes=h,end_hashes={p:digest(p) for p in h},mismatch_count=0);save(out/'source_protection.json',protection);ready=rec.get('slicer_success',False) and archive['passed'] and paths['passed'];result=dict(status='READY_FOR_BAMBU_IMPORT_PROCESS_UNCALIBRATED' if ready else 'NEEDS_ADJUSTMENT',input_audit=inputaudit,layout_audit=packing,archive_audit=archive,toolpaths=paths,slice=rec,source_unchanged=protection['all_unchanged'],no_printer_delivery=True)
    if ready:
        dest=out/FILENAME;shutil.copy2(rec['archive'],dest);result['deliverable']=str(dest);result['sha256']=digest(dest);preview(out,meshes,records,dest)
    save(out/'single_plate_summary.json',result);print('s1-10 plate RESULT',result['status'],flush=True);return result
