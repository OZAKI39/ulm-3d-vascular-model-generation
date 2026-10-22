"""Small permanent geometry fixtures for independent V2 failure-path tests."""
import json
import numpy as np
import trimesh
from .support_inputs import ROOT,save
from .casting_fixture_v2_audit import manifold,mesh

CASES=['single_port_inside_wall','two_ports_same_wall','ports_three_walls','curved_extension_rejected','short_extension_rejected','support_blade_collision','central_tree_near_wall','local_notch_trim','side_strip_trim','trim_plane_hits_wrong_branch']


def cylinder(a,b,r=.7):
    a=np.array(a,float);b=np.array(b,float);m=trimesh.creation.cylinder(radius=r,height=np.linalg.norm(b-a),sections=32);t=trimesh.geometry.align_vectors([0,0,1],(b-a)/np.linalg.norm(b-a));t[:3,3]=(a+b)/2;m.apply_transform(t);return m


def create(case):
    central=trimesh.creation.icosphere(subdivisions=2,radius=5);central.apply_translation([0,0,10]);configs=[('I1',np.array([1.,0,0]),np.array([0.,0,10]))]
    if case in ['two_ports_same_wall','side_strip_trim','local_notch_trim']:
        configs=[('O2',np.array([1.,0,0]),np.array([0.,-3,10])),('O3',np.array([1.,0,0]),np.array([0.,3,10]))]
    if case=='ports_three_walls':configs += [('O1',np.array([-1.,0,0]),np.array([0.,0,10])),('O2',np.array([0.,-1,0]),np.array([0.,0,10]))]
    ports=[];source_ports={};full=manifold(central)
    for j,(k,n,root) in enumerate(configs):
        a=root+n*10;b=root+n*19;r=.7;pm=cylinder(root,b,r);source_ports[k]=pm;full=full+manifold(pm)
        face=('+' if n[np.argmax(abs(n))]>0 else '-')+'XY'[np.argmax(abs(n))]
        ports.append(dict(port_id=k,face=face,axis=n,straight_start=a,end=b,straight_length_mm=9.,radius_mm=r,diameter_mm=2*r,attachment=root,parent_branch=k,swc_id=j+1,total_route_length_mm=19.,trim_origin=(a+b)/2,trim_normal=n,axis_plane_angle_degrees=90.,extension_wall_setback_mm=1.5,role='MASTER' if j==0 else 'SECONDARY',cradle_clearance_mm=.22 if j==0 else .33,anti_float_gap_mm=.22,axial_float_mm=0. if j==0 else 1.,usable_straight_length_after_trim_mm=4.5,distance_to_bifurcation_mm=14.5,distance_to_cap_attachment_mm=14.5,distance_to_curvature_mm=4.5,approved_contact_start=a+n*.8,approved_contact_end=b-n*.8))
    if case in ['support_blade_collision','central_tree_near_wall','trim_plane_hits_wrong_branch']:
        if case=='support_blade_collision':bad=cylinder([16.7,-2,8],[16.7,2,8],.8)
        elif case=='central_tree_near_wall':bad=cylinder([0,0,10],[20,2,10],.8)
        else:bad=cylinder([14.5,-5,12],[14.5,5,12],.5)
        central=mesh(manifold(central)+manifold(bad));full=full+manifold(bad)
    core=mesh(full);lo=core.bounds[0]-[1.5,3,0];hi=core.bounds[1]+[1.5,3,12];lo[2]=0;hi[2]=27
    if case in ['two_ports_same_wall','side_strip_trim']:lo[1]=-7.5;hi[1]=7.5
    if case=='local_notch_trim':lo[1]=-15;hi[1]=15
    for p in ports:
        n=p['axis'];ax=np.argmax(abs(n));bound=hi[ax] if n[ax]>0 else lo[ax];p['wall_distance_from_trim_mm']=float((bound-p['trim_origin'][ax])/n[ax])
    raw=[]
    for p in ports:
        r=dict(port_id=p['port_id'],face=p['face'],direction=p['axis'],target=p['straight_start'].copy(),end=p['end'].copy(),radius_mm=p['radius_mm'],accepted_attachment={'calculation_origin_mm':p['attachment']},parent_branch=p['parent_branch'],swc_id=p['swc_id'],length_mm=p['total_route_length_mm'])
        if case=='curved_extension_rejected':r['target'][1]+=1
        if case=='short_extension_rejected':r['end']=r['target']+r['direction']*4
        raw.append(r)
    return dict(core=core,central=central,source_ports=source_ports,ports=ports,low=lo,high=hi,transform=np.eye(4),master_port=ports[0]['port_id'],bifurcations=np.array([[0.,0,10]])),raw


def generate():
    root=ROOT/'tests/data/casting_fixture/v2'
    for case in CASES:
        l,raw=create(case);p=root/case;p.mkdir(parents=True,exist_ok=True);l['core'].export(p/'core.stl');l['central'].export(p/'central.stl')
        for k,v in l['source_ports'].items():v.export(p/(k+'.stl'))
        save(p/'fixture.json',dict(synthetic=True,case=case,raw_ports=raw,layout={k:v for k,v in l.items() if k not in ['core','central','source_ports']}))


def load(case):
    p=ROOT/'tests/data/casting_fixture/v2'/case;data=json.loads((p/'fixture.json').read_text());l=data['layout']
    for k in ['low','high','transform','bifurcations']:l[k]=np.array(l[k])
    for port in l['ports']:
        for k in ['axis','straight_start','end','attachment','trim_origin','trim_normal']:port[k]=np.array(port[k])
    l.update(core=trimesh.load_mesh(p/'core.stl'),central=trimesh.load_mesh(p/'central.stl'),source_ports={x['port_id']:trimesh.load_mesh(p/(x['port_id']+'.stl')) for x in l['ports']})
    return l,data['raw_ports']

if __name__=='__main__':generate()
