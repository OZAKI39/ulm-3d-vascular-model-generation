"""Correct derived endpoint velocity on rejected-step prefixes; raw data untouched."""
import shutil
from .common import *
from .audit import one

def main():
    old=read(BASE/'AUDIT.json');runs=read(BASE/'ALL_HF_RUNS.json');todo=[r for r in runs if r['status']!='PASS' and r.get('accepted_steps',0)>0];folder=BASE/'evidence/derived_view_v1';folder.mkdir(parents=True,exist_ok=False);notes=[]
    shutil.copy2(BASE/'AUDIT.json',folder/'AUDIT.json')
    for r in todo:
        p=BASE/'dataset'/(r['id']+'.h5');shutil.copy2(p,folder/p.name);oldsha=sha(p);new=one(r);old['rows']=[new if z['id']==r['id'] else z for z in old['rows']];notes.append(dict(id=r['id'],old_sha256=oldsha,new_sha256=sha(p),reason='Endpoint velocity now references last accepted step, never rejected attempted step; actual positions unchanged'))
    old['derived_view_repair']=notes;write(BASE/'AUDIT.json',old);write(BASE/'DERIVED_VIEW_REPAIR.json',dict(rows=notes,identity=identity()))
if __name__=='__main__':main()
