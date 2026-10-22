#!/usr/bin/env python3
"""Publish bounded commits to explicitly named new branches; never force push."""
import argparse, hashlib, json, os, subprocess, time
from pathlib import Path
from export_flow_projects import BASE,AUDIT

ap=argparse.ArgumentParser();ap.add_argument('project',choices=['brava','mouse']);ap.add_argument('--push',action='store_true');args=ap.parse_args()
root=BASE/'temp_storage'/f'{args.project}_flow_github_sync_20261022'
branch=f'sync/{"brava-flow" if args.project=="brava" else "mouse-vascular-flow"}-project-20261022'
def git(*cmd,input=None):
 return subprocess.run(['git','-C',str(root),'-c','gc.auto=0','-c','diff.renames=false','-c','credential.helper=','-c','credential.helper=!gh auth git-credential',*cmd],input=input,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout

assert git('branch','--show-current').decode().strip() in {branch,branch+'-publish'}
validation=json.loads((AUDIT/f'{args.project.upper()}_VALIDATION.json').read_text())
assert validation['status']=='PASS',validation
planfile=AUDIT/f'{args.project.upper()}_PUBLISH_PLAN.json'
if not planfile.exists():
 parent=git('rev-parse','HEAD').decode().strip()
 old={line.split()[2] for line in git('ls-tree','-r',parent).decode().splitlines()}
 files=[];links=[]
 for d,ds,fs in os.walk(root,followlinks=False):
  ds[:]=[n for n in ds if n!='.git']
  for n in list(ds):
   p=Path(d)/n
   if p.is_symlink():links.append(p);ds.remove(n)
  for n in fs:
   p=Path(d)/n
   if p.name=='.git':continue
   (links if p.is_symlink() else files).append(p)
 files.sort()
 names=[str(p.relative_to(root)) for p in files]
 assert not any('\n' in n for n in names)
 print('Hashing exact Git blobs',len(names),flush=True)
 oids=git('hash-object','-w','--no-filters','--stdin-paths',input=('\n'.join(names)+'\n').encode()).decode().splitlines()
 entries=[]
 for p,name,oid in zip(files,names,oids,strict=True):
  entries.append({'path':name,'oid':oid,'mode':'100755' if p.stat().st_mode&0o111 else '100644','bytes':p.stat().st_size})
 for p in links:
  data=os.readlink(p).encode();oid=git('hash-object','-w','--stdin',input=data).decode().strip()
  entries.append({'path':str(p.relative_to(root)),'oid':oid,'mode':'120000','bytes':len(data)})
 known=[e for e in entries if e['oid'] in old];new=[e for e in entries if e['oid'] not in old]
 # Existing objects are already reachable from the advertised parent branch.
 # Cap newly introduced logical bytes well below a single GitHub receive pack.
 new.sort(key=lambda e:(not(e['path']=='README.md' or e['path'].startswith('sync_20261022/')),e['bytes']>1_000_000,e['path']))
 batches=[known];size=0;seen=set(old)
 for e in new:
  extra=0 if e['oid'] in seen else e['bytes']
  if size+extra>550_000_000 and batches[-1]:batches.append([]);size=0
  batches[-1].append(e);size+=extra;seen.add(e['oid'])
 plan={'project':args.project,'branch':branch,'parent':parent,'files':len(entries),'bytes':sum(e['bytes'] for e in entries),'batches':batches}
 planfile.write_text(json.dumps(plan,indent=2)+'\n')
 print('Plan ready',len(batches),'bounded commits',len(known),'already published blobs',flush=True)
else:plan=json.loads(planfile.read_text())
if not args.push:raise SystemExit(0)

receiptfile=AUDIT/f'{args.project.upper()}_PUSH_RECEIPT.json'
receipt=json.loads(receiptfile.read_text()) if receiptfile.exists() else {'project':args.project,'branch':branch,'parent':plan['parent'],'commits':[],'status':'IN_PROGRESS'}
if not receipt['commits']:
 assert git('rev-parse','HEAD').decode().strip()==plan['parent']
 assert not git('ls-remote','origin','refs/heads/'+branch).strip(),'New remote branch unexpectedly exists'
 git('read-tree','--empty')
for i,batch in enumerate(plan['batches']):
 if i<len(receipt['commits']) and receipt['commits'][i].get('pushed'):continue
 if i==len(receipt['commits']):
  payload=b''.join(f"{e['mode']} {e['oid']}\t{e['path']}\0".encode() for e in batch)
  git('update-index','-z','--index-info',input=payload)
  message=AUDIT/f'{args.project}_commit_message.txt'
  message.write_text(f"Sync {args.project} flow project separately ({i+1}/{len(plan['batches'])})\n\nPreserve existing scientific results and raw accepted states; exclude rebuildable environments, superseded bulk outputs and redundant intermediates. No scientific computation or numerical data changes.\n")
  git('commit','--quiet','-F',str(message))
  commit=git('rev-parse','HEAD').decode().strip()
  receipt['commits'].append({'batch':i+1,'sha':commit,'files_added':len(batch),'pushed':False});receiptfile.write_text(json.dumps(receipt,indent=2)+'\n')
 log=AUDIT/f'{args.project}_push_{i+1:02d}.log'
 print('Push',args.project,i+1,'/',len(plan['batches']),receipt['commits'][i]['sha'],flush=True)
 with log.open('ab') as f:
  proc=subprocess.run(['git','-C',str(root),'-c','gc.auto=0','-c','credential.helper=','-c','credential.helper=!gh auth git-credential','-c','pack.threads=2','-c','pack.window=0','-c','pack.windowMemory=64m','-c','pack.deltaCacheSize=64m','push','--porcelain','origin',f'HEAD:refs/heads/{branch}'],stdout=f,stderr=subprocess.STDOUT)
 if proc.returncode:raise SystemExit(f'push failed; inspect {log}')
 remote=git('ls-remote','origin','refs/heads/'+branch).decode().split()[0]
 assert remote==receipt['commits'][i]['sha']
 receipt['commits'][i]['pushed']=True;receiptfile.write_text(json.dumps(receipt,indent=2)+'\n')
receipt['status']='PUSHED_AND_REMOTE_SHA_VERIFIED';receipt['head']=git('rev-parse','HEAD').decode().strip();receipt['completed_epoch']=time.time()
receiptfile.write_text(json.dumps(receipt,indent=2)+'\n')
git('fetch','origin',f'refs/heads/{branch}:refs/remotes/origin/{branch}')
git('branch','--set-upstream-to=origin/'+branch)
print('COMPLETE',args.project,receipt['head'],flush=True)
