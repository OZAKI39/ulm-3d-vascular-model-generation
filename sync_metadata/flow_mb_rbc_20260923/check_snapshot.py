"""Read-only packaging verification; does not integrate particles or rerun CFD."""
from pathlib import Path
import argparse,json,hashlib,re,subprocess,collections

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repository',type=Path,required=True);ap.add_argument('--inventory',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 root=a.repository.resolve();rows=json.loads(a.inventory.read_text());selected=[r for r in rows if r['status']=='INCLUDED'];fail=[]
 for r in selected:
  p=root/r['target']
  if not p.is_file() or sha(p)!=r['source_sha256']:fail.append(r['target'])
 base=subprocess.check_output(['git','ls-tree','-r','--name-only','c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2'],cwd=root,text=True).splitlines()
 missing_base=[p for p in base if not(root/p).is_file()]
 paths=set(base)|{r['target'] for r in selected}
 patterns={
  'private_key':rb'-----BEGIN (?:RSA |EC |OPENSSH |DSA |ENCRYPTED )?PRIVATE KEY-----',
  'github_token':rb'\b(?:gh[pousr]_[A-Za-z0-9_]{30,}|github_pat_[A-Za-z0-9_]{30,})\b',
  'aws_access_key':rb'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b',
  'openai_or_other_sk_key':rb'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_\-]{25,}\b',
  'slack_token':rb'\bxox[baprs]-[A-Za-z0-9-]{20,}\b',
  'credential_in_url':rb'https?://[^\s/:@]{1,80}:[^\s/@]{6,80}@',
 }
 broad=re.compile(rb'''(?i)["']?(?:api[_-]?key|access[_-]?token|auth[_-]?token|password|secret[_-]?key|vast[_-]?api[_-]?key)["']?\s*[:=]\s*["']([^"'\r\n]{8,180})["']''')
 hits=[];assignments=[];scanned=0
 binary={'.npz','.npy','.vtu','.vtp','.stl','.h5','.bin','.png','.jpg','.jpeg','.gif','.mp4','.pdf','.gz','.pptx','.xlsx','.ico','.woff','.ttf','.msh'}
 for name in sorted(paths):
  p=root/name
  if not p.is_file() or p.suffix.lower() in binary:continue
  data=p.read_bytes()
  if b'\x00' in data[:8192]:continue
  scanned+=1
  for label,pattern in patterns.items():
   for m in re.finditer(pattern,data):hits.append({'path':name,'line':data[:m.start()].count(b'\n')+1,'type':label,'value_sha256':hashlib.sha256(m.group()).hexdigest()})
  for m in broad.finditer(data):
   val=m.group(1)
   if val.lower() in {b'redacted',b'<redacted>',b'********',b'password',b'not-used',b'not_used',b'undefined',b'not configured'}:continue
   assignments.append({'path':name,'line':data[:m.start()].count(b'\n')+1,'value_sha256':hashlib.sha256(val).hexdigest(),'value_length':len(val)})
 oversized=[r['target'] for r in selected if r['target'] not in base and r['bytes']>25*1024**2]
 report={'status':'PASS' if not(fail or missing_base or hits or assignments or oversized) else 'REVIEW_REQUIRED','selected_files':len(selected),'verified_source_destination_hashes':len(selected)-len(fail),'hash_failures':fail,'missing_published_base_files':missing_base,'new_files_over_25_MiB':oversized,'scanned_text_files':scanned,'credential_pattern_hits':hits,'credential_assignment_review':assignments}
 a.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
