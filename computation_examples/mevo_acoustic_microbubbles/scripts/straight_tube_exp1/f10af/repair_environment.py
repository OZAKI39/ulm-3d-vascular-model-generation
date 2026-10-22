"""Relocate only this copied environment's textual launchers; no dependency change."""
from pathlib import Path
import hashlib,json,re,sys

def main():
    env=Path(sys.argv[1]).resolve();dest=Path(sys.argv[2]);text=(env/'bin/activate').read_text()
    old=re.search(r"^\s*(?:export )?VIRTUAL_ENV=['\"]?(/[^'\"\s]+)",text,re.M).group(1)
    assert env.name=='.venv' and 'f10af' in str(env)
    changed=[]
    for p in (env/'bin').iterdir():
        if not p.is_file() or p.is_symlink() or p.stat().st_size>1024*1024:continue
        try:s=p.read_text()
        except UnicodeDecodeError:continue
        if old==str(env) or old not in s:continue
        before=hashlib.sha256(p.read_bytes()).hexdigest();p.write_text(s.replace(old,str(env)))
        changed.append(dict(path=str(p),before_sha256=before,after_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(dict(previous_environment=old,environment=str(env),changed=changed,scope='Activation and executable text launchers only. Scientific jobs used explicit AF Python plus PYTHONPATH=src:. and were unaffected.'),indent=2)+'\n')
if __name__=='__main__':main()
