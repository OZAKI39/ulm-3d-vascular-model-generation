"""Executed on server: hashes and reproducible references, no source redistribution."""
import sys,json,subprocess,importlib.metadata as md
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
def main():
    server=Path('/workspace/f10w_vision_20261017');models={}
    for name in ['LRAFT-ST','co-tracker','UnLiteFlowNet-PIV']:
        repo=THIRD/name;licenses=[p for p in repo.iterdir() if p.is_file() and p.name.lower().startswith(('license','copying'))]
        models[name]=dict(url=subprocess.check_output(['git','remote','get-url','origin'],cwd=repo,text=True).strip(),commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),status=subprocess.check_output(['git','status','--short'],cwd=repo,text=True),licenses={p.name:dict(sha256=sha(p),text=p.read_text()) for p in licenses},license_status='FOUND_READ_ATTACHED' if licenses else 'NO_LICENSE_FOUND_DO_NOT_REDISTRIBUTE')
    files=[]
    for root in [server/'models',server/'runs/models_v1/data',THIRD/'UnLiteFlowNet-PIV/models']:
        for p in root.rglob('*'):
            if p.is_file():files.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
    packages={}
    for n in ['torch','torchvision','numpy','scipy','opencv-python','openpiv','POT','pyarrow','av','cupy-cuda12x','flowiz']:
        try:packages[n]=md.version(n)
        except md.PackageNotFoundError:packages[n]='NOT_INSTALLED'
    lock=subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True);(server/'setup/requirements_final.lock').write_text(lock)
    dump(OUT/'evidence/remote_assets.json',dict(identity=identity(),third_party=models,files=files,packages=packages,lock_path=str(server/'setup/requirements_final.lock'),lock_sha256=sha(server/'setup/requirements_final.lock'),restoration='SSH read-only copies from listed paths, verify SHA-256. Clone official repos at commits; apply saved UnLite compatibility patch. Do not publicly redistribute unclear-license sources or weights.'))
if __name__=='__main__':main()
