import subprocess,time
from .common import *
def main():
 t=time.perf_counter();p=subprocess.run(['.venv/bin/python','-m','pytest','-q','tests/f10aa','tests/f10ab','tests/f10ac'],capture_output=True,text=True)
 write(BASE/'evidence/TESTS.json',dict(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,wall_s=time.perf_counter()-t,identity=identity()));print(p.stdout);assert p.returncode==0
if __name__=='__main__':main()
