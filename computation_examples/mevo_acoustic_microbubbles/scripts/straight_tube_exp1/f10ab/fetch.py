"""Recover result files without deleting local or historical files."""
import subprocess,sys
from .common import BASE
from .remote import SERVER

def main():
 filters=['--include=*/','--include=INDEX.json','--include=COMPLETE.json','--include=SUMMARY.json','--include=*AUDIT.json','--exclude=*'] if '--metadata' in sys.argv else []
 subprocess.run(['rsync','-az','--partial',*filters,'-e','ssh -i /home/lzy/.ssh/vast_step3b_ed25519 -p 15283',
  'root@95.3.33.46:'+SERVER+'/results/straight_tube_exp1/f10ab/',str(BASE)+'/'],check=True)
if __name__=='__main__':main()
