"""Restricted GPU/CPU snapshots without process arguments or environment secrets."""
import subprocess,time,json,os
from .common import *
def main():
 d=dict(utc=time.time(),load=os.getloadavg(),cpu_count=os.cpu_count(),gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.used,utilization.gpu,power.draw','--format=csv,noheader'],text=True),memory=Path('/proc/meminfo').read_text(),disk=subprocess.check_output(['df','-h',str(ROOT)],text=True),identity=identity())
 write(BASE/'evidence'/('RESOURCES_'+str(time.time_ns())+'.json'),d)
if __name__=='__main__':main()
