import time,subprocess
from .common import *
def main():
    rows=[];start=time.perf_counter()
    while time.perf_counter()-start<600:
        s=subprocess.check_output(['nvidia-smi','--query-gpu=timestamp,memory.used,utilization.gpu,power.draw','--format=csv,noheader,nounits'],text=True).strip();rows.append(dict(elapsed_s=time.perf_counter()-start,nvidia=s));write(BASE/'RESOURCE_SAMPLES.json',dict(samples=rows,scope='Whole-GPU sampled usage, not exclusive per-process peak',identity=identity()));time.sleep(2)
if __name__=='__main__':main()
