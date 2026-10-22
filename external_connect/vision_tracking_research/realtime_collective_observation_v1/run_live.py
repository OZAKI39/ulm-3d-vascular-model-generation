import sys,os
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('RAYON_NUM_THREADS','4')
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from cli import main
if __name__=='__main__':main('LIVE_CAMERA',False)
