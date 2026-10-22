"""Resumeable V3.1 study, invoked with the project environment."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vascular_processing.support_v3_1_pipeline import run
if __name__ == '__main__':
    run()
