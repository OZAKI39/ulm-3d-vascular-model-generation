from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from pd_clot.lammps_fragmentation.targeted_contact_transfer import *
stage_zero()
