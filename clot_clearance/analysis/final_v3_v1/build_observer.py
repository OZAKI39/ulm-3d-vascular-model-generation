from common import *
import subprocess,time
cmd=['mpicxx','-std=c++17','-O2','-shared','-fPIC','-DLAMMPS_SMALLBIG','-DOMPI_SKIP_MPICXX','-I'+str(S/'src'),'-I'+str(S/'third_party'),str(Path(__file__).with_name('observer.cpp')),'-L'+str(B),'-llammps','-o',str(O/'observer.so')]
t=time.perf_counter();p=subprocess.run(cmd,capture_output=True,text=True);dump(O/'OBSERVER_BUILD.json',{'argv':cmd,'exit_code':p.returncode,'wall_s':time.perf_counter()-t,'stdout':p.stdout,'stderr':p.stderr,'native_library_recompiled':False});print(p.stdout,p.stderr);assert p.returncode==0
