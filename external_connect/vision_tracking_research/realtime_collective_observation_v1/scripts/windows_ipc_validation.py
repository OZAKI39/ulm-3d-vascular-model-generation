import sys,time,subprocess,threading
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,stats
from live_adapter import ConfiguredWindowsBridge,winpath
b=ConfiguredWindowsBridge(ROOT/'results/windows_ipc_mock','/mnt/f/CHRONOS_control','/mnt/f/CHRONOS_control/.venv/Scripts/python.exe',10)
try:
 b.log=(b.output/'stderr.log').open('wb');b.process=subprocess.Popen([b.python,'-B',winpath(ROOT/'scripts/windows_mock_server.py')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=b.log,bufsize=0);b.thread=threading.Thread(target=b._loop,daemon=True);start=time.monotonic();b.thread.start();time.sleep(5);b.stop();span=time.monotonic()-start;rows=[r for r in b.records if r['status']=='FRAME'];rtt=[r['get_roundtrip_ms'] for r in rows];dump(b.output/'summary.json',dict(real_windows_process=True,real_wsl_pipe=True,camera_access=False,synthetic_only=True,received=len(rows),resolution=[1280,1024],payload_bytes_per_frame=1280*1024*3,wall_s=span,throughput_fps=len(rows)/span,throughput_mib_s=len(rows)*1280*1024*3/span/2**20,rtt_ms=stats(rtt),exit_code=b.process.returncode,closed=not b.thread.is_alive(),compression=False,pass_=len(rows)>0 and b.process.returncode==0 and not b.thread.is_alive()))
except Exception as e:dump(b.output/'summary.json',dict(pass_=False,error=repr(e),camera_access=False))
