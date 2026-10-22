"""A TERM-resistant descendant is stopped even if pidfd wait APIs are unusable."""
from pathlib import Path
import subprocess
import sys
import time
import psutil

from human_vessel_flow.reference_pipe.managed_job import stop_tree


def test_kill_resistant_descendant_without_pidfd(monkeypatch):
    monkeypatch.setattr(psutil,'wait_procs',lambda *a,**k: (_ for _ in ()).throw(OSError(22,'Invalid argument')))
    code="import signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);print('ready',flush=True);time.sleep(60)"
    child_code="import subprocess,sys,time; p=subprocess.Popen([sys.executable,'-c',"+repr(code)+"],stdout=subprocess.PIPE,text=True); p.stdout.readline(); print(p.pid,flush=True); time.sleep(60)"
    proc=subprocess.Popen([sys.executable,'-c',child_code],stdout=subprocess.PIPE,text=True,start_new_session=True)
    child=psutil.Process(int(proc.stdout.readline()))
    try:
        stop_tree(proc,{},grace_seconds=.2)
        assert proc.poll() is not None
        for _ in range(20):
            if not child.is_running() or child.status()==psutil.STATUS_ZOMBIE:break
            time.sleep(.05)
        assert not child.is_running() or child.status()==psutil.STATUS_ZOMBIE
    finally:
        if proc.poll() is None:proc.kill();proc.wait()
        if child.is_running() and child.status()!=psutil.STATUS_ZOMBIE:child.kill()
