#!/usr/bin/env python3
"""Retrieve only this stage's sanitized native evidence (no historical writers)."""
import subprocess,shlex
from remote import ROOT,REMOTE,ssh_prefix
ssh=ssh_prefix()
for remote,local in [('reports','reports/sv1_3q/remote'),('logs','logs/sv1_3q/remote')]:
    (ROOT/local).mkdir(parents=True,exist_ok=True)
    subprocess.run(['rsync','-az','-e',shlex.join(ssh[:-1]),ssh[-1]+':'+REMOTE+'/'+remote+'/',str(ROOT/local)+'/'],check=True)
