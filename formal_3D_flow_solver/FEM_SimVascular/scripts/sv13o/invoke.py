"""Run a WSL-authored stage script natively on the authorized build host."""
import shlex,subprocess,sys
from remote import ROOT,REMOTE,upload,ssh_prefix
script=sys.argv[1];assert '/' not in script and script.endswith('.py')
upload(ROOT/'scripts/sv13o'/script,script)
raise SystemExit(subprocess.call(ssh_prefix()+[shlex.join(['/usr/bin/python3','-B',REMOTE+'/'+script,*sys.argv[2:]])]))
