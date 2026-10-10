"""Submit this project's finite Stage4 jobs to the instance supervisor."""

import argparse
from pathlib import Path
import shlex
import subprocess

p = argparse.ArgumentParser()
p.add_argument("name")
p.add_argument("script")
p.add_argument("args", nargs="*")
a = p.parse_args()
if not a.name.replace("_", "").isalnum():
    raise ValueError("Invalid job name")
root = Path.cwd()
name = "mevo_stage4_" + a.name
cmd = [str(root / ".venv/bin/python"), "-u", str(root / a.script), *a.args]
wrapper = Path("/opt/supervisor-scripts") / (name + ".sh")
wrapper.write_text(
    "#!/bin/bash\nset -e\n. /opt/supervisor-scripts/utils/logging.sh\n. /opt/supervisor-scripts/utils/environment.sh\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexport PYTHONPATH="
    + shlex.quote(str(root))
    + "\ncd "
    + shlex.quote(str(root))
    + "\npty "
    + shlex.join(cmd)
    + " 2>&1\n"
)
wrapper.chmod(0o755)
conf = Path("/etc/supervisor/conf.d") / (name + ".conf")
conf.write_text(
    f'[program:{name}]\nenvironment=PROC_NAME="%(program_name)s"\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nstdout_logfile=/dev/stdout\nredirect_stderr=true\nstdout_logfile_maxbytes=0\n'
)
for args in [["reread"], ["update"], ["start", name]]:
    subprocess.run(["supervisorctl", *args], check=True)
