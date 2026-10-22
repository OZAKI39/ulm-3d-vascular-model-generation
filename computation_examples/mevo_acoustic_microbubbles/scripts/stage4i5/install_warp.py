"""Install one pinned optional wheel into an isolated target, never the main venv."""
import importlib.metadata
import importlib.util
import subprocess
import sys

from scripts.stage4i5.common import ROOT,new_run,finish,write


def environment():
    return {d.metadata["Name"]:d.version for d in importlib.metadata.distributions()}


def main():
    out=new_run("p2_warp",dict(operation="optional_dependency_install",version="1.17.0"))
    target=ROOT/".optional/warp-1.17.0"
    if target.exists():
        raise FileExistsError("Optional install destination already exists; audit before reuse")
    before=environment()
    write(out/"MAIN_ENVIRONMENT_BEFORE.json",before)
    assert importlib.util.find_spec("warp") is None
    command=["uv","pip","install","--python",sys.executable,"--target",str(target),"--no-deps","warp-lang==1.17.0"]
    run=subprocess.run(command,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/"INSTALL.log").write_text(run.stdout)
    write(out/"INSTALL_COMMAND.json",command)
    after=environment()
    write(out/"MAIN_ENVIRONMENT_AFTER.json",after)
    assert before==after and importlib.util.find_spec("warp") is None
    if run.returncode:
        finish(out,dict(status="FAIL",exit_code=run.returncode,main_environment_unchanged=True))
        raise RuntimeError("Optional Warp installation failed; primary environment unchanged")
    gpu=subprocess.check_output(["nvidia-smi","--query-gpu=name,driver_version,compute_cap","--format=csv,noheader"],text=True).strip()
    finish(out,dict(status="INSTALLED_OPTIONAL_ONLY",version="1.17.0",target=str(target.relative_to(ROOT)),
        main_environment_unchanged=True,environment_diff={},GPU=gpu,
        wheel_CUDA="Official 1.17.0 PyPI wheel built with CUDA 12.9; device initialization audited at benchmark",
        license="Apache-2.0",license_source="https://github.com/NVIDIA/warp/blob/v1.17.0/LICENSE.md",
        release_source="https://github.com/NVIDIA/warp/releases/tag/v1.17.0",install_command=command,
        activation="Explicit process-local PYTHONPATH overlay only; no primary-environment Warp dependency"))


if __name__=="__main__":
    main()
