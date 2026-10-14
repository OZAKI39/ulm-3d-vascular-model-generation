"""Return only F10O artifacts; hash every remote file before and after transfer."""

import subprocess
import json
from .common import BASE, write, sha
from .remote import SERVER


def main():
    remote = SERVER + "/results/straight_tube_exp1/f10o"
    code = """import pathlib,hashlib,json
p=pathlib.Path(ROOT)
files={str(f.relative_to(p)):dict(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted(p.rglob('*')) if f.is_file()}
print(json.dumps(dict(files=files,bytes=sum(r['bytes'] for r in files.values()))))
""".replace("ROOT", repr(remote))
    data = json.loads(subprocess.check_output(["ssh", "mevo-new", "python3 -"], input=code, text=True))
    subprocess.run(["rsync", "-az", "mevo-new:" + remote + "/", str(BASE) + "/"], check=True)
    for rel, rec in data["files"].items():
        assert sha(BASE / rel) == rec["sha256"], rel
    # Append snapshots, never overwrite a previous transfer ledger.
    folder = BASE / "transfers"
    folder.mkdir(exist_ok=True)
    write(folder / f"{len(list(folder.glob('*.json'))) + 1:02d}.json", data)
    write(BASE / "REMOTE_MANIFEST.json", data)
    print("Verified", len(data["files"]), "files,", data["bytes"], "bytes")


if __name__ == "__main__":
    main()
