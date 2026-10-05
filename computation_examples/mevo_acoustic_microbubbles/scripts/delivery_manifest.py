"""Hash delivered artifacts and check local review links after all jobs close."""
import hashlib
from html.parser import HTMLParser
from pathlib import Path

from mevo_mb.io.artifacts import git_commit, write_json


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paths = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in {"src", "href"} and value and "://" not in value and not value.startswith("#"):
                self.paths.append(value)


def main():
    root = Path.cwd()
    parser = Links()
    parser.feed((root/"OPEN_RESULTS.html").read_text())
    for path in parser.paths:
        assert (root/path).is_file(), path
    skip = {".git", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache"}
    files = {}
    target = root/"docs/DELIVERY_MANIFEST.json"
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path == target or any(p in skip or p.endswith(".egg-info") for p in path.parts):
            continue
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        files[str(path.relative_to(root))] = {"bytes": path.stat().st_size, "sha256": digest}
    write_json(target, {"git_commit": git_commit(), "html_links_checked": len(parser.paths), "file_count": len(files),
                        "total_bytes": sum(v["bytes"] for v in files.values()), "files": files})
    print(f"Delivery manifest: {len(files)} files, {len(parser.paths)} valid local links", flush=True)


if __name__ == "__main__":
    main()
