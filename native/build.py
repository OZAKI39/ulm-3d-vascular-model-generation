"""Build the pinned native dependency from source, without historical build objects."""

from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--wrapper-only", action="store_true")
    args = p.parse_args()
    archive = ROOT / "vendor/lammps_core_20260930.tar.gz"
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    pin = json.loads((ROOT / "vendor/MANIFEST.json").read_text())
    if digest(archive) != pin["sha256"]:
        raise ValueError("Pinned LAMMPS source archive checksum mismatch")
    work = ROOT / "work"
    source = work / "source"
    build = work / "cmake"
    source.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        tar.extractall(source, filter="data")
    # One authoritative copy of the custom kernels lives in native/src.
    for f in (ROOT / "src").iterdir():
        if f.name != "mapped_geometry.cpp":
            shutil.copyfile(f, source / "src/BPM" / f.name)
    # Install the matching upstream wrapper as a normal third-party package.
    environment = dict(os.environ, LAMMPS_VERSION_FILE=str(source / "src/version.h"))
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "--no-deps", str(source / "python")],
        check=True,
        env=environment,
    )
    if args.wrapper_only:
        return
    build.mkdir(parents=True, exist_ok=True)
    commands = [
        [
            "cmake",
            "-S",
            str(source / "cmake"),
            "-B",
            str(build),
            "-DCMAKE_BUILD_TYPE=Release",
            "-DBUILD_SHARED_LIBS=ON",
            "-DBUILD_MPI=ON",
            "-DBUILD_OMP=ON",
            "-DPKG_BPM=ON",
            "-DPKG_MOLECULE=ON",
        ],
        ["cmake", "--build", str(build), "-j", str(args.jobs)],
    ]
    for cmd in commands:
        subprocess.run(cmd, check=True)
    destination = ROOT / "build"
    destination.mkdir(exist_ok=True)
    shutil.copyfile((build / "liblammps.so").resolve(), destination / "liblammps.so")
    mapped = destination / "libmapped_geometry.so"
    command = [
        "g++",
        "-O3",
        "-std=c++17",
        "-fPIC",
        "-shared",
        str(ROOT / "src/mapped_geometry.cpp"),
        "-o",
        str(mapped),
    ]
    subprocess.run(command, check=True)
    commands.append(command)
    mapped.with_suffix(".source.sha256").write_text(
        digest(ROOT / "src/mapped_geometry.cpp") + "\n"
    )
    report = dict(
        archive_sha256=digest(archive),
        source_sha256={f.name: digest(f) for f in (ROOT / "src").iterdir()},
        library_sha256=digest(destination / "liblammps.so"),
        mapped_geometry_sha256=digest(mapped),
        commands=commands,
    )
    (destination / "BUILD.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
