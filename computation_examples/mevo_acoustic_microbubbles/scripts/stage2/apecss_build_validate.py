"""Build unmodified APECSS v1.8 and execute relevant official validation examples."""

from pathlib import Path
import os
import shutil
import subprocess
import sys
import time

import numpy as np

from mevo_mb.io.artifacts import write_json


def main():
    root = Path.cwd()
    upstream = root / "third_party/src/apecss"
    build = root / "third_party/build/apecss"
    out = root / "results/stage2/apecss_validation" / time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    out.mkdir(parents=True, exist_ok=False)
    build.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, APECSS_DIR=str(upstream), USRLIB_DIR="/usr/lib/x86_64-linux-gnu", MPLBACKEND="Agg")
    records = []

    def run(args, directory, name, timeout=180):
        start = time.perf_counter()
        try:
            with (directory / (name + ".stdout.log")).open("w") as log:
                result = subprocess.run(
                    [str(x) for x in args],
                    cwd=directory,
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    check=False,
                )
            code, error = result.returncode, None
        except subprocess.TimeoutExpired as exc:
            code, error = 124, str(exc)
        item = {
            "command": [str(x) for x in args],
            "cwd": str(directory),
            "returncode": code,
            "elapsed_s": time.perf_counter() - start,
            "error": error,
        }
        records.append(item)
        write_json(out / "COMMANDS.json", {"commands": records})
        print(name, code, item["elapsed_s"], flush=True)
        return code

    commit = subprocess.check_output(
        ["git", "-c", f"safe.directory={upstream}", "-C", str(upstream), "rev-parse", "HEAD"], text=True
    ).strip()
    assert commit == "a0ee071ca247b5489e997d90e83b5f53c40e4ec7"
    assert (
        run(
            [
                "cmake",
                "-S",
                upstream / "lib",
                "-B",
                build,
                "-DCMAKE_BUILD_TYPE=Release",
                "-DCMAKE_C_FLAGS=-march=native",
            ],
            out,
            "configure",
        )
        == 0
    )
    assert run(["cmake", "--build", build, "--parallel", "6", "--verbose"], out, "build") == 0
    for example, stem in [("ultrasound", "ultrasound"), ("rayleighcollapse", "rayleigh")]:
        source = upstream / "examples" / example / "src" / (example + "_apecss.c")
        assert (
            run(
                [
                    "gcc",
                    "-O3",
                    "-march=native",
                    "-DNDEBUG",
                    "-I",
                    upstream / "include",
                    source,
                    build / "libapecss.a",
                    "-lm",
                    "-o",
                    build / (stem + "_apecss"),
                ],
                out,
                "build_" + stem,
            )
            == 0
        )
    probe = build / "precision.c"
    probe.write_text(
        '#include <stdio.h>\n#include "apecss.h"\nint main(void){printf("%zu\\n",sizeof(APECSS_FLOAT));return 0;}\n'
    )
    assert (
        run(["gcc", "-I", upstream / "include", probe, "-o", build / "precision"], out, "precision_build")
        == 0
    )
    precision = int(subprocess.check_output([str(build / "precision")], text=True))
    assert precision == 8
    examples = [
        (
            "ultrasound",
            "lipidcoated_simple",
            ["-freq", "2.9e6", "-amp", "130e3", "-tend", "2e-6"],
            "plot_result_lipidcoated_simple.py",
        ),
        ("rayleighcollapse", "emissions", ["-tend", "0.2"], "plot_result_emissions.py"),
        (
            "ultrasound",
            "lipidcoated_emissions",
            ["-freq", "1e6", "-amp", "600e3", "-tend", "11e-6"],
            "plot_result_lipidcoated_emissions.py",
        ),
        (
            "ultrasound",
            "sonolum_emissions",
            ["-freq", "23.5e3", "-amp", "145e3", "-tend", "40e-6"],
            "plot_result_sonolum.py",
        ),
        (
            "ultrasound",
            "powerlaw",
            ["-freq", "636619", "-amp", "25331", "-tend", "25e-6"],
            "plot_result_powerlaw.py",
        ),
        (
            "ultrasound",
            "kelvinvoigt",
            ["-freq", "1e6", "-amp", "3e6", "-tend", "6e-6"],
            "plot_result_kelvinvoigt.py",
        ),
        ("ultrasound", "zener", ["-freq", "1e6", "-amp", "1e6", "-tend", "5e-6"], "plot_result_zener.py"),
        (
            "ultrasound",
            "oldroydb",
            ["-freq", "3e6", "-amp", "400e3", "-tend", "3e-6"],
            "plot_result_oldroydb.py",
        ),
    ]
    comparisons = {}
    for directory, case, args, plot in examples:
        source = upstream / "examples" / directory
        target = out / case
        target.mkdir()
        shutil.copy(source / (case + ".apecss"), target / "config.apecss")
        for reference in source.glob("reference-data*"):
            if reference.is_dir():
                shutil.copytree(reference, target / reference.name)
        executable = build / ("rayleigh_apecss" if directory == "rayleighcollapse" else "ultrasound_apecss")
        code = run([executable, "-options", "config.apecss", *args], target, "solver")
        if code != 0:
            continue
        run([sys.executable, source / plot], target, "official_plot", timeout=120)
        if case == "emissions":
            for ref in (target / "reference-data_emissions").glob("APECSS_*.txt"):
                observed = target / ref.name.removeprefix("APECSS_")
                a, b = np.loadtxt(ref), np.loadtxt(observed)
                comparisons[observed.name] = {
                    "shape": list(b.shape),
                    "reference_shape": list(a.shape),
                    "allclose_upstream_tolerance": bool(
                        a.shape == b.shape and np.allclose(a, b, rtol=1e-6, atol=3e-16)
                    ),
                    "relative_L2": float(np.linalg.norm(a - b) / np.linalg.norm(a))
                    if a.shape == b.shape
                    else None,
                }
    info = {
        "upstream_commit": commit,
        "release": "v1.8",
        "precision_bytes": precision,
        "compiler": subprocess.check_output(["gcc", "--version"], text=True),
        "build_flags": "Release: -O3 -Wall -Werror -march=native; no fast-math; double",
        "build_directory": str(build.relative_to(root)),
        "run": str(out.relative_to(root)),
        "commands": records,
        "reference_comparisons": comparisons,
        "upstream_patches": [],
        "note": "Official lipidcoated_simple has no tabulated reference curve in repository; a successful run alone is not experimental agreement.",
    }
    info["status"] = (
        "PASS"
        if records
        and all(r["returncode"] == 0 for r in records)
        and len(comparisons) == 4
        and all(r["allclose_upstream_tolerance"] for r in comparisons.values())
        else "PARTIAL_OR_FAILED"
    )
    write_json(out / "VALIDATION.json", info)
    write_json(root / "results/stage2/apecss_validation/STATE.json", info)


if __name__ == "__main__":
    main()
