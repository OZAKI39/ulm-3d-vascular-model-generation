"""Make accidental access to either historical checkout a hard test failure."""

import os
import sys


def forbid_historical_files(event, args):
    if event != "open" or not args:
        return
    path = args[0]
    if not isinstance(path, (str, bytes, os.PathLike)):
        return
    path = os.fsdecode(path)
    if (
        "/home/lzy/projects/temp_storage/" in path
        or "/workspace/y_clot/project/" in path
    ):
        raise AssertionError(f"Test attempted to access historical file: {path}")


sys.addaudithook(forbid_historical_files)
