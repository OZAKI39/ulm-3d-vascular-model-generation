"""Consistent, dependency-free console formatting for long physics runs."""

from __future__ import annotations

from collections.abc import Iterable
import shutil
import sys
import textwrap
from typing import TextIO


PROGRESS_BAR_FORMAT = (
    "{desc:<22} {percentage:3.0f}% |{bar}| "
    "{n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]"
)


def print_banner(
    title: str,
    *,
    subtitle: str | None = None,
    file: TextIO | None = None,
) -> None:
    """Print a prominent run-level banner without relying on ANSI colors."""

    stream = _stream(file)
    width = _console_width()
    rule = "=" * width
    print(rule, file=stream)
    print(f" {title}", file=stream)
    if subtitle:
        print(f" {subtitle}", file=stream)
    print(rule, file=stream)


def print_stage(
    current: int,
    total: int,
    title: str,
    *,
    detail: str | None = None,
    file: TextIO | None = None,
) -> None:
    """Print one numbered top-level pipeline stage."""

    stream = _stream(file)
    print(file=stream)
    print(f"[{int(current):02d}/{int(total):02d}] {title}", file=stream)
    if detail:
        print(f"        {detail}", file=stream)


def print_section(title: str, *, file: TextIO | None = None) -> None:
    """Print a compact divider for a nested numerical report."""

    stream = _stream(file)
    width = _console_width()
    prefix = f"-- {title} "
    print(file=stream)
    print(prefix + "-" * max(3, width - len(prefix)), file=stream)


def print_key_values(
    items: Iterable[tuple[str, object]],
    *,
    indent: int = 2,
    file: TextIO | None = None,
) -> None:
    """Print aligned labels and values, preserving readable multiline values."""

    stream = _stream(file)
    rows = [(str(label), str(value)) for label, value in items]
    if not rows:
        return
    label_width = min(max(len(label) for label, _ in rows), 44)
    padding = " " * max(0, int(indent))
    continuation = " " * (len(padding) + label_width + 3)
    for label, value in rows:
        value_lines = value.splitlines() or [""]
        print(f"{padding}{label:<{label_width}} : {value_lines[0]}", file=stream)
        for line in value_lines[1:]:
            print(f"{continuation}{line}", file=stream)


def print_command(command: object, *, file: TextIO | None = None) -> None:
    """Print one copy-ready command without a key/value separator."""

    stream = _stream(file)
    print(f"  {command}", file=stream)


def print_warning(message: str, *, file: TextIO | None = None) -> None:
    """Print a wrapped warning that remains legible in redirected log files."""

    _print_status("WARNING", message, file=file)


def print_error(message: str, *, file: TextIO | None = None) -> None:
    """Print a wrapped fatal-error message."""

    _print_status("ERROR", message, file=sys.stderr if file is None else file)


def print_success(message: str, *, file: TextIO | None = None) -> None:
    """Print a concise successful-completion marker."""

    _print_status("OK", message, file=file)


def _print_status(level: str, message: str, *, file: TextIO | None) -> None:
    stream = _stream(file)
    prefix = f"[{level}] "
    subsequent = " " * len(prefix)
    width = max(20, _console_width() - len(prefix))
    lines = textwrap.wrap(
        " ".join(str(message).split()),
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    ) or [""]
    print(prefix + lines[0], file=stream)
    for line in lines[1:]:
        print(subsequent + line, file=stream)


def _console_width() -> int:
    columns = int(shutil.get_terminal_size(fallback=(96, 24)).columns)
    return max(72, min(columns, 112))


def _stream(file: TextIO | None) -> TextIO:
    return sys.stdout if file is None else file
