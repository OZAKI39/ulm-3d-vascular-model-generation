"""Small self-contained HTML acceptance report."""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Any

from ..reporting.acceptance import AcceptanceResult


def write_html_report(
    path: Path,
    *,
    title: str,
    summary: dict[str, Any],
    acceptance: AcceptanceResult | None,
    visualizations: list[Path],
) -> Path:
    status = acceptance.overall_status if acceptance else "INFO"
    checks = acceptance.checks if acceptance else []
    rows = "".join(
        f"<tr><td class='{check.status.lower()}'>{escape(check.status)}</td>"
        f"<td>{escape(check.name)}</td><td>{escape(check.message)}</td></tr>"
        for check in checks
    )
    summary_rows = "".join(
        f"<tr><th>{escape(str(key))}</th><td>{escape(str(value))}</td></tr>"
        for key, value in summary.items()
        if not isinstance(value, (dict, list))
    )
    images = "".join(
        f"<figure><img src='{escape(path.relative_to(path.parents[1]).as_posix())}' "
        f"alt='{escape(path.name)}'><figcaption>{escape(path.name)}</figcaption></figure>"
        for path in visualizations
        if path.is_file()
    )
    html = f"""<!doctype html>
<html lang='zh-CN'><head><meta charset='utf-8'><title>{escape(title)}</title>
<style>
body{{font-family:Arial,sans-serif;margin:2rem;color:#222}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}
th,td{{border:1px solid #ccc;padding:.45rem;text-align:left;vertical-align:top}}th{{background:#f5f5f5}}
.pass{{color:#087830;font-weight:bold}}.warning{{color:#a76500;font-weight:bold}}.fail{{color:#b00020;font-weight:bold}}
figure{{margin:1.5rem 0}}img{{max-width:100%;border:1px solid #ddd}}code{{background:#f4f4f4;padding:.1rem .25rem}}
</style></head><body><h1>{escape(title)}</h1>
<p><strong>Status: {escape(status)}</strong></p>
<p><strong>Direction convention:</strong> SWC <code>parent_id node → current node</code>. Arrows encode an inferred structural direction only; they are not measured blood-flow velocity or pressure.</p>
<h2>Summary</h2><table>{summary_rows}</table>
<h2>Acceptance checks</h2><table><tr><th>Status</th><th>Check</th><th>Evidence</th></tr>{rows}</table>
<h2>Visual evidence</h2>{images or '<p>No visualizations generated.</p>'}
</body></html>"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path
