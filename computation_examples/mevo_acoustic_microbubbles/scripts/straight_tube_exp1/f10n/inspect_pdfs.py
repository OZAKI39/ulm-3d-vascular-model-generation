"""Render actual PDF pages for subsequent visual inspection; do not assert visual PASS."""

import fitz
from .common import BASE, write, sha


def main():
    folder = BASE / "figures"
    out = folder / "pdf_previews"
    out.mkdir(exist_ok=True)
    rows = []
    for path in sorted(folder.glob("Fig_F10N_*.pdf")):
        with fitz.open(path) as pdf:
            assert len(pdf) == 1
            page = pdf[0]
            preview = out / (path.stem + ".png")
            page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False).save(preview)
            rows.append(
                dict(
                    file=path.name,
                    sha256=sha(path),
                    pages=1,
                    text_characters=len(page.get_text()),
                    font_count=len(page.get_fonts()),
                    page_size_pt=list(page.rect),
                    preview=str(preview.relative_to(folder)),
                    preview_sha256=sha(preview),
                )
            )
    assert len(rows) == 4
    write(
        folder / "PDF_INSPECTION.json",
        dict(structural_status="PASS", visual_review="Separate actual opening required", files=rows),
    )


if __name__ == "__main__":
    main()
