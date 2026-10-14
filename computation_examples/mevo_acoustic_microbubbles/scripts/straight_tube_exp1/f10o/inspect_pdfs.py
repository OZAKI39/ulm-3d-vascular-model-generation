"""Open every final PDF, inspect page text bounds and produce review images."""

import fitz
from .common import BASE, write, sha


def main():
    folder = BASE / "figures"
    out = folder / "pdf_previews"
    out.mkdir(exist_ok=True)
    records = []
    for p in sorted(folder.glob("*.pdf")):
        with fitz.open(p) as document:
            assert len(document) == 1
            page = document[0]
            rect = page.rect
            outside = []
            for block in page.get_text("dict")["blocks"]:
                for line in block.get("lines", []):
                    for span in line["spans"]:
                        if not rect.contains(fitz.Rect(span["bbox"])):
                            outside.append(span["text"])
            image = out / (p.stem + ".png")
            page.get_pixmap(matrix=fitz.Matrix(1.4, 1.4)).save(image)
            records.append(
                dict(
                    figure=p.stem,
                    pdf_sha256=sha(p),
                    preview_sha256=sha(image),
                    page_count=1,
                    text_outside_page=outside,
                    english_status_labels_present="EXPLORATORY" in page.get_text(),
                )
            )
    assert len(records) == 4 and all(
        not r["text_outside_page"] and r["english_status_labels_present"] for r in records
    )
    write(
        folder / "PDF_CHECKS.json",
        dict(
            status="PASS", figures=records, scope="Rendering and text-bounds check; manual review is separate"
        ),
    )


if __name__ == "__main__":
    main()
