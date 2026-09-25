# CF2003 corrigendum audit

Classification: **DOES_NOT_AFFECT_CLOSURE**. `CORRIGENDUM_CHECKED=YES`; `CORRIGENDUM_IMPACT=NONE`.

The complete one-page notice (QJMAM 65(4), 581, 2012; DOI [10.1093/qjmam/hbs012](https://doi.org/10.1093/qjmam/hbs012)) was obtained as the publisher's original page-preview GIF, linked directly from its [official article page](https://academic.oup.com/qjmam/article/65/4/581/1857053). This is primary publisher page imagery, not OCR. The page was inspected visually from its heading through the reference and DOI footer.

Material: `materials/corrigendum_preview.gif`, 803 × 1046 pixels, SHA256 `69c6cbb350127dd8b334998b04276641cf376b550e25c34ba2514078932265b8`. The image contains the entire notice. The direct PDF endpoint returned HTTP 403; no PDF was obtained, and no PDF hash is claimed. The source HTML, its original image URL and download receipt are retained in `materials/corrigendum_minimal.html` and `materials/third_retrieval.json`.

The notice identifies errors in the opening equations of subsection 2.1 and states: “The rest of the paper is unaffected.” It replaces the seven lines defining the perturbed translation velocity/pressure representation. Its changes concern nondimensionalization and component expressions of that field; the notice also says the cited earlier solution and calculations have no misprints. No replacement free-sphere coefficient is supplied.

| Requested closure item | Correction impact | Basis |
|---|---|---|
| Eq. 6.1 | DOES_NOT_AFFECT_CLOSURE | Outside the corrected §2.1 lines; explicit scope statement |
| Eq. 6.2 | DOES_NOT_AFFECT_CLOSURE | Same |
| Eq. 6.7 | DOES_NOT_AFFECT_CLOSURE | Same |
| Eq. 6.8 | DOES_NOT_AFFECT_CLOSURE | Same |
| Eq. 6.9 | DOES_NOT_AFFECT_CLOSURE | Same |
| Table 17 | DOES_NOT_AFFECT_CLOSURE | Same; no replacement coefficients |
| Ux-star normalization | DOES_NOT_AFFECT_CLOSURE | §6 normalization is outside correction scope |
| Omega-y-star normalization | DOES_NOT_AFFECT_CLOSURE | Same |
| Free-sphere shear results | DOES_NOT_AFFECT_CLOSURE | Explicit statement that the rest is unaffected |

This passes the correction-impact gate. The main paper was subsequently supplied by the user as a published library PDF and visually reviewed separately. No correction is applied to closure coefficients because the notice does not require one. The independent Table 17 extraction results and frozen reference are separate artifacts in this directory.
