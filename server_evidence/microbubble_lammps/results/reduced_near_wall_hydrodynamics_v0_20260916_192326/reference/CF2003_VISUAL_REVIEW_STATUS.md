# Primary-page visual review — complete

Status: **VERIFIED**. The user supplied `/home/lzy/projects/560381.pdf`, a 30-page published paper with an Oxford download watermark for Universitaetsbibliothek Bern dated 16 September 2026. The title, authors, journal, year, page range, equations and tables were visually checked. This is source hierarchy B: a user-provided institutional/library copy of the published article. The original file is read-only; its stage copy has SHA256 `2d572961b10cc7b46ae2a666f8dfc80123ee303b118d964c510ae86a739dc3e5`.

The existing PyMuPDF 1.24.10 environment rendered the actual PDF pages. No package was installed. Saved page images and their payload hashes support repeat review.

| Item | Primary page | Status |
|---|---|---|
| Paper identity | p381 | VISUALLY_VERIFIED |
| Radius a, center height l, z=0 plane | p383 | VISUALLY_VERIFIED |
| Section 6 and free/force-free/torque-free condition | pp402–403 | VISUALLY_VERIFIED |
| Eqs. 6.1 and 6.2 | p404 | VISUALLY_VERIFIED |
| Ux/(kappa*l), Omega-y/(kappa/2), both tending to 1 | p404 | VISUALLY_VERIFIED |
| Eqs. 6.5, 6.6, 6.7, 6.8 | p404 | VISUALLY_VERIFIED |
| Eq. 6.9 | p407 | VISUALLY_VERIFIED |
| Table 17, both columns, all j=0..19 | p407 | VISUALLY_VERIFIED |
| epsilon=l/a-1; small-gap approximation and crossover 0.4 | pp407–408 | VISUALLY_VERIFIED |

TRANSCRIPTION_A was entered manually from the rendered p407 image and saved before the separate primary-PDF text extraction. TRANSCRIPTION_B was parsed from PyMuPDF page text. Only decimal middle dots and Unicode minus signs were normalized. The 40 decimal literal strings match exactly, including sign, significant digits and exponent. The requested verification CSV now contains both readings and final values. The preceding incomplete table is retained only in raw/CF2003_TABLE17_PENDING_BEFORE_PRIMARY_PDF.csv. The earlier Doczz text is additional unverified historical acquisition evidence, not the authoritative numeric source.

There are 20 coefficients per column (degree 19), indexed 0–19. Negative indices: u={1,3}, omega={1,3,4}; all others positive. Scientific notation starts at u9 and omega4, with omega5–7 returning to decimal notation. Both extraction methods were checked by the same assistant; this is independent extraction, not independent reviewers or independent physical theories.

Normalization: h=l-a, epsilon=h/a, d=l=a+h, FU=Ux/(kappa*l), FOMEGA=Omega-y/(kappa/2). Eq. 6.1 uses numerator f_k*c_r + f_r*c_k*a/(2*l), and Eq. 6.2 uses f_t*c_k + f_k*c_t*2*l/a; both use denominator f_t*c_r-f_r*c_t. Eqs. 6.7/6.8 are reciprocals of polynomials in natural log epsilon, not polynomials of FU/FOMEGA themselves.

With n=+z and g=+kappa*x, n cross g points +y; rotation follows the right-hand convention. The far-field normalization is verified from the definition, but the small-gap polynomial is never extrapolated to infinity. The project certifies only [0.001,0.2], inside the paper's stated small-gap crossover at 0.4 and the RMBW lower limit of 0.001.

The complete official 2012 correction page was visually checked separately. It does not affect the above closure; see CF2003_CORRIGENDUM_AUDIT.md.
