# Zhang 2011 N=2 reference

Source: [publisher PDF](https://cpb.iphy.ac.cn/en/article/pdf/preview/10.1088/1674-1056/20/11/114302.pdf), DOI 10.1088/1674-1056/20/11/114302.
Official article page lists no supplement; author/title/DOI code searches found no usable public implementation.
This is a minimal implementation of Eqs 3–9, not a new general bubble library.
The four second-order equations are solved simultaneously as a 4×4 acceleration system.
All quantities are SI. Pdot includes both partial derivatives of Eq9. Drive is −Pa sin(2πft).
χ=1 N/m is chosen from the paper's χ≤1 statement; the figure-specific χ and averaging window are not explicit.
No buckling/rupture or calibration to SonoVue is implied by this elastic-shell model.

Run `OPENBLAS_NUM_THREADS=1 python scripts/stage2/pair_validate.py` on the server.
The run of 2026-10-01 reproduces the initial attraction/repulsion directions of Fig4 and the sampled Fig5 peak near 2.9 µm.
It **fails** Fig3's claimed positive coefficient at 120 µm for R1=3 µm, R2=6 µm.
The sign remains negative after 64 cycles; reducing the integration step does not remove the discrepancy.
The run is retained as `FAILED_TREND_GATE`. Do not tune χ, change signs, or reinterpret this as full validation.
Only the first 0.2 ms of Fig4 was computed; the final stable pair is not reproduced/validated.
The source Fig4 axes are milliseconds; its prose printing “25 µm / 0.4 µm” for elapsed time is inconsistent.
No digitized quantitative curve comparison is claimed. Figure3 sign failure blocks Stage2H/J production runs.
