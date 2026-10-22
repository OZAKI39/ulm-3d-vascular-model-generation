# Stage5A stopped by user — 2026-10-06

The user explicitly requested 停止运行. The sole new campaign service `mevo_stage5a_motion_20261006` received SIGINT and exited at 20:33 UTC. PID 911451 no longer exists. No automatic resume is authorized. No unrelated process was stopped.

Branch: `dev/stage5a-motion-constrained-20261006`; committed implementation `ce1808fd1`; later formatting, diagnostics, analysis scripts and tests remain uncommitted. No final GitHub delivery has been pushed.

B2 nominal and dt/2 both failed after 2.109375 µs (65 accepted steps each); full failure evidence and accepted trajectories are preserved. B3 nominal was interrupted after 19.296875 µs, with 263 accepted and 264 rejected intervals. Its finally block saved LAST_ACCEPTED_STATE.npz and flushed both ledgers. Its full in-memory trajectory HDF5 and SUMMARY were NOT finalized. The last coordinate snapshot is not an exact warm-cache restart checkpoint; do not fabricate the missing intermediate coordinates or call a rerun the original trajectory. B3 dt/2 was not launched. The joint 1 ms gate cannot pass because B2 already failed; 5 ms and longer runs remain closed.

The new backend is not promoted. Earlier tests recorded 97 passes; later evidence-sampling test changes and the retained real-QP xfail have not received a final combined test run. B2 sparse failed-QP reconstruction and baseline benchmark/transport audit are retained. Final figures, development report, final artifact verification and final commit/push remain unfinished.

Server evidence was copied locally without modifying the original Stage4L data. Authoritative stop record: `results/stage5/motion_solver/STOPPED_BY_USER_20261006T203352Z.json`.

On a future explicit resume, first reconcile the interrupted B3 artifact limitations and the user's desired scope. Current deliver/compare scripts still expect finalized nominal and dt/2 summaries and HDF5; do not run them blindly. Preserve this interruption and all failures under their original run identities. Analysis scripts were edited after the executed scientific source; numerical module changes are formatting-only with identical AST to ce1808fd1. Some final analysis edits still need formatting/lint checks. Do not overwrite live scientific provenance with later source hashes.
