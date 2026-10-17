# Current project handoff — 2026-10-17

The current simulation delivery is still **F10S**. Its code and scientific bytes
are identical to the 2026-10-15 public sync; start with
[F10S active code](F10S_ACTIVE_CODE_MAP.md) and
[the simulation handoff](CURRENT_PROJECT_STATE_20261015.md).

The latest related development is **F10U-Vision**, built on F10T-Vision-2 and the
human-reviewed annotation subsets. See
[the F10U report](../../../external_connect/vision_tracking_research/collective_motion_v1/FINAL_RECOMMENDATION.md),
[the annotation workbench guide](../../../external_connect/vision_tracking_research/smart_online_tracking_v2/ANNOTATION_READY.md),
and [the current Chinese progress copy](MICROBUBBLE_AGGREGATION_PROGRESS_20261017.md).

The public branch contains the source, reports, numerical evidence, saved labels,
frozen held-out data, final figures and small demonstration videos. Large raw
recordings, model weights and environments remain external. Exact scope, hashes,
restoration instructions and limitations are in
[the synchronization record](../sync/20261017/README.md).

Keep the distinction between visible image-plane distributions, geometric graph
connectivity and physically calibrated microbubble aggregation. F10U requires no
persistent bubble IDs and does not operate PZT/AFG hardware. Its difficult-case
accuracy and uncertainty checks still fail; live camera validation was unavailable.
No new scientific experiment was run during this synchronization.
