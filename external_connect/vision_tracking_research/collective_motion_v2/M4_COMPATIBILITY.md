# M4 compatibility

Original `M4-collective-1.0` is unchanged. New `configs/collective_state_v2.schema.json` uses `M4-collective-2.0` and VALID_GLOBAL_MOTION / NON_AFFINE_MOTION / NOT_OBSERVABLE. Model name, weight SHA and calibration version are required. The old-version adapter maps non-affine motion to null invalid fields with an explicit V2 reason; it never inserts unsupported enums into old messages.

Values/valid/reasons are retained for each field. This prototype's operational joint-risk policy is conservative: rejection by the six-parameter gate invalidates all six fields. The interface still supports partial fields, but the current model does not claim separately calibrated partial observability. Diagnostic `raw_estimate` is not the published valid observation.

Session, frame identity, actual PTS, fixed ROI, historical time span, interval reference, host age and no-new-frame TTL are retained. VersionedGate checks gaps; F10U ConsumerGate checks old/duplicate PTS, frame identity and stale age. CausalWindow clears on a discontinuity. No invalid observation is replaced by zero. All actions remain null; hardware_control_executed=false and training_ready=false.

Unit tests execute JSON Schema validation, rejection of new messages by the old schema, explicit old-schema mapping, gap/duplicate/TTL cases and a real local ZeroMQ MockConsumer. Causal recordings also pass through the schema and the old mock consumer. This is software observation compatibility, not hardware control admission.
