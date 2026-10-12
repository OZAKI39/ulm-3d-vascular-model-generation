"""CPU-only deterministic split controller. It never accepts rejected proposals."""


def isolate(records, parent, run, attempts=None):
    """Split a known failed batch; `run` always receives original immutable records.

    Failed singleton leaves are genuine independent failures. A successful
    child may contain several members. Execution/input errors remain unresolved.
    Original parent evidence is retained by reference and never re-executed.
    """
    attempts = [] if attempts is None else attempts
    ids = [r["sample_id"] for r in records]
    if len(set(ids)) != len(ids) or not ids:
        raise ValueError("Nonempty unique identities required")
    if len(records) == 1:
        raise ValueError("Known singleton failure does not need batch isolation")
    if any(r["events"] != records[0]["events"] for r in records):
        raise ValueError("Synchronous control family required")

    def recurse(group, provenance):
        try:
            result = run(group, provenance)
        except (OSError, ValueError) as exc:
            result = dict(status="UNRESOLVED", reason=str(exc), attempt_id=provenance + "_unresolved")
        result = dict(result, parent_attempt=provenance, sample_ids=[r["sample_id"] for r in group])
        attempts.append(result)
        if result["status"] == "COMPLETE":
            state = "RECOVERED_BY_ISOLATION"
        elif result["status"] == "FAILED" and len(group) == 1:
            state = "FAILED_INDEPENDENT"
        elif result["status"] == "FAILED":
            mid = len(group) // 2
            return recurse(group[:mid], result["attempt_id"]) + recurse(group[mid:], result["attempt_id"])
        else:
            state = "ABORTED_BATCH_UNRESOLVED"
        return [
            dict(sample_id=r["sample_id"], classification=state, attempt_id=result["attempt_id"])
            for r in group
        ]

    mid = len(records) // 2
    leaves = recurse(records[:mid], parent) + recurse(records[mid:], parent)
    return dict(parent_attempt=parent, sample_ids=ids, attempts=attempts, leaves=leaves)
