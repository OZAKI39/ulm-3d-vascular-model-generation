"""Select the latest explicit campaign, retaining superseded evidence unchanged."""

from scripts.stage5a_r3.common import BASE, read


def campaign():
    v2 = BASE / "final/CAMPAIGN_R3_V2/RESULT.json"
    if v2.exists():
        value = read(v2)
        supplement = BASE / "final/CAMPAIGN_R2_PROTECTED_REPEATS/RESULT.json"
        if supplement.exists():
            extra = read(supplement)
            value["ledger"] += extra["ledger"]
            value["spent_s"]["r2"] = extra["total_spent_r2_s"]
            value["complete"] = value["complete"] and extra["complete"]
            value["validation_supplement"] = str(supplement)
        elif (BASE / "final/PROTECTED_REPEAT_PLAN.json").exists():
            value["complete"] = False
            value["scope"] = (
                "Two presealed protected R2 B3 timing repeats pending within the original group budget"
            )
        return value
    closed = BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json"
    if closed.exists():
        value = read(closed)
        value["complete"] = False
        value["scope"] = "V1 safely closed after measured layout regression; overall validation continues"
        return value
    return read(BASE / "final/CAMPAIGN_R3_V1/RESULT.json")


def protocol_path():
    version = campaign().get("version", "R3_V1")
    return BASE / f"final/VALIDATION_PROTOCOL_{version}.json"
