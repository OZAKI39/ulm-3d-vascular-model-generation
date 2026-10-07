"""Select the latest explicit campaign, retaining superseded evidence unchanged."""

from scripts.stage5a_r3.common import BASE, read


def campaign():
    v2 = BASE / "final/CAMPAIGN_R3_V2/RESULT.json"
    if v2.exists():
        return read(v2)
    closed = BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json"
    if closed.exists():
        value = read(closed)
        value["complete"] = False
        value["scope"] = "V1 safely closed after measured layout regression; overall validation continues"
        return value
    return read(BASE / "final/CAMPAIGN_R3_V1/RESULT.json")


def protocol_path():
    v2 = BASE / "final/VALIDATION_PROTOCOL_R3_V2.json"
    return v2 if v2.exists() else BASE / "final/VALIDATION_PROTOCOL_R3_V1.json"
