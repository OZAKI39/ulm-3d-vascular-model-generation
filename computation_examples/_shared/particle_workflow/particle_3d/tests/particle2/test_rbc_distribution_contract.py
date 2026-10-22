import json
from pathlib import Path
import numpy as np


def test_contract_preserves_sources_and_model_roles(contract,p2_repo):
    c=contract
    assert c["contract_name"]=="C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0"
    assert c["species"]=="Mus musculus" and c["population"]=="mature RBC"
    assert c["exact_substrain"]=="MIXED_SOURCE_WITHIN_C57BL6_FAMILY"
    d=c["diameter"];v=c["volume"]
    assert [d[k] for k in ["mean_um","sd_um","source_mode_um","source_median_um","source_iqr_um","source_n"]]==[6.79,.93,6.67,6.67,1.33,1156720]
    assert [v[k] for k in ["mean_fL","cv","sd_fL","between_mouse_mcv_sd_fL"]]==[47.9,.148,7.0892,2.6]
    assert not v["between_mouse_sd_is_single_cell_sd"]
    assert v["cv_role"]=="MODEL_CHOICE_INFORMED_BY_C57BL6_RDW_LITERATURE"
    assert c["D_V_dependence"]=="INDEPENDENT_IN_V0_DUE_TO_MISSING_JOINT_DATA"
    assert c["D_V_dependence_role"]=="PROVISIONAL_MODEL_ASSUMPTION"
    assert "NOT_EXPERIMENTAL" in d["guard_role"] and "NOT_EXPERIMENTAL" in v["guard_role"]
    literature=json.loads((p2_repo/"particle_3d/reports/particle2/literature/sources.json").read_text())
    byid={s["id"]:s for s in literature["sources"]}
    moss=byid["MOSS_2025_PREPRINT"]
    assert moss["peer_reviewed"] is False and moss["publication_status"]=="PREPRINT"
    assert moss["doi"]=="10.1101/2025.03.05.639962" and moss["table3"]["n"]==1156720
    assert moss["table4"]["mcv_mean_fL"]==47.9 and moss["table4"]["rdw_mean_percent"]==17.3
    assert byid["RIVERA_2013"]["values"]=={"female_mcv_fL":47.8,"male_mcv_fL":48.4}
    assert byid["DE_FRANCESCHI_2005"]["values"]=={"mcv_mean_fL":49.3,"mcv_spread_fL":.9,"rdw_mean_percent":13.2,"rdw_spread_percent":.9}
    assert all(byid[k]["peer_reviewed"] for k in ["RIVERA_2013","DE_FRANCESCHI_2005"])


def test_guard_arithmetic_and_quaternion_convention_are_frozen(contract):
    for cfg,mean,sd,guard in [(contract["diameter"],"mean_um","sd_um","guard_um"),(contract["volume"],"mean_fL","sd_fL","guard_fL")]:
        np.testing.assert_allclose(cfg[guard],[cfg[mean]-3*cfg[sd],cfg[mean]+3*cfg[sd]],rtol=4*np.finfo(float).eps)
    o=contract["orientation"]
    assert o["quaternion_order"]=="wxyz" and o["rotation_map"]=="body_to_world"
    assert o["body_short_axis"]==[0,0,1] and o["world_short_axis"]=="p=R(q)@e3"
    assert not contract["production_particle_timestep_frozen"]
    assert not contract["production_orientation_distribution_frozen"]
