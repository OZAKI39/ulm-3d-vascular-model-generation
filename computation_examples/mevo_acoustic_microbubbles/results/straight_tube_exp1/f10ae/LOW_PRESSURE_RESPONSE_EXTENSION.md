# Low-pressure response extension

Final opt-in extension admitted: **True**. New lower limit 0.001 Pa at 2.5 MHz, original shell x1 and R=1.78 um only. All original >=0.1 Pa intervals and CPU queries are retained exactly. The 0.1 Pa join differs from the refined direct response by about 4.84e-6 relative. No clipping, frequency interpolation, extrapolation, or near-zero force cutoff.

Three prior domain-failed actions complete in new explicitly named joined-table trials. Earlier historical runs are untouched. A separate frozen-action followup of the F10AE final_T4_1 failure reaches 0.966 ms rather than 0.237 ms under the admitted extension, but still fails the response-domain guard; no endpoint is assigned and the original final test remains failed. See LOW_DOMAIN_FOLLOWUP.json.

The initial admission summary was conservatively false because a NumPy Boolean was tested by identity. The extra `low_admit.py` phase both fixes aggregation and adds explicit sub-0.1 Pa incident-pressure dynamics. Search/validation had consistently retained the historical production domain; no comparison was switched halfway. `LOW_PRESSURE_ADMISSION.json` is the authoritative final opt-in gate; `LOW_PRESSURE_FINAL.json` remains as earlier evidence.

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
