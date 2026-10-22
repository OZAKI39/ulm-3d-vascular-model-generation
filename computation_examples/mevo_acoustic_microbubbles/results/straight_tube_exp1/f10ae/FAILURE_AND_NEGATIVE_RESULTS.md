# Failures and negative results

Physical attempts: 230; complete: 227; incomplete: 3. Radial 0.0001 Pa finest solve timed out. All original raw attempted states, forces, exceptions and accepted prefixes remain. Failed integrations have `final=null` and no target loss.

| Target | Frozen source | Discovery aspect change | COM (nm) | Independent confirmations | Final test raw hits | Confirmed |
| --- | --- | --- | --- | --- | --- | --- |
| T1 | extension_028 | 2.9974% | 2.1325 | 2/3 | 2/3 | False |
| T2 | ax_14 | 0.4096% | 1.2628 | 2/3 | 3/3 | False |
| T3 | baseline_OFF | 0.0000% | 0.0000 | 0/3 | 0/3 | False |
| T4 | recovery_05 | 3.0208% | 0.3866 | 0/3 | 0/3 | False |

Numerical comparison rejects: []

Physical failures: [
  {
    "id": "final_T4_1",
    "status": "OUT_OF_MODEL_DOMAIN",
    "reasons": [
      "EXCEPTION: ValueError('Original acoustic/APECSS fail-stop bits=2')"
    ],
    "completed_s": 0.000237
  },
  {
    "id": "joined_followup_final_T4_1",
    "status": "OUT_OF_MODEL_DOMAIN",
    "reasons": [
      "EXCEPTION: ValueError('Original acoustic/APECSS fail-stop bits=2')"
    ],
    "completed_s": 0.000966
  },
  {
    "id": "recovery_20",
    "status": "NUMERICAL_OR_BUDGET_FAILURE",
    "reasons": [
      "CONTACT_REACTION_BALANCE_PASS"
    ],
    "completed_s": 0.002052
  }
]

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
