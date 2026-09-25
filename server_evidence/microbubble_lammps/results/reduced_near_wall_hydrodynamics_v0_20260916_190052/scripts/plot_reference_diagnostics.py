"""Small static input-reference plots; no kernel or simulation output."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[1]
data = np.genfromtxt(root / "raw/resistance_reference.csv", delimiter=",", names=True)
mask = data["epsilon"] <= .2
fig, ax = plt.subplots(2, 2, figsize=(9, 6), constrained_layout=True)
for axis, field, label in zip(ax.flat,
    ["normal_TT_over_6pi_mu_a", "tangential_TT_over_6pi_mu_a", "parallel_RR_over_8pi_mu_a3", "TR_x_Omega_y_over_6pi_mu_a2"],
    ["Normal TT / (6 pi mu a)", "Tangential TT / (6 pi mu a)", "Parallel RR / (8 pi mu a^3)", "TR x,Omega_y / (6 pi mu a^2)"]):
    axis.semilogx(data["epsilon"][mask], data[field][mask])
    axis.set(xlabel="epsilon = surface gap / radius", ylabel=label)
    axis.grid(alpha=.3)
fig.suptitle("Frozen RMBW TOTAL resistance input | RR/TR reference limitations retained\nNo new kernel qualification")
fig.savefig(root / "visualization/resistance_vs_epsilon.png", dpi=150)
plt.close(fig)

rows = list(csv.DictReader((root / "raw/gcb_transcription_differences.csv").open()))
grid = np.geomspace(.001,.05,150)
den = .6376-.2*np.log(grid)
with (root / "raw/gcb_printed_expression_diagnostic.csv").open("w",newline="") as out:
    w=csv.writer(out); w.writerow(["epsilon","FU_printed_candidate","FOmega_printed_candidate","status"])
    for e,u,r in zip(grid,.7431/den,.8436/den):w.writerow([e,u,r,"UNQUALIFIED_EXPRESSION_DIAGNOSTIC"])
for field, numerator, filename in [("FU",.7431,"gcb_translation_vs_epsilon.png"),("FOmega",.8436,"gcb_rotation_vs_epsilon.png")]:
    fig,ax=plt.subplots(figsize=(8,4.5),constrained_layout=True)
    ax.semilogx(grid,numerator/den,label="Printed asymptotic expression (unqualified)")
    for table,marker in [("Table2","o"),("Table3","x")]:
        pts=[r for r in rows if r["table"]==table and float(r["epsilon"])>=.001]
        ax.scatter([float(r["epsilon"]) for r in pts],[float(r["transcribed_"+field]) for r in pts],marker=marker,label=table+" public transcription (unverified)")
    ax.set(xlabel="epsilon = surface gap / radius",ylabel=field,
           title="BLOCKED: reference discrepancy | no certified epsilon domain")
    ax.grid(alpha=.3);ax.legend(fontsize=8)
    fig.savefig(root / "visualization" / filename,dpi=150)
    plt.close(fig)
print("3 reference-only diagnostic figures written; no shear RHS figure or kernel output.")
