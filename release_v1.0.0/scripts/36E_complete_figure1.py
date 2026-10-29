from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path("results")
OUT = ROOT / "manuscript_final"
FIG = OUT / "figures"
SRC = OUT / "source_data"

FIG.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

discovery = pd.read_csv(
    ROOT / "pseudobulk" /
    "GSE252921_endothelial_donor_metadata.csv"
)

validation = pd.read_csv(
    ROOT / "validation_GSE327289" /
    "GSE327289_endothelial_donor_metadata.csv"
)

audit = pd.read_csv(
    ROOT / "AI_DNN" /
    "32C_donor_overlap_audit.csv"
)

def summarize(df, donor_col, disease_col, dataset):
    required = [donor_col, disease_col, "n_cells"]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"{dataset}: missing {missing}")

    if df[donor_col].isna().any():
        raise ValueError(f"{dataset}: missing donor IDs")

    if df[donor_col].duplicated().any():
        raise ValueError(f"{dataset}: duplicated donor IDs")

    labels = (
        df[disease_col].astype(str)
        .str.strip().str.upper()
    )

    print(dataset, "diagnosis labels:",
          labels.value_counts().to_dict())

    # Explicitly check diagnosis labels instead of
    # silently assigning unknown categories.
    def category(value):
        if value in {"AD", "ALZHEIMER", "ALZHEIMER'S DISEASE"}:
            return "AD"
        if value in {"CT", "CTRL", "CONTROL", "HC", "HEALTHY"}:
            return "Control"
        raise ValueError(
            f"Unrecognized diagnosis label: {value}"
        )

    groups = labels.map(category)

    cells = pd.to_numeric(
        df["n_cells"], errors="raise"
    )

    if cells.isna().any() or (cells < 0).any():
        raise ValueError("Invalid cell counts")

    return {
        "dataset": dataset,
        "donors": len(df),
        "AD_donors": int((groups == "AD").sum()),
        "control_donors": int((groups == "Control").sum()),
        "endothelial_cells": int(cells.sum())
    }

d = summarize(
    discovery, "donor_ID", "disease", "GSE252921"
)

v = summarize(
    validation, "Donor", "Diagnosis", "GSE327289"
)

summary = pd.DataFrame([d, v])

print("\nCOHORT SUMMARY")
print(summary.to_string(index=False))

assert d["donors"] == 20
assert v["donors"] == 23
assert d["endothelial_cells"] == 50674
assert v["endothelial_cells"] == 6250
assert d["AD_donors"] == 10
assert d["control_donors"] == 10
assert v["AD_donors"] == 16
assert v["control_donors"] == 7

summary.to_csv(
    SRC / "Figure1_cohort_summary.csv",
    index=False
)

discovery.to_csv(
    SRC / "Figure1_discovery_donor_metadata.csv",
    index=False
)

validation.to_csv(
    SRC / "Figure1_validation_donor_metadata.csv",
    index=False
)

audit.to_csv(
    SRC / "Figure1_donor_overlap_audit.csv",
    index=False
)

fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1.05, 1],
    hspace=0.40,
    wspace=0.35
)

ax_a = fig.add_subplot(gs[0, :])
ax_b = fig.add_subplot(gs[1, 0])
ax_c = fig.add_subplot(gs[1, 1])

# Panel a: parallel datasets
ax_a.set_xlim(0, 10)
ax_a.set_ylim(0, 6)
ax_a.axis("off")

def box(ax, x, y, w, h, title, subtitle, color):
    rect = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.12",
        linewidth=1,
        edgecolor=color,
        facecolor="white"
    )
    ax.add_patch(rect)
    ax.text(
        x+w/2, y+h*0.68,
        title, ha="center",
        fontsize=13, fontweight="bold"
    )
    ax.text(
        x+w/2, y+h*0.33,
        subtitle, ha="center",
        fontsize=10
    )

box(
    ax_a, 0.7, 3.3, 3.6, 1.5,
    "GSE252921 — Discovery",
    "20 donors | 50,674 endothelial cells",
    "#277DA1"
)

box(
    ax_a, 5.7, 3.3, 3.6, 1.5,
    "GSE327289 — Validation",
    "23 donors | 6,250 endothelial cells",
    "#F8961E"
)

for x in [2.5, 7.5]:
    ax_a.annotate(
        "", xy=(x, 2.2), xytext=(x, 3.2),
        arrowprops=dict(arrowstyle="->", lw=1.5)
    )

box(
    ax_a, 2.1, 0.45, 5.8, 1.55,
    "Harmonized cross-dataset analysis",
    "10,901 shared genes | 1,168 pathways",
    "#7A5195"
)

ax_a.set_title(
    "a  Study design and analysis workflow",
    loc="left",
    fontweight="bold"
)

# Panel b: donor composition
labels = summary["dataset"].tolist()
x = np.arange(2)

control = summary["control_donors"].to_numpy()
ad = summary["AD_donors"].to_numpy()

ax_b.bar(
    x, control, color="#277DA1",
    label="Control"
)
ax_b.bar(
    x, ad, bottom=control,
    color="#F8961E", label="AD"
)

for i in range(2):
    ax_b.text(
        x[i],
        control[i] + ad[i] + 0.6,
        f"n={control[i]+ad[i]}\n"
        f"{summary.loc[i,'endothelial_cells']:,} cells",
        ha="center",
        fontsize=9
    )

ax_b.set_xticks(x)
ax_b.set_xticklabels(labels)
ax_b.set_ylim(0, 29)
ax_b.set_ylabel("Number of donors")
ax_b.legend(frameon=False)
ax_b.set_title(
    "b  Cohort composition",
    loc="left",
    fontweight="bold"
)

# Panel c: statistical specifications
ax_c.axis("off")
ax_c.set_title(
    "c  Donor-aware statistical analysis",
    loc="left",
    fontweight="bold"
)

text = (
    "Discovery model\n"
    "Disease + age + sex + batch\n\n"
    "Validation model\n"
    "Diagnosis + age + sex\n\n"
    "Sensitivity analyses\n"
    "Low-cell donor exclusion\n"
    "Endothelial composition adjustment\n"
    "Capillary endothelial subset\n\n"
    "Donor independence across datasets:\n"
    "Not established"
)

ax_c.text(
    0.04, 0.93,
    text,
    transform=ax_c.transAxes,
    va="top",
    fontsize=11,
    linespacing=1.5
)

for ax in [ax_b]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.savefig(
    FIG / "Figure1_Complete.pdf",
    bbox_inches="tight"
)

fig.savefig(
    FIG / "Figure1_Complete.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print("\nFIGURE 1 COMPLETE")
print("Saved Figure1_Complete.pdf")
print("Saved Figure1_Complete.png")
print("Source Data saved.")
