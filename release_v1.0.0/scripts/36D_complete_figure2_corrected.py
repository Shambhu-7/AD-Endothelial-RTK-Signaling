from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("results")
OUT = ROOT / "manuscript_final" / "corrected_final"
FIG = OUT / "figures"
SRC = OUT / "source_data"

FIG.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

D_FILE = ROOT / "manuscript_final" / "reproducibility" / "dispersion_check" / "primary" / "GSE252921_primary_edgeR_results.csv"

V_FILE = ROOT / "validation_GSE327289" / \
    "GSE327289_edgeR_all_genes.csv"

REQUIRED = ["gene", "logFC", "PValue", "FDR"]

def load_de(path):
    df = pd.read_csv(path)

    missing = set(REQUIRED) - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing {missing}")

    if df["gene"].isna().any():
        raise ValueError(f"{path}: missing gene IDs")

    if df["gene"].duplicated().any():
        raise ValueError(f"{path}: duplicate gene IDs")

    for col in ["logFC", "PValue", "FDR"]:
        df[col] = pd.to_numeric(
            df[col], errors="raise"
        )

    if not np.isfinite(df["logFC"]).all():
        raise ValueError("Nonfinite logFC")

    if not (
        df["PValue"].between(0, 1).all()
        and df["FDR"].between(0, 1).all()
    ):
        raise ValueError("Invalid P values or FDR")

    return df.copy()

d = load_de(D_FILE)
v = load_de(V_FILE)

print("DISCOVERY GENES:", len(d))
print("VALIDATION GENES:", len(v))

assert len(d) == 18321, "Unexpected discovery gene count"
assert len(v) == 11029, "Unexpected validation gene count"

d_sig = d["FDR"] < 0.05
v_sig = v["FDR"] < 0.05

print("DISCOVERY FDR<0.05:", int(d_sig.sum()))
print("VALIDATION FDR<0.05:", int(v_sig.sum()))

assert d_sig.sum() == 5, "Unexpected discovery DEG count"
assert v_sig.sum() == 1, "Unexpected validation DEG count"

shared = d[
    ["gene", "logFC", "PValue", "FDR"]
].merge(
    v[["gene", "logFC", "PValue", "FDR"]],
    on="gene",
    suffixes=("_discovery", "_validation"),
    validate="one_to_one"
)

shared["same_direction"] = (
    np.sign(shared["logFC_discovery"]) ==
    np.sign(shared["logFC_validation"])
)

shared["significant_both"] = (
    (shared["FDR_discovery"] < 0.05)
    & (shared["FDR_validation"] < 0.05)
)

print("SHARED GENES:", len(shared))
print(
    "SIGNIFICANT IN BOTH:",
    int(shared["significant_both"].sum())
)

assert shared["significant_both"].sum() == 0, \
    "Unexpected cross-dataset DEG replication"

# Save frozen-data-derived source tables
d.to_csv(
    SRC / "Figure2A_discovery_DE.csv",
    index=False
)
v.to_csv(
    SRC / "Figure2B_validation_DE.csv",
    index=False
)
shared.to_csv(
    SRC / "Figure2C_shared_gene_effects.csv",
    index=False
)

# Four-panel figure
fig, axes = plt.subplots(
    2, 2,
    figsize=(12, 10),
    constrained_layout=True
)

ax_a, ax_b, ax_c, ax_d = axes.flatten()

def volcano(ax, df, title, panel):
    sig = df["FDR"] < 0.05
    y = -np.log10(
        np.maximum(
            df["PValue"].to_numpy(),
            np.finfo(float).tiny
        )
    )

    ax.scatter(
        df.loc[~sig, "logFC"],
        y[~sig],
        s=5,
        alpha=0.28,
        color="#A8AFB8",
        rasterized=True
    )

    ax.scatter(
        df.loc[sig, "logFC"],
        y[sig],
        s=22,
        alpha=0.9,
        color="#C44E52",
        rasterized=True
    )

    # Annotate only FDR-significant genes.
    significant = df.loc[sig].copy()

    for _, r in significant.iterrows():
        idx = r.name
        ax.annotate(
            str(r["gene"]),
            (r["logFC"], y[df.index.get_loc(idx)]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7
        )

    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("log₂ fold change (AD vs control)")
    ax.set_ylabel("−log₁₀ nominal P value")
    ax.set_title(
        f"{panel}  {title}\n"
        f"{sig.sum()} genes at FDR < 0.05",
        loc="left",
        fontweight="bold"
    )

volcano(
    ax_a, d,
    "Discovery: GSE252921", "a"
)

volcano(
    ax_b, v,
    "Validation: GSE327289", "b"
)

# Panel c: shared gene effect sizes
ax_c.scatter(
    shared["logFC_discovery"],
    shared["logFC_validation"],
    s=6,
    alpha=0.25,
    color="#277DA1",
    rasterized=True
)

ax_c.axhline(0, color="gray", lw=0.8)
ax_c.axvline(0, color="gray", lw=0.8)

ax_c.set_xlabel("Discovery log₂ fold change")
ax_c.set_ylabel("Validation log₂ fold change")
ax_c.set_title(
    "c  Cross-dataset gene-level effects",
    loc="left",
    fontweight="bold"
)

# Panel d: significant gene counts
labels = [
    "Discovery",
    "Validation",
    "Significant\nin both"
]

counts = [
    int(d_sig.sum()),
    int(v_sig.sum()),
    int(shared["significant_both"].sum())
]

bars = ax_d.bar(
    labels,
    counts,
    color=["#277DA1", "#F8961E", "#7A5195"]
)

for bar, count in zip(bars, counts):
    ax_d.text(
        bar.get_x() + bar.get_width()/2,
        bar.get_height() + 0.12,
        str(count),
        ha="center",
        fontsize=11
    )

ax_d.set_ylim(0, 6.5)
ax_d.set_ylabel("Number of FDR-significant genes")
ax_d.set_title(
    "d  Limited gene-level replication",
    loc="left",
    fontweight="bold"
)

for ax in axes.flatten():
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.savefig(
    FIG / "Figure2_Complete.pdf",
    bbox_inches="tight"
)

fig.savefig(
    FIG / "Figure2_Complete.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print("\nFIGURE 2 COMPLETE")
print("Discovery significant genes:", int(d_sig.sum()))
print("Validation significant genes:", int(v_sig.sum()))
print(
    "Cross-dataset significant genes:",
    int(shared["significant_both"].sum())
)
print("Saved Figure2_Complete.pdf and PNG")
