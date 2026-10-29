from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
OUT = ROOT / "reproduction_test" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

def read(name):
    path = SRC / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("CREATED:", name)

# Figure 1: cohort composition
cohorts = read("Figure1_cohort_summary.csv")
required = {
    "dataset", "donors", "AD_donors",
    "control_donors", "endothelial_cells"
}
assert required.issubset(cohorts.columns)
assert (cohorts["AD_donors"] +
        cohorts["control_donors"] == cohorts["donors"]).all()

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
x = np.arange(len(cohorts))
axes[0].bar(x, cohorts["control_donors"], label="Control")
axes[0].bar(
    x, cohorts["AD_donors"],
    bottom=cohorts["control_donors"], label="AD"
)
axes[0].set_xticks(x, cohorts["dataset"], rotation=15)
axes[0].set_ylabel("Donors")
axes[0].legend()

axes[1].bar(x, cohorts["endothelial_cells"])
axes[1].set_xticks(x, cohorts["dataset"], rotation=15)
axes[1].set_ylabel("Endothelial cells")
save(fig, "Figure1_release_reconstruction")

# Figure 2: discovery versus validation effect sizes
d = read("Figure2A_discovery_DE.csv")
v = read("Figure2B_validation_DE.csv")
shared = read("Figure2C_shared_gene_effects.csv")

assert d["gene"].is_unique and v["gene"].is_unique
assert shared["gene"].is_unique
assert len(shared) == len(set(d.gene) & set(v.gene))

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(
    shared["logFC_discovery"],
    shared["logFC_validation"],
    s=8, alpha=0.35
)
lim = np.nanmax(np.abs(shared[
    ["logFC_discovery", "logFC_validation"]
].to_numpy()))
ax.plot([-lim, lim], [-lim, lim], linestyle="--")
ax.axhline(0, linewidth=0.6)
ax.axvline(0, linewidth=0.6)
ax.set(
    xlabel="Discovery log fold change",
    ylabel="Validation log fold change",
    title="Cross-cohort gene-effect comparison"
)
save(fig, "Figure2_release_reconstruction")

# Figure 3: pathway NES agreement
p = read("Figure3C_all_pathway_replication.csv")
assert p["pathway"].is_unique
assert len(p) == 1168

expected = (
    (p["padj_discovery"] < 0.05)
    & (p["padj_validation"] < 0.05)
    & (np.sign(p["NES_discovery"]) ==
       np.sign(p["NES_validation"]))
)
recorded = p["replicated_corrected"].astype(str).str.lower().eq("true")
assert (expected == recorded).all(), "Replication flags inconsistent"

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(
    p["NES_discovery"], p["NES_validation"],
    s=10, alpha=0.35
)
ax.axhline(0, linewidth=0.6)
ax.axvline(0, linewidth=0.6)
hits = p.loc[recorded]
ax.scatter(
    hits["NES_discovery"], hits["NES_validation"],
    s=45, marker="*", label="Replicated at FDR < 0.05"
)
for _, row in hits.iterrows():
    ax.annotate(
        row["pathway"],
        (row["NES_discovery"], row["NES_validation"]),
        fontsize=7
    )
ax.set(
    xlabel="Discovery normalized enrichment score",
    ylabel="Validation normalized enrichment score",
    title="Cross-cohort pathway replication"
)
ax.legend()
save(fig, "Figure3_release_reconstruction")

# Figure 4: robustness FDR across analysis models
r = read("Figure4_robustness_input.csv")
assert len(r) == 12
assert not r.duplicated(["label", "model"]).any()

fdr = pd.to_numeric(r["FDR"], errors="raise")
recorded = r["significant"].astype(str).str.lower().eq("true")
assert (recorded == (fdr < 0.05)).all(), "FDR flags inconsistent"

fig, ax = plt.subplots(figsize=(10, 5))
for label, group in r.groupby("label"):
    group = group.reset_index(drop=True)
    positions = np.arange(len(group))
    ax.scatter(
        positions,
        -np.log10(np.maximum(group["FDR"], 1e-300)),
        label=label, s=45
    )
    for i, row in group.iterrows():
        ax.annotate(
            row["model"].replace("_", "\n"),
            (i, -np.log10(max(row["FDR"], 1e-300))),
            fontsize=6, rotation=35
        )
ax.axhline(-np.log10(0.05), linestyle="--",
           label="FDR = 0.05")
ax.set_ylabel("-log10(FDR)")
ax.set_title("Pathway robustness across analysis models")
ax.legend()
save(fig, "Figure4_release_reconstruction")

print("\n=== SOURCE FILE CHECKSUMS ===")
for path in sorted(SRC.glob("*.csv")):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    print(digest, path.name)

print("\nPASS: Four diagnostic figures generated from packaged CSVs.")
print("NOTE: These are independent reconstructions, not pixel-identical originals.")
print("Output directory:", OUT)
