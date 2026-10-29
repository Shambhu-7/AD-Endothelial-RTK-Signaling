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

path = (
    ROOT / "manuscript_final" / "corrected_final" / "source_data" / "Figure4_robustness_input.csv"
)

df = pd.read_csv(path)

models = [
    "discovery_primary",
    "discovery_lowcell",
    "discovery_composition",
    "validation_primary",
    "validation_lowcell",
    "validation_composition"
]

model_names = [
    "Discovery: Primary",
    "Discovery: Low-cell exclusion",
    "Discovery: Composition adjusted",
    "Validation: Primary",
    "Validation: Low-cell exclusion",
    "Validation: Composition adjusted"
]

pathways = ["RTK", "mTORC1"]

sub = df[
    df["label"].isin(pathways) &
    df["model"].isin(models)
].copy()

assert len(sub) == 12, "Expected exactly 12 pathway/model rows"
assert not sub.duplicated(["label", "model"]).any()

sub["NES"] = pd.to_numeric(sub["NES"], errors="raise")
sub["FDR"] = pd.to_numeric(sub["FDR"], errors="raise")

assert sub["NES"].notna().all()
assert sub["FDR"].between(0, 1).all()

# Recompute significance directly from FDR.
sub["significant"] = sub["FDR"] < 0.05

assert sub["available"].astype(str).str.lower().eq("true").all()
assert sub.groupby("label")["significant"].sum().to_dict() == {
    "RTK": 6,
    "mTORC1": 4
}

sub["model"] = pd.Categorical(
    sub["model"],
    categories=models,
    ordered=True
)

sub["label"] = pd.Categorical(
    sub["label"],
    categories=pathways,
    ordered=True
)

sub = sub.sort_values(["model", "label"])

sub.to_csv(
    SRC / "Figure4_robustness_source_data.csv",
    index=False
)

nes = sub.pivot(
    index="model",
    columns="label",
    values="NES"
).reindex(index=models, columns=pathways)

fdr = sub.pivot(
    index="model",
    columns="label",
    values="FDR"
).reindex(index=models, columns=pathways)

sig = fdr < 0.05

fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1.15, 1],
    hspace=0.46,
    wspace=0.40
)

ax_a = fig.add_subplot(gs[0, :])
ax_b = fig.add_subplot(gs[1, 0])
ax_c = fig.add_subplot(gs[1, 1])

# Panel a: NES heatmap
heat = ax_a.imshow(
    nes.to_numpy(),
    aspect="auto",
    cmap="YlGnBu",
    vmin=1.4,
    vmax=2.4
)

ax_a.set_xticks(range(2))
ax_a.set_xticklabels([
    "RTK signaling",
    "mTORC1 signaling"
])

ax_a.set_yticks(range(6))
ax_a.set_yticklabels(model_names)

for i in range(6):
    for j in range(2):
        value = nes.iloc[i, j]
        q = fdr.iloc[i, j]
        symbol = "*" if q < 0.05 else "ns"

        ax_a.text(
            j, i,
            f"NES {value:.2f}\n"
            f"FDR {q:.2g} ({symbol})",
            ha="center",
            va="center",
            fontsize=9,
            color="black"
        )

ax_a.axhline(2.5, color="white", lw=3)

fig.colorbar(
    heat,
    ax=ax_a,
    fraction=0.025,
    pad=0.02,
    label="NES"
)

ax_a.set_title(
    "a  Pathway enrichment across six model specifications",
    loc="left",
    fontweight="bold"
)

# Panel b: NES across models
x = np.arange(6)

for label, color in [
    ("RTK", "#277DA1"),
    ("mTORC1", "#F8961E")
]:
    y = nes[label].to_numpy()

    ax_b.plot(
        x, y,
        marker="o",
        linewidth=1.8,
        label=label,
        color=color
    )

# Mark every nonsignificant result with a hollow marker.
for label, color in [("RTK", "#277DA1"), ("mTORC1", "#F8961E")]:
    for model in models:
        if not sig.loc[model, label]:
            ax_b.scatter([models.index(model)], [nes.loc[model, label]],
                         facecolors="white", edgecolors=color,
                         s=110, linewidths=2, zorder=5)

ax_b.set_xticks(x)
ax_b.set_xticklabels(
    ["D-P", "D-L", "D-C", "V-P", "V-L", "V-C"]
)

ax_b.set_ylabel("Normalized enrichment score")
ax_b.set_ylim(1.3, 2.5)
ax_b.legend(frameon=False)
ax_b.set_title(
    "b  Sensitivity of enrichment magnitude",
    loc="left",
    fontweight="bold"
)

ax_b.text(
    0.02, -0.23,
    "D/V: discovery/validation; "
    "P: primary; L: low-cell exclusion; "
    "C: composition adjusted",
    transform=ax_b.transAxes,
    fontsize=8
)

# Panel c: robustness counts
counts = [
    int(sig["RTK"].sum()),
    int(sig["mTORC1"].sum())
]

bars = ax_c.bar(
    ["RTK signaling", "mTORC1 signaling"],
    counts,
    color=["#277DA1", "#F8961E"]
)

for bar, n in zip(bars, counts):
    ax_c.text(
        bar.get_x() + bar.get_width()/2,
        n + 0.10,
        f"{n}/6",
        ha="center",
        fontsize=12,
        fontweight="bold"
    )

ax_c.set_ylim(0, 7)
ax_c.set_ylabel("FDR-significant specifications")
ax_c.set_title(
    "c  Robustness summary",
    loc="left",
    fontweight="bold"
)

for ax in [ax_b, ax_c]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.savefig(
    FIG / "Figure4_Complete.pdf",
    bbox_inches="tight"
)

fig.savefig(
    FIG / "Figure4_Complete.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print("\nFIGURE 4 COMPLETE")
print("RTK significant models:", counts[0], "/6")
print("mTORC1 significant models:", counts[1], "/6")
print(
    "Discovery composition-adjusted mTORC1 FDR:",
    fdr.loc["discovery_composition", "mTORC1"]
)
print("Saved Figure4_Complete.pdf and PNG")
