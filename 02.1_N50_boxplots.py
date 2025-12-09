# 02.1_N50_boxplots.py
# Generate N50 boxplots by habitat and region from QC'd contig summary stats Excel

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
from pathlib import Path

# Load QC'd data
xlsx_path = Path("/mnt/data/contig_summary_stats.xlsx") # Path to the Excel file
xls = pd.ExcelFile(xlsx_path)

def detect_columns(df: pd.DataFrame):
    n50_col = next((c for c in df.columns if "n50" in str(c).lower()), None)
    group_col = next((c for c in df.columns if any(k in str(c).lower() for k in ["habitat","group","environment","biome"])), None)
    region_col = next((c for c in df.columns if any(k in str(c).lower() for k in ["region","zone"])), None)
    return n50_col, group_col, region_col

df_base = None
for s in xls.sheet_names:
    tmp = pd.read_excel(xls, sheet_name=s)
    n50_col, group_col, region_col = detect_columns(tmp)
    if n50_col and group_col and region_col:
        df_base = tmp
        break
if df_base is None:
    raise RuntimeError("Could not auto-detect N50, habitat, and region columns.")

def map_group(val):
    v = str(val).lower()
    if "marine" in v:
        return "Marine"
    if "terrestrial" in v or "land" in v or "soil" in v:
        return "Terrestrial"
    return None

def map_region(val):
    v = str(val).lower()
    if "north" in v:
        return "North"
    if "south" in v:
        return "South"
    if "eq" in v or "equator" in v or "temperate" in v:
        return "Equator"
    return None

df = pd.DataFrame({
    "N50": pd.to_numeric(df_base[n50_col], errors="coerce"),
    "Group": df_base[group_col].map(map_group),
    "Region": df_base[region_col].map(map_region)
}).dropna()

df = df[df["N50"] > 0].copy()

order = ["North", "Equator", "South"]
palette = {"North": "#0072B2", "Equator": "#009E73", "South": "#D55E00"}  # Okabe–Ito

# Jitter defaults
JITTER_ALPHA = 0.30
JITTER_S = 6.0
JITTER_WIDTH = 0.20

# Helper to set sane log ticks given current ylim
def apply_log_ticks(ax, num_major=6):
    # Major ticks at powers of 10, limited to ~num_major
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=num_major))
    # Minor ticks between powers of 10
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=12))
    ax.yaxis.set_minor_formatter(NullFormatter())

svg_out = Path("/mnt/data/n50_boxplot_combined_FINAL_ticks.svg")
fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=150, sharey=True)

for ax, env in zip(axes, ["Terrestrial", "Marine"]):
    sub = df[df["Group"] == env].copy()
    data = [sub.loc[sub["Region"] == r, "N50"].dropna().values for r in order]
    positions = np.arange(1, len(order) + 1)

    # Boxplot
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.6,
        showfliers=False,
        patch_artist=True
    )
    for patch, r in zip(bp["boxes"], order):
        patch.set_facecolor(palette[r])
        patch.set_alpha(0.6)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.0)
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(1.5)
    for w in bp["whiskers"]:
        w.set_color("black"); w.set_linewidth(1.0)
    for cap in bp["caps"]:
        cap.set_color("black"); cap.set_linewidth(1.0)

    # Jitter
    rng = np.random.default_rng(42)
    for i, r in enumerate(order):
        vals = data[i]
        if len(vals) == 0:
            continue
        jx = rng.uniform(low=positions[i] - JITTER_WIDTH, high=positions[i] + JITTER_WIDTH, size=len(vals))
        ax.scatter(
            jx, vals,
            s=JITTER_S,
            c=palette[r],
            alpha=JITTER_ALPHA,
            edgecolors="none",
            zorder=3
        )

    # Axes
    ax.set_yscale("log")
    if len(sub):
        y_min = max(1e-6, sub["N50"].min())
        y_max = sub["N50"].max()
        ax.set_ylim(y_min / 1.4, y_max * 2.0)  # keep extra headroom
    apply_log_ticks(ax, num_major=7)

    ax.set_xlim(0.5, len(order) + 0.5)
    ax.set_xticks(positions)
    ax.set_xticklabels(order)
    ax.set_xlabel("Region")
    ax.set_title(f"{env} (n={len(sub)})", pad=10)
    if env == "Terrestrial":
        ax.set_ylabel("N50 (bp, log scale)")
    ax.grid(which="major", axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()
fig.savefig(svg_out, format="svg", bbox_inches="tight")
plt.close(fig)

svg_out
