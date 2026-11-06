"""High-visibility map figure: mean PRS313_BC raw score per 1000 Genomes
collection population, plotted at each population's real sampling location.

This is a deliberately more "shareable" companion to the panel figure in
03_figures/sensitivity_panel.png -- same underlying baseline data
(01_baseline_prs/baseline_prs_scores.tsv, Config 1 / keep_ambiguous=false),
just a different, more visually arresting chart form for a general audience.

IMPORTANT: this plots the raw, uncalibrated SUM score. It is explicitly
NOT a disease-risk map -- the whole point of this project is that this
pattern is sensitive to normalization/harmonization choices (see
01_baseline_prs/, 02_harmonization_configs/, 04_factorial/). The figure
carries a prominent on-image disclaimer for this reason; do not strip it
out when reusing the figure elsewhere.

Coordinates are approximate representative points for each population's
real 1000 Genomes collection site (city-level), not individual sample
GPS data -- see POP_COORDS below for sources (standard IGSR/1000G
population descriptions).

Output: 03_figures/world_map_prs.png (and .pdf)
"""
import logging
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from shapely.geometry import box

# Cut the map in the Atlantic instead of the Pacific: standard maps (cut at
# +-180) read AMR | AFR/EUR | ASIA left-to-right. Cutting here instead gives
# AFR/EUR | ASIA | AMR -- continents in one clean sweep, since that's the
# actual eastward order once you stop breaking the circle at the Pacific.
CUT_LON = -25.0


def recenter_lon(lon):
    return lon if lon >= CUT_LON else lon + 360


def recenter_world(world):
    """Split polygons that straddle CUT_LON before shifting the western
    piece by +360 -- a naive coordinate shift alone would draw a stray
    horizontal line across any landmass that crosses the new antimeridian
    (here, just Greenland/Iceland)."""
    east = gpd.clip(world, box(CUT_LON, -90, 180, 90))
    west = gpd.clip(world, box(-180, -90, CUT_LON, 90))
    west = west.copy()
    west["geometry"] = west.translate(xoff=360)
    return pd.concat([east, west], ignore_index=True)

logging.basicConfig(level=logging.INFO, format="[STEP] %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "03_figures"
OUTDIR.mkdir(exist_ok=True)

# lat, lon, superpop, display label -- approximate city-level point for each
# population's real 1000 Genomes collection site (official IGSR population
# descriptions, e.g. "Gujarati Indian from Houston, Texas",
# "Sri Lankan Tamil from the UK"). Diaspora populations are placed at their
# actual collection site, not an inferred ancestral homeland.
POP_COORDS = {
    "YRI": (7.3775, 3.9470, "AFR", "Yoruba, Ibadan (Nigeria)"),
    "LWK": (0.6101, 34.7676, "AFR", "Luhya, Webuye (Kenya)"),
    "GWD": (13.4549, -16.5790, "AFR", "Gambian, Western Division (Gambia)"),
    "MSL": (7.9644, -11.7383, "AFR", "Mende (Sierra Leone)"),
    "ESN": (6.7333, 6.1333, "AFR", "Esan (Nigeria)"),
    "ASW": (33.4484, -112.0740, "AFR", "African ancestry, SW USA"),
    "ACB": (13.1132, -59.5988, "AFR", "African Caribbean (Barbados)"),
    "CEU": (40.7608, -111.8910, "EUR", "N/W European ancestry, Utah (USA)"),
    "TSI": (43.7696, 11.2558, "EUR", "Toscani (Italy)"),
    "FIN": (60.1699, 24.9384, "EUR", "Finnish (Finland)"),
    "GBR": (51.5072, -0.1276, "EUR", "British (England/Scotland)"),
    "IBS": (40.4168, -3.7038, "EUR", "Iberian (Spain)"),
    # STU/ITU's official descriptor is just "from the UK" (no specific city),
    # so Cardiff/Edinburgh are illustrative stand-ins chosen to keep this UK
    # trio legible on a world-scale map -- GBR itself keeps its real London point.
    "CHB": (39.9042, 116.4074, "EAS", "Han Chinese, Beijing (China)"),
    "JPT": (35.6762, 139.6503, "EAS", "Japanese, Tokyo (Japan)"),
    "CHS": (23.1291, 113.2644, "EAS", "Southern Han Chinese (China)"),
    "CDX": (22.0017, 100.7971, "EAS", "Chinese Dai, Xishuangbanna (China)"),
    "KHV": (10.8231, 106.6297, "EAS", "Kinh, Ho Chi Minh City (Vietnam)"),
    "GIH": (29.7604, -95.3698, "SAS", "Gujarati Indian, Houston (USA)"),
    "PJL": (31.5497, 74.3436, "SAS", "Punjabi, Lahore (Pakistan)"),
    "BEB": (23.8103, 90.4125, "SAS", "Bengali (Bangladesh)"),
    "STU": (51.4816, -3.1791, "SAS", "Sri Lankan Tamil, UK"),
    "ITU": (55.9533, -3.1883, "SAS", "Indian Telugu, UK"),
    "MXL": (34.0522, -118.2437, "AMR", "Mexican ancestry, Los Angeles (USA)"),
    "PUR": (18.4655, -66.1057, "AMR", "Puerto Rican (Puerto Rico)"),
    "CLM": (6.2442, -75.5812, "AMR", "Colombian, Medellín (Colombia)"),
    "PEL": (-12.0464, -77.0428, "AMR", "Peruvian, Lima (Peru)"),
}

SUPERPOP_MARKERS = {"AFR": "o", "AMR": "s", "EAS": "^", "EUR": "D", "SAS": "*"}
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"


def load_scores():
    log.info("Loading baseline scores (Config 1, keep_ambiguous=false) and aggregating per population")
    df = pd.read_csv(ROOT / "01_baseline_prs/baseline_prs_scores.tsv", sep="\t")
    agg = df.groupby("pop")["raw_score_sum"].agg(["mean", "count"]).reset_index()
    agg.columns = ["pop", "mean_sum", "n"]
    rows = []
    for _, r in agg.iterrows():
        if r["pop"] not in POP_COORDS:
            log.info(f"  skipping {r['pop']} -- no coordinate on file")
            continue
        lat, lon, sp, label = POP_COORDS[r["pop"]]
        rows.append({"pop": r["pop"], "superpop": sp, "label": label,
                      "lat": lat, "lon": lon, "mean_sum": r["mean_sum"], "n": int(r["n"])})
    return pd.DataFrame(rows)


def main():
    data = load_scores()
    log.info(f"{len(data)} populations plotted, global mean SUM = {data['mean_sum'].mean():.4f}")

    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    world = recenter_world(world)
    data["lon"] = data["lon"].apply(recenter_lon)

    fig = plt.figure(figsize=(16, 10.5))
    fig.patch.set_facecolor("#fcfcfb")
    gs = fig.add_gridspec(2, 3, height_ratios=[11, 1], hspace=0.06, wspace=0.15,
                           left=0.02, right=0.98, top=0.90, bottom=0.04)
    ax = fig.add_subplot(gs[0, :])
    ax_shape = fig.add_subplot(gs[1, 0])
    ax_size = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 2])
    ax.set_facecolor("#fcfcfb")

    world.plot(ax=ax, color="#eceae3", edgecolor="#d8d5cb", linewidth=0.5, zorder=1)

    vmax = float(np.abs(data["mean_sum"] - data["mean_sum"].mean()).max())
    norm = plt.Normalize(vmin=data["mean_sum"].mean() - vmax, vmax=data["mean_sum"].mean() + vmax)
    cmap = plt.get_cmap("RdBu_r")

    for sp, marker in SUPERPOP_MARKERS.items():
        sub = data[data["superpop"] == sp]
        ax.scatter(
            sub["lon"], sub["lat"],
            s=sub["n"] * 3.2, c=sub["mean_sum"], cmap=cmap, norm=norm,
            marker=marker, edgecolors="white", linewidths=1.1, zorder=3, alpha=0.92,
        )

    # per-population label offsets (points), to dodge the tight UK and
    # Nigeria clusters -- default is straight above the marker
    LABEL_OFFSETS = {
        "YRI": (-14, 4), "ESN": (14, -6),
        "GBR": (0, -13), "STU": (0, 8), "ITU": (0, 8),
        "ASW": (0, -13), "MXL": (-4, 8),
    }
    for _, r in data.iterrows():
        dx, dy = LABEL_OFFSETS.get(r["pop"], (0, 7))
        ax.annotate(
            r["pop"], (r["lon"], r["lat"]), textcoords="offset points", xytext=(dx, dy),
            fontsize=7.5, ha="center", color=TEXT_SECONDARY, zorder=4,
        )

    ax.set_xlim(CUT_LON - 3, CUT_LON + 335)
    ax.set_ylim(-45, 72)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Legends live in their own horizontal strip below the map, not on top
    # of it, so the map stays the visual focus: superpop/shape, then size,
    # then the color gradient, left to right.
    for a in (ax_shape, ax_size, ax_cbar):
        a.axis("off")

    shape_handles = [
        Line2D([0], [0], marker=m, color="w", markerfacecolor="#898781", markeredgecolor="#898781", markersize=10, label=sp)
        for sp, m in SUPERPOP_MARKERS.items()
    ]
    ax_shape.legend(handles=shape_handles, title="Superpopulation", loc="center", ncol=len(shape_handles),
                     frameon=False, fontsize=9.5, title_fontsize=9.5, handletextpad=0.5, columnspacing=1.2)

    size_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#898781", markeredgecolor="#898781",
               alpha=0.7, markersize=np.sqrt(n * 3.2), label=f"n={n}")
        for n in [60, 100, 113]
    ]
    ax_size.legend(handles=size_handles, title="Sample size", loc="center", ncol=len(size_handles),
                    frameon=False, fontsize=9.5, title_fontsize=9.5, handletextpad=0.8, columnspacing=1.6,
                    borderpad=1.4)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_inset = ax_cbar.inset_axes([0.08, 0.5, 0.84, 0.28])
    cbar = fig.colorbar(sm, cax=cbar_inset, orientation="horizontal")
    cbar.set_label("Mean PRS313_BC raw score (SUM)", color=TEXT_PRIMARY, fontsize=9.5)
    cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)

    fig.suptitle(
        "PRS313_BC (PGS000004, breast cancer) — mean raw score across 1000 Genomes populations\n"
        "Populations per the 1000 Genomes Project's own classification — full analysis: "
        "github.com/alisongonpereira/1000genomes_pgs_portability",
        fontsize=14, color=TEXT_PRIMARY, ha="center", y=0.985, linespacing=1.6,
    )

    out_png = OUTDIR / "world_map_prs.png"
    out_pdf = OUTDIR / "world_map_prs.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, bbox_inches="tight", facecolor=fig.get_facecolor())
    log.info(f"Figure saved to {out_png} and {out_pdf}")


if __name__ == "__main__":
    main()
