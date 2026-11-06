"""Stage B -- factorial analysis: harmonization (2 levels) x normalization (7 levels).

The "missing genotypes" axis is excluded: 1000G Phase 3 has no missing
genotype at all in any of the ~59k extracted positions (see
02_harmonization_configs/config.yaml), so that axis only has 1 observable
level here, not 2 -- an axis with no variation cannot be factorialized.

Question: does the effect of a normalization choice on between-population
separation (D_g,h) depend on which harmonization configuration was used?
Tested as an interaction: how much "how much normalization changes D_g,h"
differs between Configuration 1 (ambiguous variants excluded) and
Configuration 2 (ambiguous variants included).

Inputs:
  01_baseline_prs/pgsc_calc_run/kgpphase3/score/kgpphase3_pgs.txt.gz   (harmonization Config 1)
  02_harmonization_configs/keep_ambiguous/kgpphase3/score/kgpphase3_pgs.txt.gz  (Config 2)

Outputs (04_factorial/):
  factorial_Dgh.tsv         -- D_g,h for all 2x7 combinations, per population pair
  factorial_interaction.tsv -- harmonization x normalization interaction, per pair and metric
  interaction_summary.tsv   -- mean |interaction| per metric, ranked
"""
import csv
import gzip
import logging
import statistics
from itertools import combinations
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[STEP] %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "04_factorial"
OUTDIR.mkdir(exist_ok=True)

METRICS = ["SUM", "Z_global", "Z_EUR", "Z_percop", "Z_MostSimilarPop", "Z_norm1", "Z_norm2"]
HARM_CONFIGS = {
    "config1_ambiguous_excluded": ROOT / "01_baseline_prs/pgsc_calc_run/kgpphase3/score/kgpphase3_pgs.txt.gz",
    "config2_ambiguous_included": ROOT / "02_harmonization_configs/keep_ambiguous/kgpphase3/score/kgpphase3_pgs.txt.gz",
}


def load_panel():
    panel = {}
    with open(ROOT / "resources/integrated_call_samples_v3.20130502.ALL.panel") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            panel[row["sample"]] = (row["pop"], row["super_pop"])
    return panel


def load_pgs(path, panel):
    rows = []
    with gzip.open(path, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["sampleset"] != "kgpphase3":
                continue
            iid = row["IID"]
            pop, sp = panel.get(iid, ("NA", "NA"))
            rows.append({
                "sample_id": iid, "pop": pop, "superpop": sp,
                "SUM": float(row["SUM"]),
                "Z_MostSimilarPop": float(row["Z_MostSimilarPop"]),
                "Z_norm1": float(row["Z_norm1"]),
                "Z_norm2": float(row["Z_norm2"]),
            })
    return rows


def add_affine_normalizations(rows):
    """Z_global, Z_EUR, Z_percop -- simple affine transforms (not from pgsc_calc,
    these are the ones from PROJECT_PLAN.md Section 4), computed here."""
    sums = [r["SUM"] for r in rows]
    mu_g, sd_g = statistics.mean(sums), statistics.stdev(sums)
    eur = [r["SUM"] for r in rows if r["superpop"] == "EUR"]
    mu_e, sd_e = statistics.mean(eur), statistics.stdev(eur)
    by_pop = {}
    for r in rows:
        by_pop.setdefault(r["superpop"], []).append(r["SUM"])
    stats_pop = {sp: (statistics.mean(v), statistics.stdev(v)) for sp, v in by_pop.items()}
    for r in rows:
        r["Z_global"] = (r["SUM"] - mu_g) / sd_g
        r["Z_EUR"] = (r["SUM"] - mu_e) / sd_e
        mu_p, sd_p = stats_pop[r["superpop"]]
        r["Z_percop"] = (r["SUM"] - mu_p) / sd_p
    return rows


def cohens_d(a, b):
    na, nb = len(a), len(b)
    va, vb = statistics.variance(a), statistics.variance(b)
    pooled = (((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)) ** 0.5
    return (statistics.mean(a) - statistics.mean(b)) / pooled


def build_D(rows, superpops):
    by_sp = {sp: {m: [r[m] for r in rows if r["superpop"] == sp] for m in METRICS} for sp in superpops}
    return {
        f"{g}-{h}": {m: cohens_d(by_sp[g][m], by_sp[h][m]) for m in METRICS}
        for g, h in combinations(superpops, 2)
    }


def main():
    panel = load_panel()
    log.info("Loading both harmonization configurations")
    configs = {}
    for name, path in HARM_CONFIGS.items():
        rows = add_affine_normalizations(load_pgs(path, panel))
        configs[name] = rows
        log.info(f"  {name}: {len(rows)} samples")

    superpops = sorted(set(r["superpop"] for r in configs["config1_ambiguous_excluded"]))
    D = {name: build_D(rows, superpops) for name, rows in configs.items()}
    pairs = list(D["config1_ambiguous_excluded"].keys())

    log.info("Writing factorial_Dgh.tsv (full 2-harmonization x 7-normalization grid)")
    with open(OUTDIR / "factorial_Dgh.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["harmonization", "pair"] + METRICS)
        for cfg in HARM_CONFIGS:
            for p in pairs:
                w.writerow([cfg, p] + [round(D[cfg][p][m], 4) for m in METRICS])

    log.info("Computing harmonization x normalization interaction")
    # effect of normalization within each harmonization config: D[m] - D[SUM]
    # interaction: difference of that effect between the two harmonization configs
    interaction_rows = []
    for p in pairs:
        for m in METRICS:
            if m == "SUM":
                continue
            eff1 = D["config1_ambiguous_excluded"][p][m] - D["config1_ambiguous_excluded"][p]["SUM"]
            eff2 = D["config2_ambiguous_included"][p][m] - D["config2_ambiguous_included"][p]["SUM"]
            interaction_rows.append({
                "pair": p, "metric": m,
                "normalization_effect_config1": round(eff1, 4),
                "normalization_effect_config2": round(eff2, 4),
                "interaction": round(eff2 - eff1, 4),
            })
    with open(OUTDIR / "factorial_interaction.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(interaction_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(interaction_rows)

    log.info("Summarizing mean interaction per metric (mean |interaction| over the 10 pairs)")
    summary = {}
    for row in interaction_rows:
        summary.setdefault(row["metric"], []).append(abs(row["interaction"]))
    summary_rows = sorted(
        ({"metric": m, "mean_abs_interaction": round(statistics.mean(v), 4),
          "max_abs_interaction": round(max(v), 4)} for m, v in summary.items()),
        key=lambda r: -r["mean_abs_interaction"],
    )
    with open(OUTDIR / "interaction_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["metric", "mean_abs_interaction", "max_abs_interaction"], delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    log.info("Summary (mean |interaction|, largest first):")
    for row in summary_rows:
        log.info(f"  {row['metric']:20s} mean={row['mean_abs_interaction']:.4f}  max={row['max_abs_interaction']:.4f}")


if __name__ == "__main__":
    main()
