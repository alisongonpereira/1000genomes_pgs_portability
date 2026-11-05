"""Etapa B -- analise fatorial: harmonizacao (2 niveis) x normalizacao (7 niveis).

O eixo "genotipos ausentes" nao entra: 1000G Phase 3 nao tem genotipo ausente
em nenhuma das ~59k posicoes extraidas (ver 02_harmonization_configs/config.yaml),
entao so tem 1 nivel observavel ali, nao 2 -- nao da pra fatorializar um eixo
sem variacao.

Pergunta: o efeito de uma escolha de normalizacao sobre a separacao entre
populacoes (D_g,h) depende de qual configuracao de harmonizacao foi usada?
Isso e testado como interacao: quanto o "quanto a normalizacao muda D_g,h"
difere entre Config 1 (ambiguas excluidas) e Config 2 (ambiguas incluidas).

Entradas:
  01_baseline_prs/pgsc_calc_run/kgpphase3/score/kgpphase3_pgs.txt.gz   (harmonizacao Config 1)
  02_harmonization_configs/keep_ambiguous/kgpphase3/score/kgpphase3_pgs.txt.gz  (Config 2)

Saidas (04_factorial/):
  factorial_Dgh.tsv       -- D_g,h para as 2x7 combinacoes, por par de populacao
  factorial_interaction.tsv -- interacao harmonizacao x normalizacao, por par e metrica
  interaction_summary.tsv -- |interacao| media por metrica, ranqueada
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
    "config1_ambig_excl": ROOT / "01_baseline_prs/pgsc_calc_run/kgpphase3/score/kgpphase3_pgs.txt.gz",
    "config2_ambig_incl": ROOT / "02_harmonization_configs/keep_ambiguous/kgpphase3/score/kgpphase3_pgs.txt.gz",
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
    """Z_global, Z_EUR, Z_percop -- transformacoes afins simples, calculadas aqui
    (nao vem do pgsc_calc, sao as da Secao 4 do PROJECT_PLAN.md)."""
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
    log.info("Carregando as duas configuracoes de harmonizacao")
    configs = {}
    for name, path in HARM_CONFIGS.items():
        rows = add_affine_normalizations(load_pgs(path, panel))
        configs[name] = rows
        log.info(f"  {name}: {len(rows)} amostras")

    superpops = sorted(set(r["superpop"] for r in configs["config1_ambig_excl"]))
    D = {name: build_D(rows, superpops) for name, rows in configs.items()}
    pairs = list(D["config1_ambig_excl"].keys())

    log.info("Escrevendo factorial_Dgh.tsv (grade completa 2 harmonizacoes x 7 normalizacoes)")
    with open(OUTDIR / "factorial_Dgh.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["harmonizacao", "pair"] + METRICS)
        for cfg in HARM_CONFIGS:
            for p in pairs:
                w.writerow([cfg, p] + [round(D[cfg][p][m], 4) for m in METRICS])

    log.info("Calculando interacao harmonizacao x normalizacao")
    # efeito da normalizacao dentro de cada harmonizacao: D[m] - D[SUM]
    # interacao: diferenca desse efeito entre as duas harmonizacoes
    interaction_rows = []
    for p in pairs:
        for m in METRICS:
            if m == "SUM":
                continue
            eff1 = D["config1_ambig_excl"][p][m] - D["config1_ambig_excl"][p]["SUM"]
            eff2 = D["config2_ambig_incl"][p][m] - D["config2_ambig_incl"][p]["SUM"]
            interaction_rows.append({
                "pair": p, "metric": m,
                "efeito_normalizacao_config1": round(eff1, 4),
                "efeito_normalizacao_config2": round(eff2, 4),
                "interacao": round(eff2 - eff1, 4),
            })
    with open(OUTDIR / "factorial_interaction.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(interaction_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(interaction_rows)

    log.info("Resumindo interacao media por metrica (|interacao| media sobre os 10 pares)")
    summary = {}
    for row in interaction_rows:
        summary.setdefault(row["metric"], []).append(abs(row["interacao"]))
    summary_rows = sorted(
        ({"metric": m, "mean_abs_interaction": round(statistics.mean(v), 4),
          "max_abs_interaction": round(max(v), 4)} for m, v in summary.items()),
        key=lambda r: -r["mean_abs_interaction"],
    )
    with open(OUTDIR / "interaction_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["metric", "mean_abs_interaction", "max_abs_interaction"], delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    log.info("Resumo (|interacao| media, maior primeiro):")
    for row in summary_rows:
        log.info(f"  {row['metric']:20s} mean={row['mean_abs_interaction']:.4f}  max={row['max_abs_interaction']:.4f}")


if __name__ == "__main__":
    main()
