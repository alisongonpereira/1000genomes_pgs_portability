# LinkedIn post — final, grounded in real results

Supersedes the hypothetical draft in `linkedin_2026-11-05.md` (written before
the analysis existed). Every number below is pulled directly from this
repo's committed outputs — see the file path noted after each claim.

**Attach, in this order:**
1. `03_figures/world_map_prs.png` (hero image)
2. `03_figures/sensitivity_panel.png` (supporting evidence)

---

## Post text

Does a polygenic risk score mean the same thing in every population? I ran a small, fully reproducible audit to find out — PRS313_BC (PGS000004, breast cancer), scored on all 2,504 samples of 1000 Genomes Phase 3.

Not testing predictive accuracy. Just one question: how much does the *distribution* of the score move when you change processing choices that are usually invisible in a methods section?

Three findings.

**1. Normalization choice matters less than you'd expect — mathematically.**
Standardizing by the global mean/SD, or by one reference population's mean/SD, gives the exact same standardized difference between any two populations — to the 4th decimal place, across all 10 population pairs. That's not a result, it's algebra: both are the same linear transform applied to everyone, so they can't change a standardized contrast.

**2. Harmonization choice matters more than you'd expect — empirically.**
PRS313_BC has 30 strand-ambiguous SNPs (A/T, C/G) that a strict pipeline drops by default. I tested the assumption needed to keep them instead (same strand orientation) against real allele frequencies — not an assumption, a check: 30/30 consistent with 1000G's European subset, largest deviation 0.039. Keeping those 30 variants flips the sign of the EUR–SAS contrast: from essentially nothing (standardized difference 0.008) to a real, opposite-direction separation (−0.218).

**3. High correlation is not the same thing as individual stability.**
Across those same two harmonization choices: Spearman ρ = 0.954 on the full cohort — and 52% of the 2,504 individuals change decile.

A flattened map isn't a cured map, either: ancestry-aware normalization (PCA regression or population-specific standardization) pulls every group's mean back toward zero, but *how much* flattening it does depends on the harmonization decision made upstream — tested directly as an interaction, mean |effect| ≈ 0.16, up to 0.35, for every ancestry-aware method. The simple global/reference-population normalizations show zero interaction, for the same algebraic reason as point 1.

None of this is a claim about who has "more risk." The map shows a raw, uncalibrated score at each population's real 1000 Genomes collection site — not disease incidence, and not something you can read off and compare across ancestries without recalibration. If anything, it's an argument for being precise about what a PRS number actually encodes before it gets anywhere near a clinical decision.

Full pipeline, data, and code: github.com/alisongonpereira/1000genomes_pgs_portability

#StatisticalGenetics #PolygenicRiskScores #Bioinformatics #ComputationalGenomics

---

## Source for each number

| Claim | File |
|---|---|
| 4th-decimal match, Z_global/Z_EUR vs raw SUM, all 10 pairs | `01_baseline_prs/sensitivity_Dgh.tsv` |
| 30 ambiguous SNPs, 30/30 frequency-consistent, max delta 0.039 | `02_harmonization_configs/ambiguous_strand_validation.tsv` |
| EUR-SAS: 0.008 -> -0.218 | `02_harmonization_configs/delta_Dgh_keepambig.tsv` |
| Spearman rho = 0.954, 51.9% decile change | `03_figures/sensitivity_panel.png` (Panel C), produced by `scripts/script_1_sensitivity_figure.py` |
| Interaction means (0.156-0.158 ancestry-aware, 0.0 affine) | `04_factorial/interaction_summary.tsv` |
