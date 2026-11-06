# Factorial analysis — Stage B

Design: 2 (harmonization) x 7 (normalization), per population pair (10 pairs).
The "missing genotypes" axis is left out: it only has 1 observable level in
this data (see `02_harmonization_configs/config.yaml`), and an axis with no
variation can't be factorialized.

## Question

Does how much a normalization reduces the separation between populations
depend on which harmonization configuration was used upstream? Tested as an
interaction:

```
normalization_effect(config) = D_g,h[normalization] − D_g,h[raw SUM]  (within the same config)
interaction = normalization_effect(Config 2) − normalization_effect(Config 1)
```

## Result

| normalization | mean \|interaction\| | max \|interaction\| |
|---|---|---|
| `Z_global` | 0.0000 | 0.0000 |
| `Z_EUR` | 0.0000 | 0.0000 |
| `Z_percop` | 0.1560 | 0.3361 |
| `Z_MostSimilarPop` | 0.1564 | 0.3361 |
| `Z_norm1` | 0.1582 | 0.3491 |
| `Z_norm2` | 0.1583 | 0.3491 |

**The simple affine transforms (`Z_global`, `Z_EUR`) have exactly zero
interaction** — this re-confirms, from a different angle, the mathematical
property from Section 4: since they are the same linear transformation
applied to everyone, "how much they change D_g,h" cannot depend on anything
beyond their own scale, so it is always zero, under any harmonization
configuration.

**Ancestry-aware normalizations (pgsc_calc's native ones, and `Z_percop`)
show a real interaction.** The most extreme case is the EUR-SAS pair: at
baseline the separation was already nearly zero (D=0.008), so normalization
had little to "flatten"; under Config 2 (ambiguous included) a real
separation appears (D=-0.218), and there normalization does real work. The
same mechanism shows up in AFR-EAS and AMR-EAS.

## Interpretation for the project

This refines the conclusion from `PROJECT_PLAN.md` Section 8: ancestry-aware
normalization is not a universal correction independent of harmonization —
**how much it "hides" a prior technical choice depends on that choice**. Two
people who disagree only on whether to include 30 ambiguous SNPs can reach
different conclusions about "how much the normalization fixed the problem",
even using the same normalization method.

## Files

- `factorial_Dgh.tsv` — full 2x7 grid, D_g,h per population pair.
- `factorial_interaction.tsv` — interaction row by row (10 pairs x 6 metrics, SUM excluded since it's the reference).
- `interaction_summary.tsv` — the table above.
