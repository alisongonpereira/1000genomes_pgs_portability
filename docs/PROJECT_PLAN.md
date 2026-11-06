# PRS Portability Audit — Revised Project

## Central scientific question

To what extent do technical processing and normalization choices alter the observed differences in polygenic risk score distributions across populations?

This reformulation makes the project methodologically more interesting. Instead of simply showing that populations exhibit distinct PRS distributions, we investigate how much these differences change when we vary specific analytical decisions, holding the data and the genetic model constant.

There is, however, a fundamental distinction: we can quantify the sensitivity of the distributions to technical choices, but we cannot necessarily decompose the observed difference into one percentage that is purely biological and another that is purely technical.

This distinction will be central to the project's credibility.

## 1. Scientific objective

### General objective

Quantify the influence of different genomic processing and normalization strategies on the distributions of the same PRS across genetically diverse populations.

### Specific objectives

1. Assess how decisions about allelic harmonization, variant inclusion, and missing-genotype handling modify individual scores.
2. Compare the effects of global normalization, normalization based on a reference population, and population-specific normalization.
3. Quantify the stability of between-population differences under different analytical configurations.
4. Determine under which circumstances normalization alters the interpretation of results without modifying the underlying genetic information.

**Central hypothesis:** part of the observed variation in PRS distributions depends on analytical decisions, but the effect of these decisions is neither uniform nor necessarily equivalent across populations.

## 2. Data and scope

We keep the originally proposed data.

| Component | Choice |
|---|---|
| PRS | PRS313_BC — PGS000004 |
| Associated phenotype | Breast cancer |
| PRS source | PGS Catalog |
| Genomic data | 1000 Genomes Project — Phase 3 |
| Populations | AFR, AMR, EAS, EUR and SAS |
| Language | Python |
| Analysis | Comparison of analytical configurations |
| Product | Reproducible repository + scientific figure |

The experiment will not evaluate the PRS's ability to predict breast cancer.

The variable of interest is the stability of score distributions under controlled changes to processing.

## 3. Experimental design

I propose organizing the study into three stages.

### Stage A — Build a reference distribution

First, we compute the PRS using a previously defined analytical configuration.

This configuration is our baseline.

The conditions are:

- Same weights file.
- Same set of individuals.
- Same reference genome.
- Validated allelic harmonization.
- Explicit inclusion/exclusion criteria for variants.
- Logging of missing variants and no-call genotypes.

The result is an initial score distribution for each population.

This stage matters because we cannot compare the effect of different pipelines without first defining a methodological reference.

### Stage B — Vary technical decisions one at a time

This is the core of the project.

**Experimental matrix**

| Dimension | Configuration 1 | Configuration 2 |
|---|---|---|
| Harmonization | Non-ambiguous variants only | Additional inclusion of variants resolved by validated criteria |
| Missing genotypes | Exclusion of incomplete loci | Imputation by expected dosage, with an explicit reference |
| Normalization | Global mean and SD | EUR mean and SD |
| Additional normalization | — | Population-specific mean and SD |

Harmonization must preserve correct allele orientation. We should not compare a correct pipeline against a deliberately incorrect one as if they were two equally valid alternatives.

Population-specific normalization is treated separately, as a diagnostic transformation.

**Why?**

Because standardizing each population individually to mean zero and SD one eliminates, by construction, differences in means and variances.

This does not demonstrate that the prior differences were artifacts.

It only demonstrates that the transformation removed those features from the distributions.

This distinction will be one of the project's most important conceptual results.

### Stage C

*(Reserved for the factorial analysis of relevant combinations — see Section 5.)*

## 4. How do we answer the question quantitatively?

We need to define what it means for a distribution to have changed.

I would use four indicators.

| Indicator | What it quantifies | Interpretation |
|---|---|---|
| Standardized difference between means | Separation between populations | Sensitivity of population contrasts |
| Wasserstein distance | Distance between distributions | Changes in the distribution as a whole |
| Spearman correlation | Stability of individual ranking | How much decisions alter individuals' ranking |
| Percentile change | Relative reclassification | How many individuals change position depending on processing |

Distances between distributions should be presented on a scale that is comparable across configurations.

### An important mathematical property

Consider a raw score $S_i$.

Under global normalization:

$$Z_i=\frac{S_i-\mu_{\mathrm{global}}}{\sigma_{\mathrm{global}}}$$

Under normalization based on the EUR population:

$$Z_{i,\mathrm{EUR}}=\frac{S_i-\mu_{\mathrm{EUR}}}{\sigma_{\mathrm{EUR}}}$$

Both are affine transformations common to every individual.

Consequently, when applied to the same set of scores, they preserve individual ranking and the standardized differences between populations.

Therefore, if we use the standardized difference between means as our indicator, these two normalizations will produce the same result.

This is a mathematical property, not an empirical discovery.

By contrast, population-specific normalization:

$$Z_{i,g}=\frac{S_i-\mu_g}{\sigma_g}$$

modifies individuals' relative position when comparing different populations.

This transformation changes the reference frame of the comparison.

The project needs to separate real changes driven by genomic processing from changes that are mathematically expected purely from the choice of scale.

## 5. How do we quantify the contribution of technical decisions?

I propose comparing each configuration against the baseline.

For two populations $g$ and $h$, define:

$$D_{g,h}^{(k)}$$

as the standardized distance between their distributions under analytical configuration $k$.

Sensitivity to the configuration can be expressed as:

$$\Delta D_{g,h}^{(k)} = D_{g,h}^{(k)} - D_{g,h}^{(\mathrm{baseline})}$$

This lets us answer concrete questions:

- How much did including additional, correctly harmonized variants change the separation between AFR and EUR?
- How much did the handling of missing genotypes change the score distribution in AMR?
- What proportion of individuals changed decile when we altered the processing strategy?

To assess the uncertainty of these estimates, we use bootstrap with population-stratified resampling.

**Caveat:** if we modify several parameters simultaneously, their contributions may interact. For this reason, the main analysis should vary one factor at a time, followed by a factorial analysis of the combinations that show relevant effects.

### What can we conclude?

One possible conclusion, if the results support it, would be:

> Technical processing choices accounted for substantial variation in the observed separation between PRS distributions.

But I would avoid:

> Technical artifacts explain X% of ancestry-related PRS differences.

The second statement would require a causal decomposition that this design, on its own, does not allow us to establish.

## 6. Main scientific figure

My proposal is a figure made of three panels.

### Panel A — PRS distributions

Violin plots or density plots showing the distributions by population across different configurations.

The comparison should use appropriate scales to avoid a unit change being interpreted as a substantive change.

### Panel B — Sensitivity of population contrasts

Heatmap with:

- Rows: between-population comparisons.
- Columns: analytical configurations.
- Values: changes in standardized differences relative to the baseline.

This will be the main panel for answering the scientific question.

### Panel C — Individual-level stability

Scatter plot comparing individual scores between two processing configurations.

We complement it with the Spearman correlation and the proportion of individuals who changed decile.

This panel will show whether the technical choices only changed the scale or also altered individuals' ordering.

## 7. English-language publication — revised version

**Suggested title:** How much of the variation in PRS distributions is driven by analytical choices?

**Draft, contingent on completing the project:**

> Polygenic risk scores are often compared across genetically diverse populations. But before interpreting differences in their distributions, we should ask a more fundamental question:
>
> How sensitive are those differences to the way the scores are processed and normalized?
>
> I designed a reproducible analysis using a published breast cancer polygenic risk score (PGS000004) and public genomic data from the 1000 Genomes Project.
>
> Rather than evaluating predictive performance, the analysis focuses on the technical decisions that precede clinical interpretation.
>
> Specifically, I examined how allele harmonization, variant inclusion, missing-genotype handling, and normalization strategies influence score distributions and between-population comparisons.
>
> [Insert the main quantitative finding after completing the analysis.]
>
> One distinction is particularly important: changing the normalization reference does not necessarily change the underlying separation between populations. Some transformations preserve standardized differences, while others remove distributional differences by construction.
>
> This means that changes in the appearance of PRS distributions should not automatically be interpreted as changes in genetic information or predictive performance.
>
> The objective is not to determine which population has higher genetic risk, but to quantify how robust our analytical conclusions are to technical decisions.
>
> [Insert GitHub repository link.]
>
> #StatisticalGenetics #Bioinformatics #PolygenicRiskScores

## 8. Strategic assessment

I consider this question stronger for scientific positioning than the original proposal.

The first version would demonstrate that you know how to compute a PRS and interpret distributions across populations.

The revised version lets you demonstrate something additional: the ability to critically assess the robustness of a computational procedure and identify which analytical decisions change its results.

This is relevant to three professional profiles:

| Career | Competency demonstrated |
|---|---|
| Academic faculty | Formulating methodological hypotheses, experimental design, and sensitivity analysis |
| Senior Bioinformatician | Pipeline validation, QC, harmonization, reproducibility, and robustness assessment |
| Computational Genomics / Industry | Characterizing sources of variability, standardization, and analytical reliability |

### Recommendation: an additional refinement

I would refine the scientific question to:

> How much do technical processing choices contribute to the observed differences in PRS distributions across populations, and which differences are introduced or removed by normalization?

This version separates two questions that are frequently conflated:

1. **Genomic processing:** can change the computed score values and, consequently, the observed differences between populations.
2. **Normalization:** can change the scale, the reference frame, and the interpretation of the scores without necessarily changing the underlying genetic information.

The first requires empirical evaluation.

The second combines known mathematical properties with an empirical evaluation of their consequences.

The most interesting outcome will not necessarily be finding large differences between pipelines. Demonstrating that certain conclusions remain stable despite the technical decisions is also a relevant methodological result.

That should be the project's goal: quantify the robustness of contrasts between PRS distributions, without conflating analytical sensitivity with predictive validity or genetic risk differences.
