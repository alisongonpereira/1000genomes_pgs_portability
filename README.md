# PRS Portability Audit

Quanto das diferenças observadas nas distribuições de Polygenic Risk Scores entre
populações pode ser atribuído às escolhas técnicas de processamento e normalização?

*To what extent do technical processing and normalization choices alter the
observed differences in polygenic risk score distributions across populations?*

Não avaliamos a capacidade do PRS de predizer a doença. A variável de interesse é
a **estabilidade da distribuição do score** diante de escolhas analíticas
controladas — harmonização alélica e normalização — mantendo fixos os dados e o
modelo genético. Plano científico completo em [`docs/PROJECT_PLAN.md`](docs/PROJECT_PLAN.md).

## Dados

| Componente | Escolha |
|---|---|
| PRS | PRS313_BC — [PGS000004](https://www.pgscatalog.org/score/PGS000004/) (Mavaddat et al. 2019), câncer de mama |
| Genoma | GRCh37, harmonizado pelo PGS Catalog |
| Genótipos | 1000 Genomes Phase 3 (2504 amostras, 5 superpopulações: AFR/AMR/EAS/EUR/SAS) |
| Pipeline | [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) v2.3.0 (nf-core/Nextflow), perfil conda |

## Estrutura do repositório

```
01_baseline_prs/        baseline: harmonização conservadora (SNPs ambíguos excluídos),
                         SUM bruto + normalizações nativas do pgsc_calc (Z_norm1/Z_norm2/
                         Z_MostSimilarPop via PCA+RF contra painel de ancestralidade dedicado)
02_harmonization_configs/ Configuração 2 (SNPs ambíguos incluídos) + validação da suposição
                         de mesma fita por frequência alélica (EUR, 1000G vs. GWAS original)
03_figures/              figura de 3 painéis (distribuições / sensibilidade / estabilidade individual)
04_factorial/            interação harmonização × normalização sobre a separação entre populações
scripts/                 scripts numerados que gera cada saída acima
docs/                    plano científico completo e material de divulgação
```

Dados brutos baixados e saídas intermediárias do pgsc_calc (`resources/`, `*/pgsc_calc_run/`,
`*/keep_ambiguous/`) não são versionados — são regenerados pelos scripts, não são fonte.

## Principais achados até agora

- **Transformações afins comuns preservam a separação entre populações, exatamente.**
  `Z_global` e `Z_EUR` (normalização por média/DP global ou só de EUR) reproduzem a
  diferença padronizada (Cohen's d) do `SUM` bruto até a 4ª casa decimal, em todos os
  10 pares de população — confirmação empírica de uma propriedade matemática, não uma
  descoberta do dado.
- **Normalização por ancestralidade reduz mas não zera a separação.** `Z_norm1`/`Z_norm2`
  (regressão contínua em PCs, nativas do pgsc_calc) e `Z_MostSimilarPop` (padronização
  discreta pela população mais similar) aproximam as médias de zero em todas as
  superpopulações, mas não de forma idêntica entre si.
- **A escolha de harmonização muda a interpretação, não só a escala.** Incluir as 30
  variantes ambíguas (A/T, C/G) do PGS313_BC — decisão validada empiricamente por
  frequência alélica (30/30 coerentes contra EUR do 1000G) — altera a separação entre
  pares de população de forma não trivial: o par EUR-SAS **inverte de sinal**
  (D = 0,008 → −0,218).
- **Correlação de ranking alta não implica estabilidade individual.** Entre as duas
  configurações de harmonização, Spearman ρ = 0,954, mas **51,9% das amostras mudam
  de decil**.
- **Normalização e harmonização interagem.** O quanto uma normalização por
  ancestralidade "esconde" uma diferença populacional depende de qual configuração de
  harmonização foi usada antes (interação média |Δ| ≈ 0,16, até 0,35) — as
  transformações afins simples (`Z_global`/`Z_EUR`) não têm essa interação (é zero,
  exatamente, pela mesma propriedade matemática acima).
- O eixo "genótipos ausentes" (exclusão vs. imputação) não é testável com esta fonte
  de dados: o VCF integrado do 1000G Phase 3 não tem nenhum genótipo ausente.

Detalhes, números completos e a discussão de cada achado estão nos `README`/`config.yaml`
de cada pasta numerada.

## Reprodução

Cada pasta numerada documenta, no seu `config.yaml`/`README.md`, exatamente qual comando
do pgsc_calc (ou script em `scripts/`) gerou aquela saída, incluindo os parâmetros que
divergem do default. Ambiente: Python 3 (pandas/matplotlib/seaborn/scipy) + Nextflow +
pgsc_calc via conda/micromamba.
