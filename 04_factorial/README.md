# Análise fatorial — Etapa B

Desenho 2 (harmonização) × 7 (normalização), por par de população (10 pares).
O eixo "genótipos ausentes" ficou fora: só tem 1 nível observável nesses dados
(ver `02_harmonization_configs/config.yaml`), não dá pra fatorializar um eixo
sem variação.

## Pergunta

O quanto uma normalização reduz a separação entre populações depende de qual
configuração de harmonização foi usada antes? Testado como interação:

```
efeito_normalização(config) = D_g,h[normalização] − D_g,h[SUM bruto]  (dentro da mesma config)
interação = efeito_normalização(Config 2) − efeito_normalização(Config 1)
```

## Resultado

| normalização | \|interação\| média | \|interação\| máxima |
|---|---|---|
| `Z_global` | 0,0000 | 0,0000 |
| `Z_EUR` | 0,0000 | 0,0000 |
| `Z_percop` | 0,1560 | 0,3361 |
| `Z_MostSimilarPop` | 0,1564 | 0,3361 |
| `Z_norm1` | 0,1582 | 0,3491 |
| `Z_norm2` | 0,1583 | 0,3491 |

**As transformações afins simples (`Z_global`, `Z_EUR`) têm interação exatamente
zero** — confirma de novo, agora sob um ângulo diferente, a propriedade
matemática da Seção 4: como são a mesma transformação linear aplicada a todo
mundo, o "quanto elas mudam D_g,h" não pode depender de nada além da própria
escala, então é sempre zero, em qualquer configuração de harmonização.

**As normalizações por ancestralidade (nativas do pgsc_calc, e `Z_percop`)
têm interação real.** O caso mais extremo é o par EUR-SAS: no baseline a
separação já era quase nula (D=0,008), então a normalização tem pouco a
"achatar"; na Config 2 (ambíguas incluídas) surge separação real (D=-0,218),
e aí a normalização faz trabalho de verdade. O mesmo mecanismo aparece em
AFR-EAS e AMR-EAS.

## Interpretação para o projeto

Isso refina a conclusão da Seção 8 do `PROJECT_PLAN.md`: normalização por
ancestralidade não é uma correção universal e independente da harmonização —
**quanto ela "esconde" de uma escolha técnica anterior depende dessa escolha**.
Duas pessoas que discordam só na inclusão de 30 SNPs ambíguos podem chegar a
conclusões diferentes sobre "quanto a normalização resolveu o problema",
mesmo usando o mesmo método de normalização.

## Arquivos

- `factorial_Dgh.tsv` — grade completa 2×7, D_g,h por par de população.
- `factorial_interaction.tsv` — interação linha a linha (10 pares × 6 métricas, SUM excluído por ser a referência).
- `interaction_summary.tsv` — tabela acima.
