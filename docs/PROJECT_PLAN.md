# PRS Portability Audit — Projeto Revisado

## Pergunta científica central

Quanto das diferenças observadas nas distribuições de Polygenic Risk Scores entre populações pode ser atribuído às escolhas técnicas de processamento e normalização?

**Em inglês:**
To what extent do technical processing and normalization choices alter the observed differences in polygenic risk score distributions across populations?

Essa reformulação torna o projeto metodologicamente mais interessante. Em vez de simplesmente mostrar que populações apresentam distribuições distintas de PRS, vamos investigar quanto essas diferenças mudam quando modificamos decisões analíticas específicas, mantendo constantes os dados e o modelo genético.

Há, porém, uma distinção fundamental: podemos quantificar a sensibilidade das distribuições às escolhas técnicas, mas não necessariamente decompor a diferença observada em uma porcentagem exclusivamente biológica e outra exclusivamente técnica.

Essa distinção será central para a credibilidade do projeto.

## 1. Objetivo científico

### Objetivo geral

Quantificar a influência de diferentes estratégias de processamento genômico e normalização sobre as distribuições de um mesmo PRS em populações geneticamente diversas.

### Objetivos específicos

1. Avaliar como decisões de harmonização alélica, inclusão de variantes e tratamento de genótipos ausentes modificam os scores individuais.
2. Comparar os efeitos da normalização global, da normalização baseada em uma população de referência e da normalização específica por população.
3. Quantificar a estabilidade das diferenças entre populações sob diferentes configurações analíticas.
4. Determinar em quais circunstâncias a normalização altera a interpretação dos resultados sem modificar a informação genética utilizada.

**Hipótese central:** parte da variação observada nas distribuições de PRS depende das decisões analíticas, mas o efeito dessas decisões não é uniforme nem necessariamente equivalente entre populações.

## 2. Dados e escopo

Manteremos os dados propostos originalmente.

| Componente | Escolha |
|---|---|
| PRS | PRS313_BC — PGS000004 |
| Fenótipo associado | Câncer de mama |
| Fonte do PRS | PGS Catalog |
| Dados genômicos | 1000 Genomes Project — Phase 3 |
| Populações | AFR, AMR, EAS, EUR e SAS |
| Linguagem | Python |
| Análise | Comparação de configurações analíticas |
| Produto | Repositório reproduzível + figura científica |

O experimento não avaliará a capacidade do PRS de predizer câncer de mama.

A variável de interesse será a estabilidade das distribuições dos scores diante de alterações controladas no processamento.

## 3. Desenho experimental

Proponho organizar o estudo em três etapas.

### Etapa A — Construir uma distribuição de referência

Primeiro, calcularemos o PRS utilizando uma configuração analítica previamente definida.

Essa configuração será o nosso baseline.

As condições serão:

- Mesmo arquivo de pesos.
- Mesmo conjunto de indivíduos.
- Mesmo genoma de referência.
- Harmonização alélica validada.
- Critérios explícitos de inclusão e exclusão de variantes.
- Registro de variantes ausentes e genótipos não chamados.

O resultado será uma distribuição inicial dos scores para cada população.

Essa etapa é importante porque não podemos comparar o efeito de diferentes pipelines sem definir uma referência metodológica.

### Etapa B — Modificar decisões técnicas isoladamente

Aqui está o núcleo do projeto.

**Matriz experimental**

| Dimensão | Configuração 1 | Configuração 2 |
|---|---|---|
| Harmonização | Apenas variantes não ambíguas | Inclusão adicional de variantes resolvidas por critérios validados |
| Genótipos ausentes | Exclusão de loci incompletos | Imputação por dosagem esperada, com referência explicitada |
| Normalização | Média e DP globais | Média e DP de EUR |
| Normalização adicional | — | Média e DP específicos por população |

A harmonização deve preservar a orientação correta dos alelos. Não devemos comparar um pipeline correto com outro deliberadamente incorreto como se fossem duas alternativas igualmente válidas.

A normalização específica por população será tratada separadamente, como uma transformação diagnóstica.

**Por quê?**

Porque normalizar individualmente cada população para média zero e desvio-padrão um elimina, por construção, as diferenças entre médias e variâncias.

Isso não demonstra que as diferenças anteriores eram artefatos.

Demonstra apenas que a transformação removeu essas características das distribuições.

Essa distinção será um dos resultados conceituais mais importantes do projeto.

### Etapa C

*(Reservado para a análise fatorial das combinações relevantes — ver seção 5.)*

## 4. Como responder quantitativamente à pergunta?

Precisamos definir o que significa uma distribuição ter mudado.

Eu utilizaria quatro indicadores.

| Indicador | O que quantifica | Interpretação |
|---|---|---|
| Diferença padronizada entre médias | Separação entre populações | Sensibilidade dos contrastes populacionais |
| Distância de Wasserstein | Distância entre distribuições | Mudanças na distribuição como um todo |
| Correlação de Spearman | Estabilidade da ordenação individual | Quanto as decisões alteram o ranking dos indivíduos |
| Mudança de percentil | Reclassificação relativa | Quantos indivíduos mudam de posição conforme o processamento |

As distâncias entre distribuições devem ser apresentadas em uma escala comparável entre configurações.

### Uma propriedade matemática importante

Considere um score bruto $S_i$.

Na normalização global:

$$Z_i=\frac{S_i-\mu_{\mathrm{global}}}{\sigma_{\mathrm{global}}}$$

Na normalização baseada na população EUR:

$$Z_{i,\mathrm{EUR}}=\frac{S_i-\mu_{\mathrm{EUR}}}{\sigma_{\mathrm{EUR}}}$$

Ambas são transformações afins comuns a todos os indivíduos.

Consequentemente, quando aplicadas ao mesmo conjunto de scores, preservam a ordenação individual e as diferenças padronizadas entre populações.

Portanto, se utilizarmos a diferença padronizada entre médias como indicador, essas duas normalizações produzirão o mesmo resultado.

Isso é uma propriedade matemática, não uma descoberta empírica.

Por outro lado, a normalização específica por população:

$$Z_{i,g}=\frac{S_i-\mu_g}{\sigma_g}$$

modifica a posição relativa dos indivíduos quando comparamos populações distintas.

Essa transformação altera o referencial da comparação.

O projeto precisa separar mudanças reais decorrentes do processamento genômico de mudanças esperadas matematicamente pela escolha da escala.

## 5. Como quantificar a contribuição das decisões técnicas?

Proponho comparar cada configuração com o baseline.

Para duas populações $g$ e $h$, definimos:

$$D_{g,h}^{(k)}$$

como a distância padronizada entre suas distribuições sob a configuração analítica $k$.

A sensibilidade à configuração pode ser expressa por:

$$\Delta D_{g,h}^{(k)} = D_{g,h}^{(k)} - D_{g,h}^{(\mathrm{baseline})}$$

Assim, poderemos responder perguntas concretas:

- Quanto a inclusão de variantes adicionais, corretamente harmonizadas, modificou a separação entre AFR e EUR?
- Quanto o tratamento de genótipos ausentes modificou a distribuição dos scores em AMR?
- Qual proporção dos indivíduos mudou de decil quando alteramos a estratégia de processamento?

Para avaliar a incerteza dessas estimativas, utilizaremos bootstrap com reamostragem estratificada por população.

**Ressalva:** se modificarmos vários parâmetros simultaneamente, suas contribuições poderão interagir. Por isso, a análise principal deverá alterar um fator por vez, seguida de uma análise fatorial das combinações que apresentarem efeitos relevantes.

### O que poderemos concluir?

Uma conclusão possível, caso os resultados a sustentem, seria:

> Technical processing choices accounted for substantial variation in the observed separation between PRS distributions.

Mas eu evitaria:

> Technical artifacts explain X% of ancestry-related PRS differences.

A segunda afirmação exigiria uma decomposição causal que esse desenho, isoladamente, não permite estabelecer.

## 6. Figura científica principal

Minha proposta é uma figura composta por três painéis.

### Painel A — Distribuições dos PRS

Violin plots ou density plots mostrando as distribuições por população em diferentes configurações.

A comparação deverá utilizar escalas adequadas para evitar que uma mudança de unidade seja interpretada como mudança substantiva.

### Painel B — Sensibilidade dos contrastes populacionais

Heatmap com:

- Linhas: comparações entre populações.
- Colunas: configurações analíticas.
- Valores: mudanças nas diferenças padronizadas em relação ao baseline.

Esse será o painel principal para responder à pergunta científica.

### Painel C — Estabilidade individual

Scatter plot comparando os scores individuais entre duas configurações de processamento.

Complementaremos com a correlação de Spearman e a proporção de indivíduos que mudaram de decil.

Esse painel mostrará se as escolhas técnicas apenas modificaram a escala ou também alteraram a ordenação dos indivíduos.

## 7. Publicação em inglês — versão revisada

**Título sugerido:** How much of the variation in PRS distributions is driven by analytical choices?

**Rascunho condicionado à execução do projeto:**

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

## 8. Avaliação estratégica

Considero essa pergunta mais forte para posicionamento científico do que a proposta original.

A primeira versão demonstraria que você sabe calcular PRS e interpretar distribuições entre populações.

A versão revisada permite demonstrar algo adicional: capacidade de avaliar criticamente a robustez de um procedimento computacional e identificar quais decisões analíticas modificam seus resultados.

Isso é relevante para três perfis profissionais:

| Carreira | Competência demonstrada |
|---|---|
| Academic faculty | Formulação de hipóteses metodológicas, desenho experimental e análise de sensibilidade |
| Senior Bioinformatician | Validação de pipelines, QC, harmonização, reprodutibilidade e avaliação de robustez |
| Computational Genomics / Industry | Caracterização de fontes de variabilidade, padronização e confiabilidade analítica |

### Recomendação: uma delimitação adicional

Eu refinaria a pergunta científica para:

> How much do technical processing choices contribute to the observed differences in PRS distributions across populations, and which differences are introduced or removed by normalization?

Essa versão separa duas questões que frequentemente são confundidas:

1. **Processamento genômico:** pode modificar os valores calculados dos scores e, consequentemente, as diferenças observadas entre populações.
2. **Normalização:** pode modificar a escala, o referencial e a interpretação dos scores sem necessariamente modificar a informação genética subjacente.

A primeira exige avaliação empírica.

A segunda combina propriedades matemáticas conhecidas com uma avaliação empírica de suas consequências.

O resultado mais interessante não será necessariamente encontrar grandes diferenças entre pipelines. Demonstrar que determinadas conclusões permanecem estáveis apesar das decisões técnicas também constitui um resultado metodológico relevante.

Esse deve ser o objetivo do projeto: quantificar a robustez dos contrastes entre distribuições de PRS, sem confundir sensibilidade analítica com validade preditiva ou diferenças de risco genético.
