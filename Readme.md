[![R-CMD-check](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/zh1peng/BrainEnrich/graph/badge.svg)](https://app.codecov.io/gh/zh1peng/BrainEnrich)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/zh1peng/BrainEnrich/total?logo=plume&label=Download&color=%23328da8)

<p align="left">
  <img src="images/sticker.png" alt="BrainEnrich sticker" align="right" height="160" style="margin-right: 15px;">
  <span style="vertical-align: top; font-size: 16px; text-align: justify;">
    <strong>BrainEnrich</strong> is an R package for transcriptomic enrichment analysis of brain imaging phenotypes. It supports group-level enrichment, individual-level gene-set scoring, null-model based inference, and simulation workflows with native BrainEnrich result objects and plotting utilities.
  </span>
</p>

<p align="center">
  <img src="images/workflow1.png" alt="BrainEnrich workflow">
  <br>
  <em>Overview of the BrainEnrich analysis workflow.</em>
</p>

## Features

- Group-level enrichment with `brainenrich()`
- Individual-level scoring with `brainscore()`
- Linear-model testing with `brainscore.lm_test()`
- Simulation workflows with `brainscore.simulate()`
- Built-in gene-set resources including SynGO, GO, KEGG, Reactome, WikiPathways, MeSH, DisGeNET, and cell-type collections
- Native result plotting via `plot_terms()`, `plot_core_genes()`, `plot_heatmap_terms()`, and `plot_term_network()`

## Native outputs

- `brainenrich()` and `brainscore.lm_test()` return `EnrichRes`
- `get_annoData()` returns `EnrichAnno`
- `brainscore()` and `brainscore.simulate()` keep native score/simulation list outputs

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("zh1peng/BrainEnrichData")
remotes::install_github("zh1peng/BrainEnrich")
```

Optional vignette and data-preparation workflows may require extra packages such as `ggplot2`, `ggseg`, `kableExtra`, or Bioconductor annotation packages, but BrainEnrich no longer depends on external enrichment or plotting packages at runtime.

Legacy projects can stay on the pre-refactor DOSE-based line with:

```r
remotes::install_github("zh1peng/BrainEnrich@v1.0.0")
```

## Quick start

```r
library(BrainEnrich)

data(brain_data)
data(perm_id_dk_lh_5000)

gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
annoData <- get_annoData(type = "SynGO")

res <- brainenrich(
  brain_data = brain_data,
  gene_data = gene_data,
  annoData = annoData,
  perm_id = perm_id_dk_lh_5000,
  cor_method = "pearson",
  aggre_method = "mean",
  null_model = "spin_brain",
  n_perm = 10,
  n_cores = 1,
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 1
)

plot_terms(res, type = "dot")
```

## Documentation

Tutorials are available at [github-pages](https://zh1peng.github.io/BrainEnrich/).

- Group-level enrichment: `vignettes/brainenrich_demo.Rmd`
- Individual-level enrichment: `vignettes/brainscore_demo.Rmd`
- Simulation analysis: `vignettes/brainscore.simulate_demo.Rmd`
- Gene-set resources: `vignettes/annoData_demo.Rmd`
- HPC workflows: `vignettes/hpc_demo.Rmd`

## Notes on data resources

BrainEnrich now keeps the full curated expression and gene-set resources in the companion package `BrainEnrichData`. The loaders prefer:

1. An installed `BrainEnrichData` package
2. A sibling development checkout at `../BrainEnrichData`
3. Local BrainEnrich package resources, when present
4. A user cache download from the matching `BrainEnrichData` release

## License

This project is licensed under the GNU Affero General Public License. See [LICENSE.md](LICENSE.md) for details.
