[![R-CMD-check](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/zh1peng/BrainEnrich/graph/badge.svg)](https://app.codecov.io/gh/zh1peng/BrainEnrich)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/zh1peng/BrainEnrich/total?logo=plume&label=Download&color=%23328da8)


# BrainEnrich: Revealing Biological Insights from Imaging-Derived Features through Transcriptomic Enrichment 🧠🧬

## Aim of the Toolbox 🎯
**BrainEnrich** is an R package designed to facilitate the integration of brain imaging data with transcriptomic profiles. It enables researchers to explore the molecular underpinnings of brain phenotypes by performing enrichment analysis of predefined gene sets. Whether working at the group or individual level, the package offers a flexible and powerful tool for examining associations between brain imaging phenotypes (e.g., cortical thickness) and gene expression, using a variety of statistical models, null models, and aggregation methods. 

![BrainEnrich](images/workflow.png)  
_Overview of the BrainEnrich package and its analysis workflow._

## Features 🚀

- **Group-Level Enrichment Analysis**: Correlate group-level imaging phenotypes (e.g., cortical thickness effect maps) with transcriptional profiles from the Allen Human Brain Atlas (AHBA).
- **Individual-Level Enrichment Analysis**: Apply enrichment analysis to individual imaging-derived phenotypes (IDPs) to explore personalized brain signatures.
- **Transcriptional Profiles**: Preprocessed gene expression data from the AHBA, Desikan-Killiany, and Schaefer atlases.
- **Predefined Gene Sets**: Includes Gene Ontology (GO), DisGeNet, KEGG, WikiPathways, Reactome, MeSH, SynGO, and Cell Type gene sets for flexible pathway analysis.
- **Multiple Association Methods**: Pearson and Spearman correlations, Partial Least Squares (PLS) regression, and user-defined methods for exploring gene-imaging associations.
- **Aggregation Methods**: Multiple options for aggregating gene set scores, including mean, median, and Kolmogorov-Smirnov (KS)-based statistics.
- **Null Models**: Includes both self-contained and competitive null models to assess the significance of gene set associations.
- **Core Genes Identification**: Identify core genes within gene sets using a leave-one-out procedure for more interpretable results.
- **Simulation Studies**: Type I error and power simulations to assess the reliability of the analysis methods.



## Installation 💾

```r
# Install remotes if you haven't already
if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("DOSE")
# Install brainEnrich from GitHub
remotes::install_github("zh1peng/BrainEnrich")
```


## Usage 🔬

### Group-Level Enrichment Analysis
The group-level analysis involves correlating brain imaging phenotypes (e.g., cortical thickness) with gene expression data. The results are summarized in **GS scores** representing the association between each predefined gene set and the imaging phenotype.

### Individual-Level Enrichment Analysis
For individual participants, **BrainEnrich** calculates **GS scores** for each participant, enabling the exploration of personalized molecular signatures and their relationships with phenotypic outcomes.

### Association Methods
You can use **Pearson's** or **Spearman's** correlation, or **Partial Least Squares (PLS) regression** to investigate gene-imaging associations. The package also supports user-defined methods for advanced analyses.

### Predefined Gene Sets
The package includes predefined gene sets from various sources like **Gene Ontology**, **DisGeNet**, **KEGG**, **WikiPathways**, **Reactome**, **MeSH**, **SynGO**, and **Cell Type** gene sets. These can be used for enrichment analysis with a variety of biological contexts.

### Statistical Tests & Null Models
The package provides both **self-contained** and **competitive null models** to assess statistical significance and control false positives. It also includes co-expression-matched null models for more advanced analysis.

### Simulations
**BrainEnrich** includes simulations to evaluate **Type I error rates** and **statistical power** for enrichment analyses, helping researchers ensure their results are reliable.


## Contributing 🤝
We welcome contributions to **BrainEnrich**! Feel free to fork the repository, and create pull requests.

## Versioning 🏷️
We use git for versioning. For the versions available, see the [tags on this repository](https://github.com/zh1peng/brainEnrich/tags).



## Authors 👩‍💻👨‍💻

* **Zhipeng Cao @ Xuhui Mental Health Center, Shanghai** - *Initial work* - [zh1peng](https://github.com/zh1peng)

## License 📜

This project is licensed under the GNU Affero General Public (AGP) License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments 👏

Thanks to my family for their support.