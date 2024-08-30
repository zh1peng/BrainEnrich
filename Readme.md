[![R-4.4-build](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-4.4-build.yml/badge.svg)](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-4.4-build.yml)
[![Codecov test coverage](https://codecov.io/gh/zh1peng/BrainEnrich/graph/badge.svg)](https://app.codecov.io/gh/zh1peng/BrainEnrich)
[![R-CMD-check](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/zh1peng/BrainEnrich/actions/workflows/R-CMD-check.yml)

![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/zh1peng/BrainEnrich/total)


# BrainEnrich: Revealing Biological Insights from Imaging-Derived Features through Transcriptomic Enrichment 🧠🧬

## Aim of the Toolbox 🎯

`BrainEnrich` is an R package designed to facilitate the correlation of imaging-derived phenotypes with transcriptional profiles. This toolbox aims to provide researchers and clinicians with robust statistical tools to uncover molecular architectures associated with cognitive functions, brain development, and disorders.

## Timeline of Development 🗓️

- **Q4 2023**: Initial conceptualization and development of core functions. 🛠️
- **Q1 2024**: Implementation of competitive null models and self-contained null models. 🧪
- **Q1 2024**: Testing with simulated datasets and refinement of statistical tests. 🔬
- **Q2 2024**: Beta release for community feedback and additional testing. 🔄
- **Q3 2024**: Incorporation of feedback and preparation for CRAN submission. ✍️
- **Q4 2024**: Submission to CRAN and publication of accompanying paper. 📰

## Installation 💾

*Please note that `brainEnrich` is currently in development and not yet available for installation.*

Once available, it can be installed from GitHub via the `devtools` package:

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

## Usage 📖

*Instructions on how to use the toolbox will be provided here, including example code.*

## To-Do List 📋
- [x] Initialize the project 2023/11/04
- [x] Finalize manuscript revision. 🔧
- [x] Development of core functions🔧
- [x] Create detailed vignettes for each major function. 📚
- [x] Optimize performance for large datasets. ⚡
- [ ] Conduct extensive testing with real-world data. 🌏
- [ ] Develop a comprehensive test suite. ✅
- [x] Set up continuous integration for automated testing. 🔄
- [ ] Prepare documentation for public release. 📄


## Contributing 🤝


## Versioning 🏷️
We use git for versioning. For the versions available, see the [tags on this repository](https://github.com/zh1peng/brainEnrich/tags).



## Authors 👩‍💻👨‍💻

* **Zhipeng Cao @ Xuhui Mental Health Center, Shanghai** - *Initial work* - [zh1peng](https://github.com/zh1peng)

## License 📜

This project is licensed under the GNU Affero General Public (AGP) License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments 👏
