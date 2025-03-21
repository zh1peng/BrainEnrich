# .onLoad <- function(libname, pkgname) {
#   ensure_bioc_packages()
# }


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ", pkgname, "!")
  # packageStartupMessage(
  #   "To cite the package, please use the following reference:\n",
  #   "1. Cao, Z., Zhan, G., Qin, J., Cupertino, R. B., Ottino-Gonzalez, J., Murphy, A., ... & Garavan, H. (2024). ",
  #   "Unraveling the molecular relevance of brain phenotypes: A comparative analysis of null models and test statistics. Neuroimage, 293, 120622.\n",
  #   "2. Cao Z (2024). BrainEnrich: Revealing Biological Insights from Imaging-Derived Features through Transcriptomic Enrichment. ",
  #   "R package version 1.0, <https://zh1peng.github.io/BrainEnrich/>."
  # )
  # Uncomment the following lines if you want to dynamically fetch the citation
  # citation_info <- citation(pkgname)
  # packageStartupMessage(capture.output(print(citation_info, style = "text")))
}


# # Ensure Bioconductor Dependencies are Installed
# #' @importFrom utils install.packages
# ensure_bioc_packages <- function() {
#   if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
#   }
#   bioc_packages <- c("DOSE")
#   for (pkg in bioc_packages) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       BiocManager::install(pkg)
#     }
#   }
# }
