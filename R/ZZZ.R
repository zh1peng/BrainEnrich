.onLoad <- function(libname, pkgname) {
  ensure_bioc_packages()
}


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ", pkgname, "!")
  packageStartupMessage("To cite the package, please use the following reference:")
  citation_info <- citation(pkgname)
  packageStartupMessage(capture.output(print(citation_info, style = "text")))
}


# Ensure Bioconductor Dependencies are Installed
ensure_bioc_packages <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  bioc_packages <- c("DOSE")
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }
}
