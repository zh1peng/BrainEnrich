# install necessary packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}

if (!requireNamespace("roxygen2", quietly = TRUE)) {
  install.packages("roxygen2")
}

if (!requireNamespace("testthat", quietly = TRUE)) {
  install.packages("testthat")
}

devtools::has_devel()

devtools::document()
devtools::build()
devtools::install()

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DOSE", "clusterProfiler", "org.Hs.eg.db"))


# action on github uses pak to install the package
# install.packages("pak")
# pak::pkg_install("BrainEnrich")



usethis::use_mit_license()
usethis::use_citation()
usethis::use_agpl3_license()

uesthis::check_travis() # Check travis

usethis::use_data()

usethis::use_testthat() # Create testthat infrastructure
usethis::use_test('corr_brain_gene') # Create test file for corr_brain_gene function
devtools::test(filter = "corr_brain_gene") # Run tests for corr_brain_gene function
# testthat::test_file("tests/testthat/test-corr_brain_gene.R") # Run tests for corr_brain_gene function

usethis::use_test('aggregate_geneSet') # Create test file for aggregate_gene_data function
devtools::test(filter = "aggregate_geneSet") 

usethis::use_test('get_geneSets.R')
usethis::use_test('brainenrich.R')
usethis::use_test('brainscore.R')

usethis::use_vignette('brainenrich')
usethis::use_vignette('visualization')
usethis::use_vignette('brainscore')

usethis::use_pkgdown()
usethis::edit_r_environ()
usethis::browse_github_pat()

usethis::use_pkgdown()
pkgdown::build_site()
pkgdown::deploy_to_branch()
usethis::use_pkgdown_github_pages()

usethis::use_github_action()
# prepare example data



install.packages("styler")

styler::style_pkg()


usethis::use_github_action("lint")

usethis::use_github_action()

devtools::install_github("zh1peng/BrainEnrich")