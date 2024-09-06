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
#devtools::check()
devtools::check(run_tests = FALSE)
styler::style_pkg()
devtools::build()
devtools::install()


devtools::test()
usethis::use_testthat()
usethis::use_test("get_annoData")
usethis::use_test("get_geneSetList")
usethis::use_test("split_Anno")
usethis::use_test("filter_geneSetList")
usethis::use_test("get_termDescription")
usethis::use_test("get_geneExp")
usethis::use_test("brainscore.simulate.R")

usethis::use_test("brainenrich")
usethis::use_test("rotate_parcellation")




install.packages("covr")
covr::package_coverage()
covr::codecov(token = "bf94b382-482f-4c28-9ced-e988216dde4a")



detach("package:BrainEnrich", unload = TRUE)
library(BrainEnrich)

devtools::build_vignettes()
^vignettes$
^.*\.Rmd$
^.*\.Rnw$


# install pdflatex for building manual
devtools::build_manual()


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DOSE", "clusterProfiler", "org.Hs.eg.db"))


# action on github uses pak to install the package
# install.packages("pak")
# pak::pkg_install("BrainEnrich")



usethis::use_mit_license()
usethis::use_citation()
usethis::use_agpl3_license()


usethis::use_pkgdown()
usethis::edit_r_environ()
usethis::browse_github_pat()

usethis::use_pkgdown_github_pages()

usethis::use_pkgdown()
pkgdown::build_site()
pkgdown::deploy_to_branch()


usethis::use_github_action('check-standard')
# prepare example data



install.packages("styler")

styler::style_pkg()

devtools::build_vignettes()

usethis::use_github_action("lint")


usethis::gh_token_help()
usethis::use_github_action()
usethis::use_github_action('test-coverage')
usethis::use_pkgdown_github_pages()

devtools::install_github("zh1peng/BrainEnrich")


# fix push issue
git config --global --unset http.proxy
git config --global http.proxy http://127.0.0.1:7890
git push origin

R CMD Rd2pdf .