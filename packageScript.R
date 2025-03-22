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

# check if devtools is ready
devtools::has_devel()


#
devtools::document()
#devtools::check()
devtools::check(args="--no-tests --no-vignettes")
devtools::check(args="--no-tests --no-vignettes --as-cran")
styler::style_pkg()
devtools::build()
devtools::install()



# do coverage test
install.packages("covr")
# create test files for functions
usethis::use_test("get_annoData")
usethis::use_test("get_geneSetList")
usethis::use_test("split_Anno")
usethis::use_test("filter_geneSetList")
usethis::use_test("get_termDescription")
usethis::use_test("get_geneExp")
usethis::use_test("brainscore.simulate.R")
usethis::use_test("brainenrich")
usethis::use_test("rotate_parcellation")
usethis::use_test("job_splitter")
usethis::use_test("job_cat")
usethis::use_test("brainscore")
usethis::use_test("misc")
usethis::use_test("brainscore.lm_test.R")
usethis::use_test("plot_functions")
usethis::use_test("corr_brain_gene")

# run coverage test and upload to codecov
covr::codecov(token = "bf94b382-482f-4c28-9ced-e988216dde4a")

# Report coverage for a specific function (e.g., `job_splitter`)
coverage <- covr::function_coverage(
  fun = job_cat,               # The function to test
  test_file("tests/testthat/test-job_cat.R")  # The specific test file(s)
)
# test a specific function
test_file("tests/testthat/test-misc.R")

# build vignettes
devtools::build_vignettes()
# install pdflatex for building manual
devtools::build_manual()



usethis::use_citation()

# options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# install.packages("pkgdown")
#usethis::use_pkgdown_github_pages()
usethis::use_pkgdown()
pkgdown::build_site(lazy=TRUE)
pkgdown::clean_site()
pkgdown::build_articles()
pkgdown::deploy_to_branch()
pkgdown::build_home()
pkgdown::build_home_index()


usethis::use_github_action('check-standard')
usethis::use_github_action('pr-commands')
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


library(hexSticker)
library(ggplot2)

p <- ggplot(aes(x = mpg, y = wt), data = mtcars) + geom_point(fill="white", color="white") 
p <- p + theme_void() + theme_transparent()
sticker(p, package="BrainEnrich", p_size=20, p_y=0.52, p_color="#000000", 
        s_x=1, s_y=.75, s_width=1.3, s_height=1,
        h_fill="#025a63", h_color="#000000", 
        url="https://zh1peng.github.io/BrainEnrich/",u_size=3, u_color="#000000",
        filename="inst/figures/ggplot2.png")


imgurl <- "C:/Users/Zhipeng/Desktop/tmp2add2021/test1.png"

sticker(subplot=imgurl, asp=0.4,  package="BrainEnrich", p_size=20, p_y=0.52, p_color="#ffffff", 
        s_x=1, s_y=1, s_width=1.6, s_height=1.2,
        h_fill="#003c61", h_color="#000000", 
        url="https://zh1peng.github.io/BrainEnrich/",u_size=3, u_color="#ffffff",
        filename="inst/figures/url.png")