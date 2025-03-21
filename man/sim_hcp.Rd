% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sim_hcp}
\alias{sim_hcp}
\title{Simulated data from HCP data}
\format{
A data frame with rows as subjects and columns as age, sex, BMI and regional cortical thickness values
}
\source{
mvrnorm(n = 100, mu = mean_vals, Sigma = cov_matrix); sample(df.hcp$Age_in_Yrs, 100, replace = FALSE); sample(df.hcp$Sex, 100, replace = FALSE); sample(df.hcp$BMI, 100, replace = FALSE)
}
\usage{
sim_hcp
}
\description{
Simulated data from HCP data
}
\keyword{datasets}
