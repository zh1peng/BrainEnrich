% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_core_genes.R
\name{find_core_genes}
\alias{find_core_genes}
\title{Find Core Genes Influencing Aggregated Score or LM Coefficients between molecular profile and behavioral data}
\usage{
find_core_genes(
  geneList,
  geneSetList,
  pred_df = NULL,
  cov_df = NULL,
  aggre_method,
  n_cores = 1,
  threshold_type = c("sd", "percentile"),
  threshold_value = 1
)
}
\arguments{
\item{geneList}{A matrix of genes by subs, each column representing a subject / a group-level result.}

\item{geneSetList}{A list of gene sets, each containing names of genes.}

\item{pred_df}{Optional data frame of a predictor. If NULL, it is perfomred for group-level enrichment.}

\item{cov_df}{Optional data frame of covariates. If NULL, it is perfomred for group-level enrichment.}

\item{aggre_method}{The aggregation method used to compute the scores.}

\item{n_cores}{The number of cores to use for parallel processing; defaults to 1.
Uses all available cores minus one if set to 0.}

\item{threshold_type}{The method to determine significance ('sd' for standard deviation, 'percentile' for percentile threshold).}

\item{threshold_value}{Numeric value specifying the threshold level; meaning depends on \code{threshold_type}.}
}
\value{
A list of core genes for each gene set.
}
\description{
This function performs a Leave-One-Out (LOO) analysis on gene sets to determine core genes
that influence the aggregated score. It can utilize parallel processing
to enhance computation efficiency and supports two types of analysis: one that considers
only gene sets and another that includes predictor and covariate data frames.
}
