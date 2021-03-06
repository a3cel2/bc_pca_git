% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_connectivity_permutation.R
\name{connectivity_histogram}
\alias{connectivity_histogram}
\title{A histogram of edgewise-shuffled largest component size iterations}
\usage{
connectivity_histogram(pca_universe, pca_enhanced, pca_depleted, condition,
  iterations = 10000, metric = get_largest_component_from_edgelist,
  seed = 123, text_size = 1.5, xlim = c(0, 25), breakstep = 1,
  xlab = "Simulated largest component", sampling_mode = "edgewise",
  prob_node = 0)
}
\arguments{
\item{pca_universe}{A PCA file corresponding to all protein-protein interactions
available in a given condition}

\item{pca_enhanced}{A PCA file corresponding to all significantly enhanced
protein-protein interactions}

\item{pca_depleted}{A PCA file corresponding to all significantly depleted
protein-protein interactions}

\item{condition}{condition, a string}

\item{iterations}{number of iterations to plot in the histogram}

\item{metric}{the metric to be tested, defaults to largest connected component size}

\item{seed}{seed for the random number generator}
}
\value{
plots a histogram
}
\description{
A histogram of edgewise-shuffled largest component size iterations
}

