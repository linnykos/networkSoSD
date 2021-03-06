% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/greedy_refinement.R
\name{greedy_refinement}
\alias{greedy_refinement}
\title{Greedy clustering}
\usage{
greedy_refinement(adj_list, K, init_clust = NA, iter_max = 100, verbose = F)
}
\arguments{
\item{adj_list}{list of adjacency matrices}

\item{K}{positive integer}

\item{init_clust}{initial membership vector. If \code{NA}, we matrixify the tensor and perform K-means}

\item{iter_max}{positive integer}

\item{verbose}{boolean}
}
\value{
membership vector
}
\description{
Greedy clustering
}
