% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symmetric.R
\name{bounds.symmetric.internal}
\alias{bounds.symmetric.internal}
\title{Computes upper identification interval with symmetry constraints.}
\usage{
bounds.symmetric.internal(X, sampling.ratio = 5, xmin = NULL, xmax = NULL,
  buckets = 1000, alpha = 1/sqrt(length(X)))
}
\arguments{
\item{X}{The observed data.}

\item{sampling.ratio}{Bound on the sampling weights gamma.}

\item{xmin}{Used to construct histogram representation.}

\item{xmax}{Used to construct histogram representation.}

\item{buckets}{Used to construct histogram representation.}

\item{alpha}{Significance level used for KS bounds.}
}
\value{
mu.bound The upper bound for mu(x).

Xhat Unweighted empirical CDF of the data.

xvals Points at which Xhat is evaluated.

Xhat.weighted Weighted version of Xhat that maximizes mu, subject to symmetry.
}
\description{
Computes upper identification interval with symmetry constraints.
}

