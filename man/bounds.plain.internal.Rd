% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plain.R
\name{bounds.plain.internal}
\alias{bounds.plain.internal}
\title{Computes upper identification interval without shape constraints, as
in Aronow & Lee (2013).}
\usage{
bounds.plain.internal(X, sampling.ratio = 5, xmin = NULL, xmax = NULL,
  buckets = 1000)
}
\arguments{
\item{X}{The observed data.}

\item{sampling.ratio}{Bound on the sampling weights gamma.}

\item{xmin}{Used to construct histogram representation.}

\item{xmax}{Used to construct histogram representation.}

\item{buckets}{Used to construct histogram representation.}
}
\value{
mu.bound The upper bound for mu(x).

Xhat Unweighted empirical CDF of the data.

xvals Points at which Xhat is evaluated.

Xhat.weighted Weighted version of Xhat that maximizes mu, subject to sampling ratio constraint.
}
\description{
Computes upper identification interval without shape constraints, as
in Aronow & Lee (2013).
}

