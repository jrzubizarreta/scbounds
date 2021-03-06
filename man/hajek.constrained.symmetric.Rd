% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hajek.constrained.R
\name{hajek.constrained.symmetric}
\alias{hajek.constrained.symmetric}
\title{Computational workhorse function. Maximizes the value of hat{mu} among all Hajek
ratio estimators that satisfy the following properties:
- The ratio of the sampling weights gamma is bounded
- The resulting CDF is symmetric to within tolerance delta.}
\usage{
hajek.constrained.symmetric(Xcdf, xvals, sampling.ratio, center, delta)
}
\arguments{
\item{Xcdf}{The raw, unweighted, empirical CDF.}

\item{xvals}{The points at which Xcdf is specified.}

\item{sampling.ratio}{The bound on the sampling ratio.}

\item{center}{The center of symmetry.}

\item{delta}{The tolerance parameter.}
}
\value{
Xcdf.weighted A weighted version of Xcdf that maximizes hat{mu}
}
\description{
Computational workhorse function. Maximizes the value of hat{mu} among all Hajek
ratio estimators that satisfy the following properties:
- The ratio of the sampling weights gamma is bounded
- The resulting CDF is symmetric to within tolerance delta.
}

