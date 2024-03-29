% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indices.R
\name{maybe_dither_lx}
\alias{maybe_dither_lx}
\title{maybe_dither_lx}
\usage{
maybe_dither_lx(lx, pert = 1e-04)
}
\arguments{
\item{lx}{numeric. vector of the lifetable survivorship.}

\item{pert}{numeric scalar. The maximum size of the perturbation, default value of 0.0001}
}
\value{
a vector of lx, either perturbed or not.
}
\description{
Force \code{lx} to be monotonically decreasing. Ties in consecutive values of \code{lx} produce problems to estimate quantiles from \code{lx} in various of our functions. These either result in warnings or errors. This function checks for ties, and if found, it adds a small amount of monotonically decreasing noise to lx.
}
\details{
This function is here for robustness, as this situation often occurs in lifetables with a radix of 100000 and integer expressed output, especially around closeout ages. This perturbation is only applied if needed. The perturbation has a trivial effect of estimated quantile ages along most of the age range. The argument \code{pert} should be a small positive amount.
}
\examples{
lx <- 10:0 / 10
age <- 0:10
all(maybe_dither_lx(lx) == lx)
lx2 <- c(10,9,8,7,7,5,4,3,2,1,0) / 10
lx2 - maybe_dither_lx(lx2)
lx3 <- c(10:2,0,0)
lx3 - maybe_dither_lx(lx3)
}
