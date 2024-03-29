% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indices.R
\name{ineq_iqr}
\alias{ineq_iqr}
\title{ineq_iqr}
\usage{
ineq_iqr(age, lx, upper = 0.75, lower = 0.25)
}
\arguments{
\item{age}{numeric. vector of lower age bounds.}

\item{lx}{numeric. vector of the lifetable survivorship.}

\item{upper}{numeric. upper survival quantile, default .75}

\item{lower}{numeric. lower survival quantile, defauly .25}
}
\description{
Calculate the interquartile range survivorship age from a lifetable. Other quantile ranges can also be calculated.
}
\details{
All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in \code{MortalityLaws::MortalityLaw}). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package \code{ungroup} or with a penalized B-spline approach in package \code{MortalitySmooth}).
}
\examples{

data(LT)
# The iqr range of survival age
iqr <- ineq_iqr(age=LT$Age,lx=LT$lx)
iqr
# age distance between 10th and 90th percentiles
idr <- ineq_iqr(age=LT$Age,lx=LT$lx,upper=.9,lower=.1)
idr
}
\seealso{
\code{MortalityLaws::\link[MortalityLaws]{MortalityLaw}}

\code{ungroup::\link[ungroup]{pclm}}

\code{MortalitySmooth::\link[MortalitySmooth]{Mort1Dsmooth}}
}
