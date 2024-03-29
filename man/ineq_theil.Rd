% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indices.R
\name{ineq_theil}
\alias{ineq_theil}
\title{ineq_theil}
\usage{
ineq_theil(age, dx, ex, ax, distribution_type = c("aad", "rl"), check = TRUE)
}
\arguments{
\item{age}{numeric. vector of lower age bounds.}

\item{dx}{numeric. vector of the lifetable death distribution.}

\item{ex}{numeric. vector of remaining life expectancy.}

\item{ax}{numeric. vector of the average time spent in the age}

\item{distribution_type}{character. Either \code{"aad"} (age at death) or \code{"rl"} (remaining life)}

\item{check}{logical. Shall we perform basic checks on input vectors? Default TRUE}
}
\description{
Calculate a lifetable column for the conditional Theil index of inequality in survivorship
}
\details{
All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in \code{MortalityLaws::MortalityLaw}). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package \code{ungroup} or with a penalized B-spline approach in package \code{MortalitySmooth}).
}
\examples{
data(LT)
# A vector containing the conditional Theil indices
Ta = ineq_theil(age=LT$Age,dx=LT$dx,ex=LT$ex,ax=LT$ax, distribution_type = "aad")
# The Theil index from birth
Ta[1]
# The Theil index conditional upon survival to age 10
Ta[11]

# A shortfall (remaining years) version of the same:
Tr = ineq_theil(age=LT$Age,dx=LT$dx,ex=LT$ex,ax=LT$ax, distribution_type = "rl")
Tr[1]
Tr[11]
\dontrun{
plot(0:110, Tr, type='l',col="red",ylab="conditional Theil",xlab="Age")
lines(0:110, Ta, col = "blue")
legend("topleft",col = c("red","blue"), lty=1,legend = c("remaining life","age at death"))
}
}
\references{
\insertRef{theil1967economics}{LifeIneq}
\insertRef{vanraalte2012}{LifeIneq}
\insertRef{hakkert1987}{LifeIneq}
\insertRef{cowell1980}{LifeIneq}
}
\seealso{
\code{MortalityLaws::\link[MortalityLaws]{MortalityLaw}}

\code{ungroup::\link[ungroup]{pclm}}

\code{MortalitySmooth::\link[MortalitySmooth]{Mort1Dsmooth}}
}
