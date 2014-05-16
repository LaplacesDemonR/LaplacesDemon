\name{log-log}
\alias{cloglog}
\alias{invcloglog}
\alias{invloglog}
\alias{loglog}
\title{The log-log and complementary log-log functions}
\description{
  The log-log and complementary log-log functions, as well as the
  inverse functions, are provided.
}
\usage{
cloglog(p)
invcloglog(x)
invloglog(x)
loglog(p)
}
\arguments{
  \item{x}{This is a vector of real values that will be transformed to
    the interval [0,1].}
  \item{p}{This is a vector of probabilities p in the interval [0,1]
    that will be transformed to the real line.}
}
\details{
  The logit and probit links are symmetric, because the probabilities
  approach zero or one at the same rate. The log-log and complementary
  log-log links are asymmetric. Complementary log-log links approach
  zero slowly and one quickly. Log-log links approach zero quickly and
  one slowly. Either the log-log or complementary log-log link will tend
  to fit better than logistic and probit, and are frequently used when
  the probability of an event is small or large. A mixture of the two
  links, the log-log and complementary log-log is often used, where each
  link is weighted. The reason that logit is so prevalent is because
  logistic parameters can be interpreted as odds ratios.
}
\value{
  \code{cloglog} returns \code{x},
  \code{invcloglog} and \code{invloglog} return probability \code{p},
  and \code{loglog} returns \code{x}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{\code{\link{LaplacesDemon}}}
\examples{
library(LaplacesDemon)
x <- -5:5
p <- invloglog(x)
x <- loglog(p)
}
\keyword{Complementary log-log}
\keyword{Link Function}
\keyword{log-log}
\keyword{Transformation}
