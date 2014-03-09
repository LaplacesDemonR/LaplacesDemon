\name{BayesianBootstrap}
\alias{BayesianBootstrap}
\title{The Bayesian Bootstrap}
\description{
  This function performs the Bayesian Bootstrap of Rubin (1981),
  returning either bootstrapped weights or statistics.
}
\usage{
BayesianBootstrap(X, n=1000, Method="weights", Status=NULL)
}
\arguments{
  \item{X}{This is a vector or matrix of data. When a matrix is
    supplied, sampling is based on the first column.}
  \item{n}{This is the number of bootstrapped replications.}
  \item{Method}{When \code{Method="weights"} (which is the default),
    a matrix of row weights is returned. Otherwise, a function is
    accepted. The function specifies the statistic to be
    bootstrapped. The first argument of the function should be a matrix
    of data, and the second argument should be a vector of weights.}
  \item{Status}{This determines the periodicity of status messages. When
    \code{Status=100}, for example, a status message is displayed every
    100 replications. Otherwise, \code{Status} defaults to \code{NULL},
    and status messages are not displayed.}
}
\details{
  The term, `bootstrap', comes from the German novel \emph{Adventures of
  Baron Munchausen} by Rudolph Raspe, in which the hero saves himself
  from drowning by pulling on his own bootstraps. The idea of the
  statistical bootstrap is to evaluate properties of an estimator
  through the empirical, rather than theoretical, CDF.
  
  Rubin (1981) introduced the Bayesian bootstrap. In contrast to the
  frequentist bootstrap which simulates the sampling distribution of a
  statistic estimating a parameter, the Bayesian bootstrap simulates the
  posterior distribution.

  The data, \eqn{\textbf{X}}{X}, are assumed to be independent and
  identically distributed (IID), and to be a representative sample of
  the larger (bootstrapped) population. Given that the data has \eqn{N}
  rows in one bootstrap replication, the row weights are sampled from a
  Dirichlet distribution with all \eqn{N} concentration parameters equal
  to \eqn{1} (a uniform distribution over an open standard \eqn{N-1}
  simplex). The distributions of a parameter inferred from considering
  many samples of weights are interpretable as posterior distributions
  on that parameter.

  Just as with the frequentist bootstrap, inappropriate use of the
  Bayesian bootstrap can lead to inappropriate inferences. The Bayesian
  bootstrap violates the likelihood principle, because the evaluation of
  a statistic of interest depends on data sets other than the observed
  data set. For more information on the likelihood principle, see
  \url{http://www.bayesian-inference.com/likelihood#likelihoodprinciple}.

  The \code{BayesianBootstrap} function has many uses, including
  creating test statistics on the population data given the observed
  data (supported here), imputation (with this variation:
  \code{\link{ABB}}), validation, and more.
}
\value{
  When \code{Method="weights"}, this function returns a
  \eqn{N \times n}{S x J} matrix of weights, where the number of rows
  \eqn{N} is equal to the number of rows in \code{X}.

  For statistics, a matrix or array is returned, depending on the number
  of dimensions. The replicates are indexed by row in a matrix or in
  the first dimension of the array.
}
\references{
  Rubin, D.B. (1981). "The Bayesian Bootstrap". \emph{The Annals of
  Statistics}, 9(1), p. 130--134.
}
\author{Bogumil Kaminski, \email{bkamins@sgh.waw.pl} and
  Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{ABB}},
  \code{\link{LaplaceApproximation}}, and
  \code{\link{LaplacesDemon}}.
}
\examples{
#Example 1: Samples
x <- 1:2
BB <- BayesianBootstrap(X=x, n=100, Method="weights"); BB

#Example 2: Mean, Univariate
x <- 1:2
BB <- BayesianBootstrap(X=x, n=100, Method=weighted.mean); BB

#Example 3: Mean, Multivariate
data(demonsnacks)
BB <- BayesianBootstrap(X=demonsnacks, n=100,
     Method=function(x,w) apply(x, 2, weighted.mean, w=w)); BB

#Example 4: Correlation
dye <- c(1.15, 1.70, 1.42, 1.38, 2.80, 4.70, 4.80, 1.41, 3.90)
efp <- c(1.38, 1.72, 1.59, 1.47, 1.66, 3.45, 3.87, 1.31, 3.75)
X <- matrix(c(dye,efp), length(dye), 2)
colnames(X) <- c("dye","efp")
BB <- BayesianBootstrap(X=X, n=100,
     Method=function(x,w) cov.wt(x, w, cor=TRUE)$cor); BB
}
\keyword{Utility}
