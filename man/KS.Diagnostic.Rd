\name{KS.Diagnostic}
\alias{KS.Diagnostic}
\title{Kolmogorov-Smirnov Convergence Diagnostic}
\description{
  The Kolmogorov-Smirnov test is a nonparametric test of stationarity
  that has been applied as an MCMC diagnostic (Brooks et al, 2003), such
  as to the posterior samples from the \code{\link{LaplacesDemon}}
  function. The first and last halves of the chain are compared. This
  test assumes IID, which is violated in the presence of
  autocorrelation.

  The \code{KS.Diagnostic} is a univariate diagnostic that is usually
  applied to each marginal posterior distribution. A multivariate form
  is not included. By chance alone due to multiple independent tests,
  5\% of the marginal posterior distributions should appear
  non-stationary when stationarity exists. Assessing multivariate 
  convergence is difficult.
}
\usage{
KS.Diagnostic(x)
}
\arguments{
  \item{x}{This is a vector of posterior samples for which a
    Kolmogorov-Smirnov test will be applied that compares the first and
    last halves for stationarity.}
}
\details{
  There are two main approaches to using the Kolmogorov-Smirnov test as
  an MCMC diagnostic. There is a version of the test that has
  been adapted to account for autocorrelation (and is not included
  here). Otherwise, the chain is thinned enough that autocorrelation is
  not present or is minimized, in which case the two-sample
  Kolmogorov-Smirnov test is applied. The CDFs of both samples are
  compared. The \code{ks.test} function in base R is used.

  The advantage of the Kolmogorov-Smirnov test is that it is easier and
  faster to calculate. The disadvantages are that autocorrelation biases
  results, and the test is generally biased on the conservative side
  (indicating stationarity when it should not).
}
\value{
  The \code{KS.Diagnostic} function returns a frequentist p-value, and
  stationarity is indicated when p > 0.05.
}
\references{
  Brooks, S.P., Giudici, P., and Philippe, A. (2003). "Nonparametric
  Convergence Assessment for MCMC Model Selection". \emph{Journal of
  Computational and Graphical Statistics}. 12(1), p. 1--22.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{is.stationary}},
  \code{\link{ks.test}}, and
  \code{\link{LaplacesDemon}}.
}
\examples{
library(LaplacesDemon)
x <- rnorm(1000)
KS.Diagnostic(x)
}
\keyword{Diagnostic}
\keyword{MCMC}
