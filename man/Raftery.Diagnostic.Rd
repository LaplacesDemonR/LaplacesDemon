\name{Raftery.Diagnostic}
\alias{Raftery.Diagnostic}
\title{Raftery and Lewis's diagnostic}
\description{
  Raftery and Lewis (1992) introduced an MCMC diagnostic that estimates
  the number of iterations needed for a given level of precision in
  posterior samples, as well as estimating burn-in, when quantiles are
  the posterior summaries of interest.
}
\usage{
Raftery.Diagnostic(x, q=0.025, r=0.005, s=0.95, eps=0.001)
}
\arguments{
  \item{x}{This required argument accepts an object of class
    \code{demonoid}. It attempts to use \code{Posterior2}, but when this
    is missing it uses \code{Posterior1}.}
  \item{q}{This is the quantile to be estimated.}
  \item{r}{This is the desired margin of error of the estimate, also
    called the accuracy.}
  \item{s}{This is the probability of obtaining an estimate in the
    interval (q-r, q+r).}
  \item{eps}{This is the precision required for the estimate of time to
    convergence.}
}
\details{
  In this MCMC diagnostic, a posterior quantile \eqn{q} of interest is
  specified. Next, an acceptable tolerance \eqn{r} is specified for
  \eqn{q}, which means that it is desired to measure \eqn{q} with an
  accuracy of +/- \eqn{r}. Finally, the user selects a probability
  \eqn{s}, which is the probability of being within the interval
  \eqn{(q-r, q+r)}. The \code{Raftery.Diagnostic} then estimates the
  number \eqn{N} of iterations and the number \eqn{M} of burn-in
  iterations that are necessary to satisfy the specified conditions
  regarding quantile \eqn{q}.

  The diagnostic was designed to test a short, initial update, in which
  the chains were called pilot chains, and the application was later
  suggested for iterative use after any update as a general method for
  pursuing convergence (Raftery and Lewis, 1996).

  Results of the \code{Raftery.Diagnostic} differ depending on the
  chosen quantile \eqn{q}. Estimates are conservative, so more
  iterations are suggested than necessary.
}
\note{
  The \code{Raftery.Diagnostic} function was adapted from the
  \code{raftery.diag} function in the coda package, which was adapted
  from the FORTRAN program `gibbsit', written by Steven Lewis.
}
\value{
  The \code{Raftery.Diagnostic} function returns an object of class
  \code{raftery} that is list. A print method is available for objects
  of this class. The list has the following components:
  \item{tspar}{These are the time-series parameters of the posterior
    samples in \code{x}.}
  \item{params}{This is a vector containing the parameters \code{q},
    \code{r}, and \code{s}.}
  \item{Niters}{This is the number of iterations in the posterior
    samples in \code{x}.}
  \item{resmatrix}{This is a 3-dimensional array containing the
    results: \eqn{M} is the suggested burn-in, \eqn{N} is the suggested
    number of iterations, \eqn{Nmin} is the suggested number of
    iterations based on zero autocorrelation, and
    \eqn{I = (M+N)/Nmin} is the "dependence factor". The dependence
    factor is interpreted as the proportional increase in the number of
    iterations attributable to autocorrelation. Highly autocorrelated
    chains (> 5) are worrisome, and may be due to influential initial
    values, parameter correlations, or poor mixing.}
}
\references{
  Raftery, A.E. and Lewis, S.M. (1992). "How Many Iterations in the
  Gibbs Sampler?" In \emph{Bayesian Statistics}, 4 (J.M. Bernardo, J.O.
  Berger, A.P. Dawid and A.F.M. Smith, eds.). Oxford, U.K.: Oxford
  University Press, p. 763--773.

  Raftery, A.E. and Lewis, S.M. (1992). "One Long Run with
  Diagnostics: Implementation Strategies for Markov chain Monte Carlo".
  \emph{Statistical Science}, 7, p. 493--497.

  Raftery, A.E. and Lewis, S.M. (1996). "Implementing MCMC".
  \emph{In} Practical Markov Chain Monte Carlo (W.R. Gilks,
  D.J. Spiegelhalter and S. Richardson, eds.). Chapman and Hall:
  Baton Rouge, FL.
}
 \seealso{
  \code{\link{burnin}},
  \code{\link{LaplacesDemon}},
  \code{\link{print.raftery}}, and
  \code{\link{Thin}}.
}
\examples{
#library(LaplacesDemon)
###After updating with LaplacesDemon, do:
#rd <- Raftery.Diagnostic(Fit)
#print(rd)
}
\keyword{Diagnostic}
\keyword{Gibbsit}
\keyword{MCMC}
