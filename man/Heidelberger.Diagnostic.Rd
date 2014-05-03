\name{Heidelberger.Diagnostic}
\alias{Heidelberger.Diagnostic}
\title{Heidelberger and Welch's MCMC Convergence Diagnostic}
\description{
  Heidelberger and Welch (1981; 1983) proposed a two-part MCMC
  convergence diagnostic that calculates a test statistic (based on the
  Cramer-von Mises test statistic) to accept or reject the null
  hypothesis that the Markov chain is from a stationary distribution.
}
\usage{Heidelberger.Diagnostic(x, eps=0.1, pvalue=0.05)}
\arguments{
  \item{x}{This required argument accepts an object of class
    \code{demonoid}. It attempts to use \code{Posterior2}, but when this
    is missing it uses \code{Posterior1}.}
  \item{eps}{This argument specifies the target value for the ratio of
    halfwidth to sample mean.}
  \item{pvalue}{This argument specifies the level of statistical
    significance.}
}
\details{
  The Heidelberg and Welch MCMC convergence diagnostic consists of two
  parts:

  First Part
  1. Generate a chain of \eqn{N} iterations and define an alpha level.
  2. Calculate the test statistic on the whole chain. Accept or reject
  the null hypothesis that the chain is from a stationary distribution.
  3. If the null hypothesis is rejected, then discard the first 10\% of
  the chain. Calculate the test statistic and accept or reject the null
  hypothesis.
  4. If the null hypothesis is rejected, then discard the next 10\% and
  calculate the test statistic.
  5. Repeat until the null hypothesis is accepted or 50\% of the chain
  is discarded. If the test still rejects the null hypothesis, then
  the chain fails the test and needs to be run longer.

  Second Part
  If the chain passes the first part of the diagnostic, then the part of
  the chain that was not discarded from the first part is used to test
  the second part.

  The halfwidth test calculates half the width of the (1 - alpha)\%
  probability interval (credible interval) around the mean.

  If the ratio of the halfwidth and the mean is lower than \code{eps},
  then the chain passes the halfwidth test. Otherwise, the chain fails
  the halfwidth test and must be updated for more iterations until
  sufficient accuracy is obtained. In order to avoid problems caused by
  sequential testing, the test should not be repeated too frequently.
  Heidelberger and Welch (1981) suggest increasing the run length by a
  factor I > 1.5, each time, so that estimate has the same, reasonably
  large, proportion of new data.

  The Heidelberger and Welch MCMC convergence diagnostic conducts
  multiple hypothesis tests. The number of potentially wrong results
  increases with the number of non-independent hypothesis tests
  conducted.

  The \code{Heidelberger.Diagnostic} is a univariate diagnostic that is
  usually applied to each marginal posterior distribution. A
  multivariate form is not included. By chance alone due to multiple
  independent tests, 5\% of the marginal posterior distributions should
  appear non-stationary when stationarity exists. Assessing multivariate
  convergence is difficult.
}
\note{
  The \code{Heidelberger.Diagnostic} function was adapted from the
  \code{heidel.diag} function in the coda package.
}
\value{
  The \code{Heidelberger.Diagnostic} function returns an object of class
  \code{heidelberger}. This object is a \eqn{J \times 6}{J x 6}
  matrix, and it is intended to be summarized with the
  \code{\link{print.heidelberger}} function. Nonetheless, this object of
  class \code{heidelberger} has \eqn{J} rows, each of which corresponds
  to a Markov chain. The column names are \code{stest}, \code{start},
  \code{pvalue}, \code{htest}, \code{mean}, and \code{halfwidth}. The
  \code{stest} column indicates convergence with a one, and
  non-convergence with a zero, regarding the stationarity test. When
  non-convergence is indicated, the remaining columns have missing
  values. The \code{start} column indicates the starting iteration, and
  the \code{pvalue} column shows the p-value associated with the first
  test. The \code{htest} column indicates convergence for the halfwidth
  test. The \code{mean} and \code{halfwidth} columns report the mean and
  halfwidth.
}
\references{
Heidelberger, P. and Welch, P.D. (1981). "A Spectral Method for
Confidence Interval Generation and Run Length Control in Simulations".
\emph{Comm. ACM.}, 24, p. 233--245.

Heidelberger, P. and Welch, P.D. (1983). "Simulation Run Length Control
in the Presence of an Initial Transient". \emph{Opns Res.}, 31,
p. 1109--1144.

Schruben, L.W. (1982). "Detecting Initialization Bias in Simulation
Experiments". \emph{Opns. Res.}, 30, p. 569--590.
}
\seealso{
  \code{\link{burnin}},
  \code{\link{is.stationary}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{print.heidelberger}}.
}
\examples{
#library(LaplacesDemon)
###After updating with LaplacesDemon, do:
#hd <- Heidelberger.Diagnostic(Fit)
#print(hd)
}
\keyword{Diagnostic}
\keyword{MCMC}
