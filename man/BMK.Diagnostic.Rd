\name{BMK.Diagnostic}
\alias{BMK.Diagnostic}
\title{BMK Convergence Diagnostic}
\description{
  Given a matrix of posterior samples from MCMC, the
  \code{BMK.Diagnostic} function calculates Hellinger distances between
  consecutive batches for each chain. This is useful for monitoring
  convergence of MCMC chains.
}
\usage{
BMK.Diagnostic(X, batches=10)
}
\arguments{
  \item{X}{This required argument accepts a matrix of posterior
    samples or an object of class \code{demonoid}, in which case it uses
    the posterior samples in \code{X$Posterior1}.}
  \item{batches}{This is the number of batches on which the convergence
    diagnostic will be calculated. The \code{batches} argument defaults
    to 10.}
}
\details{
  Hellinger distance is used to quantify dissimilarity between two
  probability distributions. It is based on the Hellinger integral,
  introduced by Hellinger (1909). Traditionally, Hellinger distance is
  bound to the interval [0,1], though another popular form occurs in
  the interval [0,\eqn{\sqrt{2}}{sqrt(2)}]. A higher value of Hellinger
  distance is associated with more dissimilarity between the
  distributions.

  Convergence is assumed when Hellinger distances are below a threshold,
  indicating that posterior samples are similar between consecutive
  batches. If all Hellinger distances beyond a given batch of samples is
  below the threshold, then \code{burnin} is suggested to occur
  immediately before the first batch of satisfactory Hellinger
  distances.

  As an aid to interpretation, consider a matrix of 1,000 posterior
  samples from three chains: \code{beta[1]}, \code{beta[2]}, and
  \code{beta[3]}. With 10 batches, the column names are: 100, 200,
  \dots, 900. A Hellinger distance for the chain \code{beta[1]} at 100
  is the Hellinger distance between two batches: samples 1-100, and
  samples 101:200.

  A benefit to using \code{BMK.Diagnostic} is that the resulting
  Hellinger distances may easily be plotted with the \code{plotMatrix}
  function, allowing the user to see quickly which consecutive batches
  of which chains were dissimilar. This makes it easier to find
  problematic chains.

  The \code{BMK.Diagnostic} is calculated automatically in the
  \code{\link{LaplacesDemon}} function, and is one of the criteria in
  the \code{\link{Consort}} function regarding the recommendation of
  when to stop updating the Markov chain Monte Carlo (MCMC) sampler in
  \code{\link{LaplacesDemon}}.

  For more information on the related topics of burn-in and
  stationarity, see the \code{burnin} and \code{\link{is.stationary}}
  functions, and the accompanying vignettes.
}
\value{
  The \code{BMK.Diagnostic} function returns an object of class
  \code{bmk} that is a \eqn{J \times B}{J x B} matrix of Hellinger
  distances between consecutive batches for \eqn{J} parameters of
  posterior samples. The number of columns, \eqn{B} is equal to the
  number of batches minus one.

  The \code{BMK.Diagnostic} function is similar to the
  \code{bmkconverge} function in package BMK.
}
\references{
  Boone, E.L., Merrick, J.R. and Krachey, M.J. (2013). "A Hellinger
  Distance Approach to MCMC Diagnostics". \emph{Journal of Statistical
  Computation and Simulation}, in press.

  Hellinger, E. (1909). "Neue Begrundung der Theorie quadratischer
  Formen von unendlichvielen Veranderlichen" (in German). \emph{Journal
  fur die reine und angewandte Mathematik}, 136, p. 210--271.
}
\seealso{
  \code{\link{burnin}},
  \code{\link{Consort}},
  \code{\link{is.stationary}}, and 
  \code{\link{LaplacesDemon}}.
}
\examples{
library(LaplacesDemon)
N <- 1000 #Number of posterior samples
J <- 10 #Number of parameters
Theta <- matrix(runif(N*J),N,J)
colnames(Theta) <- paste("beta[", 1:J, "]", sep="")
for (i in 2:N) {Theta[i,1] <- Theta[i-1,1] + rnorm(1)}
HD <- BMK.Diagnostic(Theta, batches=10)
plot(HD, title="Hellinger distance between batches")
}
\keyword{Diagnostic}
\keyword{MCMC}
\keyword{Stationarity}
\keyword{Utility}
