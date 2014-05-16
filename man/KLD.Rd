\name{KLD}
\alias{KLD}
\title{Kullback-Leibler Divergence (KLD)}
\description{
  This function calculates the Kullback-Leibler divergence (KLD) between
  two probability distributions, and has many uses, such as in lowest
  posterior loss probability intervals, posterior predictive checks,
  prior elicitation, reference priors, and Variational Bayes.
}
\usage{
KLD(px, py, base)
}
\arguments{
   \item{px}{This is a required vector of probability densities,
     considered as \eqn{p(\textbf{x})}{p(x)}. Log-densities are also
     accepted, in which case both \code{px} and \code{py} must be
     log-densities.}
   \item{py}{This is a required vector of probability densities,
     considered as \eqn{p(\textbf{y})}{p(y)}. Log-densities are also
     accepted, in which case both \code{px} and \code{py} must be
     log-densities.}
   \item{base}{This optional argument specifies the logarithmic base,
     which defaults to \code{base=exp(1)} (or \eqn{e}) and represents
     information in natural units (nats), where \code{base=2} represents
     information in binary units (bits).}
}
\details{
  The Kullback-Leibler divergence (KLD) is known by many names, some of
  which are Kullback-Leibler distance, K-L, and logarithmic divergence.
  KLD is an asymmetric measure of the difference, distance, or direct
  divergence between two probability distributions
  \eqn{p(\textbf{y})}{p(y)} and \eqn{p(\textbf{x})}{p(x)} (Kullback and
  Leibler, 1951). Mathematically, however, KLD is not a distance,
  because of its asymmetry.
  
  Here, \eqn{p(\textbf{y})}{p(y)} represents the
  ``true'' distribution of data, observations, or theoretical
  distribution, and \eqn{p(\textbf{x})}{p(x)} represents a theory,
  model, or approximation of \eqn{p(\textbf{y})}{p(y)}.

  For probability distributions \eqn{p(\textbf{y})}{p(y)} and
  \eqn{p(\textbf{x})}{p(x)} that are discrete (whether the underlying
  distribution is continuous or discrete, the observations themselves
  are always discrete, such as from \eqn{i=1,\dots,N}{i=1,...,N}),

  \deqn{\mathrm{KLD}[p(\textbf{y}) || p(\textbf{x})] = \sum^N_i
  p(\textbf{y}_i)
  \log\frac{p(\textbf{y}_i)}{p(\textbf{x}_i)}}{KLD[p(y)||p(x)] = sum of
  p(y[i]) log(p(y[i]) / p(x[i]))}

  In Bayesian inference, KLD can be used as a measure of the information
  gain in moving from a prior distribution, \eqn{p(\theta)}{p(theta)},
  to a posterior distribution, \eqn{p(\theta | \textbf{y})}{p(theta |
  y)}. As such, KLD is the basis of reference priors and lowest
  posterior loss intervals (\code{\link{LPL.interval}}), such as in
  Berger, Bernardo, and Sun (2009) and Bernardo (2005). The intrinsic
  discrepancy was introduced by Bernardo and Rueda (2002). For more
  information on the intrinsic discrepancy, see
  \code{\link{LPL.interval}}.
}
\value{
  \code{KLD} returns a list with the following components:
  \item{KLD.px.py}{This is \eqn{\mathrm{KLD}_i[p(\textbf{x}_i) ||
      p(\textbf{y}_i)]}{KLD[i](p(x[i]) || p(y[i]))}.}
  \item{KLD.py.px}{This is \eqn{\mathrm{KLD}_i[p(\textbf{y}_i) ||
      p(\textbf{x}_i)]}{KLD[i](p(y[i]) || p(x[i]))}.}
  \item{mean.KLD}{This is the mean of the two components above. This is
      the expected posterior loss in \code{\link{LPL.interval}}.}
  \item{sum.KLD.px.py}{This is \eqn{\mathrm{KLD}[p(\textbf{x}) ||
      p(\textbf{y})]}{KLD(p(x) || p(y))}. This is a directed
      divergence.}
  \item{sum.KLD.py.px}{This is \eqn{\mathrm{KLD}[p(\textbf{y}) ||
      p(\textbf{x})]}{KLD(p(y) || p(x))}. This is a directed divergence.}
  \item{mean.sum.KLD}{This is the mean of the two components above.}
  \item{intrinsic.discrepancy}{This is minimum of the two directed
      divergences.}
}
\references{
  Berger, J.O., Bernardo, J.M., and Sun, D. (2009). "The Formal
  Definition of Reference Priors". \emph{The Annals of Statistics},
  37(2), p. 905--938.

  Bernardo, J.M. and Rueda, R. (2002). "Bayesian Hypothesis Testing: A
  Reference Approach". \emph{International Statistical Review}, 70,
  p. 351--372.

  Bernardo, J.M. (2005). "Intrinsic Credible Regions: An Objective
  Bayesian Approach to Interval Estimation". \emph{Sociedad de
  Estadistica e Investigacion Operativa}, 14(2), p. 317--384.

  Kullback, S. and Leibler, R.A. (1951). "On Information and
  Sufficiency". \emph{The Annals of Mathematical Statistics}, 22(1),
  p. 79--86.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LPL.interval}} and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
px <- dnorm(runif(100),0,1)
py <- dnorm(runif(100),0.1,0.9)
KLD(px,py)
}
\keyword{Distribution}
\keyword{Elicitation}
\keyword{Reference Priors}