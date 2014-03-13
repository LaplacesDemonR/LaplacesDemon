\name{WAIC}
\alias{WAIC}
\title{Widely Applicable Information Criterion}
\description{
  This function calculates the Widely Applicable Information Criterion
  (WAIC), also known as the Widely Available Information Criterion or
  the Watanable-Akaike, of Watanabe (2010).
}
\usage{
WAIC(x)
}
\arguments{
  \item{x}{This required argument accepts a \eqn{N \times S}{N x S}
    matrix of log-likelihood (LL) for \eqn{N} records and \eqn{S}
    samples.}
}
\details{
  WAIC is an extension of the Akaike Information Criterion (AIC) that is
  more fully Bayesian than the Deviance Information Criterion (DIC).

  Like DIC, WAIC estimates the effective number of parameters to adjust
  for overfitting. Two adjustments have been proposed. pWAIC1 is similar
  to pD in the original DIC. In contrast, pWAIC2 is calculated with
  variance more similarly to pV, which Gelman proposed for DIC. Gelman
  et al. (2014, p.174) recommends pWAIC2 because its results are closer
  in practice to the results of leave-one-out cross-validation
  (LOO-CV). pWAIC is considered an approximation to the number of
  unconstrained and uninformed parameters, where a parameter counts as 1
  when estimated without contraint or any prior information, 0 if fully
  constrained or all information comes from the prior distribution, or
  an intermediate number if both the data and prior are informative,
  which is usually the case.

  Gelman et al. (2014, p. 174) scale the WAIC of Watanabe (2010) by a
  factor of 2 so that it is comparable to AIC and DIC. WAIC is then
  reported as \eqn{-2(lppd - pWAIC)}. Gelman et al. (2014) prefer WAIC
  to AIC or DIC when feasible, which is less often than AIC or DIC.
  The \code{\link{LaplacesDemon}} function requires the model
  specification function to return the model-level deviance, which is
  \eqn{-2(LL)}, where \eqn{LL} is the sum of the record-level
  log-likelihood. Therefore, if the user desires to calculate WAIC, then
  the record-level log-likelihood must be monitored.
}
\value{
  The \code{WAIC} argument returns a list with four components:
  \item{WAIC}{This is the Widely Applicable Information Criterion
  (WAIC), which is \eqn{-2(lppd - pWAIC)}.}
  \item{lppd}{This is the log pointwise predictive density. For more
  information, see Gelman et al. (2014, p. 168).}
  \item{pWAIC}{This is the effective number of parameters preferred by
    Gelman et al. (2014).}
  \item{pWAIC1}{This is the effective number of parameters, is
    calculated with an alternate method, and is included here for
    completeness. It is not used to calculate WAIC in the \code{WAIC}
    function.}
}
\references{
  Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., and
  Rubin, D.B. (2014). "Bayesian Data Analysis, 3rd ed.". CRC Press:
  Boca Raton, FL.
  
  Watanabe, S. (2010). "Asymptotic Equivalence of Bayes Cross Validation
  and Widely Applicable Information Criterion in Singular Learning
  Theory". \emph{Journal of Machine Learning Research}, 11,
  p. 3571--3594.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}}
}
\examples{
#library(LaplacesDemon)
#N <- 10 #Number of records
#S <- 1000 #Number of samples
#LL <- t(rmvn(S, -70+rnorm(N),
#     as.positive.definite(matrix(rnorm(N*N),N,N))))
#WAIC(LL)
### Compare with DIC:
#Dev <- -2*colSums(LL)
#DIC <- list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)
#DIC
}
\keyword{Diagnostic}
\keyword{Information Criterion}
\keyword{MCMC}
\keyword{Monte Carlo}
\keyword{Utility}
