\name{PMC}
\alias{PMC}
\title{Population Monte Carlo}
\description{
  The \code{PMC} function updates a model with Population Monte Carlo.
  Given a model specification, data, and initial values, \code{PMC}
  maximizes the logarithm of the unnormalized joint
  posterior density and provides samples of the marginal
  posterior distributions, deviance, and other monitored variables.
}
\usage{
PMC(Model, Data, Initial.Values, Covar=NULL, Iterations=10, Thinning=1,
alpha=NULL, M=1, N=1000, nu=9, CPUs=1, Type="PSOCK")
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{Initial.Values}{This is either a vector initial values, one for
    each of \eqn{K} parameters, or in the case of a mixture of \eqn{M}
    components, this is a \eqn{M \times K}{M x K} matrix of initial
    values. If all initial values are zero in this vector, or in the
    first row of a matrix, then \code{\link{LaplaceApproximation}} is
    used to optimize initial values, in which case all mixture
    components receive the same initial values and covariance matrix
    from the object of class \code{laplace}. Parameters must be
    continuous.}
  \item{Data}{This is a list of data. For more information, see
    \code{\link{LaplacesDemon}}.}
  \item{Covar}{This is a \eqn{K \times K}{K x K} covariance matrix for
    \eqn{K} parameters, or for multiple mixture components, this is a
    \eqn{K \times K \times M}{K x K x M} array of \eqn{M} covariance
    matrices, where \eqn{M} is the number of mixture components.
    \code{Covar} defaults to \code{NULL}, in which case a scaled
    identity matrix (with the same scale as in
    \code{\link{LaplacesDemon}}) is applied to all mixture components.}
  \item{Iterations}{This is the number of iterations during which PMC
    will update the model. Updating the model for only one iteration is
    the same as applying non-adaptive importance sampling.}
  \item{Thinning}{This is the number by which the posterior is
    thinned. To have 1,000 posterior samples with \code{M=3} mixture
    components and \code{N=10000} samples each, \code{Thinning=30}. For
    more information, see the \code{\link{Thin}} function.}
  \item{alpha}{This is a vector of length \eqn{M}, the number of mixture
    components. \eqn{\alpha}{alpha} is the probability of each mixture
    component. The default value is \code{NULL}, which assigns an equal
    probability to each component.}
  \item{M}{This is the number \eqn{M} of multivariate t distribution
    mixture components.}
  \item{N}{This is the number \eqn{N} of samples per mixture
    component. The required number of samples increases with the number
    \eqn{K} of parameters. These samples are also called walkers or
    particles.}
  \item{nu}{This is the degrees of freedom parameter \eqn{\nu}{nu} for
    the multivariate t distribution for each mixture component. If a
    multivariate normal distribution is preferred, then set
    \eqn{\nu > 1e4}{nu > 1e4}.}
  \item{CPUs}{This argument is required for parallel processing, and
    and indicates the number of central processing units (CPUs) of the
    computer or cluster. For example, when a user has a quad-core
    computer, \code{CPUs=4}.}
  \item{Type}{This argument defaults to \code{"PSOCK"} and uses the
    Simple Network of Workstations (SNOW) for parallelization.
    Alternatively, \code{Type="MPI"} may be specified to use Message
    Passing Interface (MPI) for parallelization.}
}
\details{
  The \code{PMC} function uses the adaptive importance sampling
  algorithm of Wraith et al. (2009), also called Mixture PMC or M-PMC
  (Cappe et al., 2008). Iterative adaptive importance sampling was
  introduced in the 1980s. Modern PMC was introduced (Cappe et al.,
  2004), and extended to multivariate Gaussian or t-distributed mixtures
  (Cappe et al., 2008). This version uses a multivariate t distribution
  for each mixture component, and also allows a multivariate normal
  distribution when the degrees of freedom, \eqn{\nu > 1e4}{nu >
    1e4}. At each iteration, a mixture distribution is sampled with
  importance sampling, and the samples (or populations) are adapted to
  improve the importance sampling. Adaptation is a variant of EM
  (Expectation-Maximization). The sample is self-normalized, and is an
  example of self-normalized importance sampling (SNIS), or
  self-importance sampling. The vector \eqn{\alpha}{alpha} contains the
  probability of each mixture component. These, as well as multivariate
  t distribution mixture parameters (except \eqn{\nu}{nu}), are adapted
  each iteration.

  Advantages of PMC over MCMC include:
  \itemize{
    \item It is difficult to assess convergence of MCMC chains, and this
  is not necessary in PMC (Wraith et al., 2009).  
    \item MCMC chains have autocorrelation that effectively reduces
  posterior samples. PMC produces independent samples that are not
  reduced with autocorrelation.
    \item PMC has been reported to produce samples with less variance
  than MCMC.
    \item It is difficult to parallelize MCMC. Posterior samples from
  parallel chains can be pooled when all chains have converged, but
  until this occurs, parallelization is unhelpful. PMC, on the other
  hand, can parallelize the independent, Monte Carlo samples during each
  iteration and reduce run-time as the number of processors
  increases. Currently, PMC is not parallelized here.
    \item The multivariate mixture in PMC can represent a multimodal
  posterior, where MCMC with parallel chains may be used to identify a
  multimodal posterior, but probably will not yield combined samples
  that proportionally represent it.
  }

  Disadvantages of PMC, compared to MCMC, include:
  \itemize{
    \item In PMC, the required number of samples at each iteration
  increases quickly with respect to an increase in parameters. MCMC is
  more suitable for models with large numbers of parameters, and
  therefore, MCMC is more generalizable.
    \item PMC is more sensitive to initial values than MCMC, especially
  as the number of parameters increases.
    \item PMC is more sensitive to the initial covariance matrix (or
  matrices for mixture components) than adaptive MCMC. PMC requires more
  information about the target distributions before updating. The
  covariance matrix from a converged iterative quadrature algorithm,
  Laplace Approximation, or Variational Bayes may be required (see
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  or \code{\link{VariationalBayes}} for more information).
  }

  Since PMC requires better initial information than iterative
  quadrature, Laplace Approximation, MCMC, and Variational Bayes, it is
  not recommended to begin updating a model that has little prior
  information with PMC, especially when the model has more than a few
  parameters. Instead, iterative quadrature, Laplace Approximation,
  MCMC, or Variational Bayes should be used. However, once convergence
  is found or assumed, it is recommended to attempt to update the model
  with PMC, given the latest parameters and convariance matrix from
  iterative quadrature, Laplace Approximation, MCMC, or Variational
  Bayes. Used in this way, PMC may improve the model fit obtained with
  MCMC and should reduce the variance of the marginal posterior
  distributions,  which is desirable for predictive modeling.
  
  Convergence is assessed by observing two outputs: normalized effective
  sample size (\code{ESSN}) and normalized perplexity
  (\code{Perplexity}). These are described below. PMC is considered to
  have converged when these diagnostics stabilize (Wraith et al., 2009),
  or when the normalized perplexity becomes sufficiently close to 1
  (Cappe et al., 2008). If they do not stabilize, then it is suggested
  to begin PMC again with a larger number \code{N} of samples, and
  possibly with different initial values and covariance matrix or
  matrices. \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, or \code{\link{VariationalBayes}}
  may be helpful to provide better starting values for \code{PMC}.

  If a message appears that warns about `bad weights', then \code{PMC}
  is attempting to work with an iteration in which importance weights
  are problematic. If this occurs in the first iteration, then all
  importance weights are set to \eqn{1/N}. If this occurs in other
  iterations, then the information from the previous iteration is used
  instead and different draws are made from that importance
  distribution. This may allow \code{PMC} to eventually adapt
  successfully to the target. If not, the user is advised to begin again
  with a larger number \eqn{N} of samples, and possibly different
  initial values and covariance matrix or matrices, as above. PMC can
  experience difficulty when it begins with poor initial conditions.

  The user may combine samples from previous iterations with samples
  from the latest iteration for inference, if the algorithm converged
  before the last iteration. Currently, a function is not provided for
  combining previous samples.
}
\value{
  The returned object is an object of class \code{pmc} with the
  following components:
  \item{alpha}{This is a \eqn{M \times T}{M x T} matrix of the
    probabilities of mixture components, where \eqn{M} is the number of
    mixture components and \eqn{T} is the number of iterations.}
  \item{Call}{This is the matched call of \code{PMC}.}
  \item{Covar}{This stores the \eqn{K \times K \times T \times M}{K x K
    x T x M} proposal covariance matrix in an array, where \eqn{K} is
    the dimension or number of parameters or initial values, \eqn{T} is
    the number of iterations, and \eqn{M} is the number of mixture
    components. If the model is updated in the future, then the latest
    covariance matrix for each mixture component can be extracted and
    used to start the next update where the last update left off.}
  \item{Deviance}{This is a vector of the deviance of the model, with a
    length equal to the number of thinned samples that were retained.
    Deviance is useful for considering model fit, and is equal to the
    sum of the log-likelihood for all rows in the data set, which is
    then multiplied by negative two.}
  \item{DIC}{This is a vector of three values: Dbar, pD, and DIC. Dbar
    is the mean deviance, pD is a measure of model complexity indicating
    the effective number of parameters, and DIC is the Deviance
    Information Criterion, which is a model fit statistic that is the
    sum of Dbar and pD. \code{DIC} is calculated over the thinned
    samples. Note that pD is calculated as \code{var(Deviance)/2} as in
    Gelman et al. (2004).}
  \item{ESSN}{This is a vector of length \eqn{T} that contains the
    normalized effective sample size (ESSN) per iteration across \eqn{T}
    iterations. ESSN is used as a convergence diagnostic. ESSN is
    normalized between zero and one, and can be interpreted as the
    proportion of samples with non-zero weight. Higher is better.}
  \item{Initial.Values}{This is the vector or matrix of
    \code{Initial.Values}.}
  \item{Iterations}{This reports the number of \code{Iterations} for
    updating.}
  \item{LML}{This is an approximation of the logarithm of the marginal
    likelihood of the data (see the \code{\link{LML}} function for more
    information). \code{LML} is estimated with nonparametric
    self-normalized importance sampling (NSIS), given LL and the
    marginal posterior samples of the parameters. \code{LML} is useful
    for comparing multiple models with the \code{\link{BayesFactor}}
    function.}
  \item{M}{This reports the number of mixture components.}
  \item{Minutes}{This indicates the number of minutes that \code{PMC}
    was running, and this includes the initial checks as well as time it
    took to perform final sampling and create summaries.}
  \item{Model}{This contains the model specification \code{Model}.}
  \item{N}{This is the number of un-thinned samples per mixture
    component.}
  item{nu}{This is the degrees of freedom parameter \eqn{\nu}{nu} for
  each multivariate t distribution in each mixture component.}
  \item{Mu}{This is a \eqn{T \times K \times M}{T x K x M} array of
    means for the importance sampling distribution across \eqn{T}
    iterations, \eqn{K} parameters, and \eqn{M} mixture components.}
  \item{Monitor}{This is a \eqn{S \times J}{S x J} matrix of
    thinned samples of monitored variables, where \eqn{S} is the number
    of thinned samples and \eqn{J} is the number of monitored
    variables.}
  \item{Parameters}{This reports the number \eqn{K} of parameters.}
  \item{Perplexity}{This is a vector of length \eqn{T} that contains the
    normalized perplexity per iteration across \eqn{T} iterations, and
    is used as a convergence diagnostic. Perplexity is an approximation
    of the negative of the Kullback-Leibler divergence (see
    \code{\link{KLD}}) between the target and the importance function.
    Perplexity is normalized between zero and one, and a higher
    normalized perplexity relates to less divergence, so higher is
    better. A normalized perplexity that is close to one indicates good
    agreement between the target density and the importance
    function. This is based on the Shannon entropy of the normalized
    importance weights, which is used frequently to measure the quality
    of importance samples.}
  \item{Posterior1}{This is an \eqn{N \times K \times T \times M}{N x K
    x T x M} array of un-thinned posterior samples across \eqn{N}
    samples, \eqn{K} parameters, \eqn{T} iterations, and \eqn{M} mixture
    components.}
  \item{Posterior2}{This is a \eqn{S \times K}{S x K} matrix of
    thinned posterior samples, where \eqn{S} is the number of thinned
    samples and \eqn{K} is the number of parameters.}
  \item{Summary}{This is a matrix that summarizes the marginal
    posterior distributions of the parameters, deviance, and monitored
    variables from thinned samples. The following summary statistics are
    included: mean, standard deviation, MCSE (Monte Carlo Standard
    Error), ESS is the effective sample size due to autocorrelation, and
    finally the 2.5\%, 50\%, and 97.5\% quantiles are reported. MCSE is
    essentially a standard deviation around the marginal posterior mean
    that is due to uncertainty associated with using Monte Carlo
    sampling. The acceptable size of the MCSE depends on the acceptable
    uncertainty associated around the marginal posterior mean. The
    default \code{IMPS} method is used. Next, the desired precision of
    ESS depends on the user's goal.}
  \item{Thinned.Samples}{This is the number of thinned samples in
    \code{Posterior2}.}
  \item{Thinning}{This is the amount of thinning requested by the user.}
  \item{W}{This is a \eqn{N \times T}{N x T} matrix of normalized
    importance weights, where \eqn{N} is the number of un-thinned
    samples per mixture component and \eqn{T} is the number of
    iterations. Computationally, the algorithm uses the logarithm of the
    weights.}
}
\references{
  Cappe, O., Douc, R., Guillin, A., Marin, J.M., and Robert, C. (2008).
  "Adaptive Importance Sampling in General Mixture Classes".
  \emph{Statistics and Computing}, 18, p. 587--600.

  Cappe, O., Guillin, A., Marin, J.M., and Robert, C. (2004). "Population 
  Monte Carlo". \emph{Journal of Computational and Graphical
  Statistics}, 13, p. 907--929.

  Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
  Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
  Hall, London.
  
  Wraith, D., Kilbinger, M., Benabed, K., Cappe, O., Cardoso, J.F.,
  Fort, G., Prunet, S., and Robert, C.P. (2009). "Estimation of
  Cosmological Parameters Using Adaptive Importance Sampling".
  \emph{Physical Review D}, 80(2), p. 023507.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BayesFactor}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LML}},
  \code{\link{PMC.RAM}},
  \code{\link{Thin}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
# The accompanying Examples vignette is a compendium of examples.
####################  Load the LaplacesDemon Library  #####################
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

#########################  Data List Preparation  #########################
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
     beta <- rnorm(Data$J)
     sigma <- runif(1)
     return(c(beta, sigma))
     }
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

##########################  Model Specification  ##########################
Model <- function(parm, Data)
     {
     ### Parameters
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
     parm[Data$pos.sigma] <- sigma
     ### Log-Prior
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     ### Log-Likelihood
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     ### Log-Posterior
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     return(Modelout)
     }

set.seed(666)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

########################  Population Monte Carlo  #########################
Fit <- PMC(Model, MyData, Initial.Values, Covar=NULL, Iterations=5,
     Thinning=1, alpha=NULL, M=1, N=100, CPUs=1)
Fit
print(Fit)
PosteriorChecks(Fit)
caterpillar.plot(Fit, Parms="beta")
plot(Fit, BurnIn=0, MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData, CPUs=1)
summary(Pred, Discrep="Chi-Square")
plot(Pred, Style="Covariates", Data=MyData)
plot(Pred, Style="Density", Rows=1:9)
plot(Pred, Style="ECDF")
plot(Pred, Style="Fitted")
plot(Pred, Style="Jarque-Bera")
plot(Pred, Style="Predictive Quantiles")
plot(Pred, Style="Residual Density")
plot(Pred, Style="Residuals")
Levene.Test(Pred)
Importance(Fit, Model, MyData, Discrep="Chi-Square")

#End
}
\keyword{High Performance Computing}
\keyword{Importance Sampling}
\keyword{Monte Carlo}