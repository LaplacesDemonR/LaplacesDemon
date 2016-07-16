\name{BayesFactor}
\alias{BayesFactor}
\title{Bayes Factor}
\description{
  This function calculates Bayes factors for two or more fitted objects
  of class \code{demonoid}, \code{iterquad}, \code{laplace}, \code{pmc},
  or \code{vb} that were estimated respectively with the
  \code{\link{LaplacesDemon}}, \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, \code{\link{PMC}}, or
  \code{\link{VariationalBayes}} functions, and indicates the strength
  of evidence in favor of the hypothesis (that each model,
  \eqn{\mathcal{M}_i}{M[i]}, is better than another model,
  \eqn{\mathcal{M}_j}{M[j]}).
}
\usage{
BayesFactor(x)
}
\arguments{
  \item{x}{This is a list of two or more fitted objects of class
    \code{demonoid}, \code{iterquad}, \code{laplace}, \code{pmc}, or
    \code{vb}. The components are named in order beginning with model 1,
    \code{M1}, and \eqn{k} models are usually represented as
    \eqn{\mathcal{M}_1,\dots,\mathcal{M}_k}{M[1],...,M[k]}.}
}
\details{
  Introduced by Harold Jeffreys, a 'Bayes factor' is a Bayesian
  alternative to frequentist hypothesis testing that is most often used
  for the comparison of multiple models by hypothesis testing, usually
  to determine which model better fits the data (Jeffreys, 1961). Bayes
  factors are notoriously difficult to compute, and the Bayes factor is
  only defined when the marginal density of \eqn{\textbf{y}}{y} under
  each model is proper (see \code{\link{is.proper}}). However, the Bayes
  factor is easy to approximate with the Laplace-Metropolis estimator
  (Lewis and Raftery, 1997) and other methods of approximating the
  logarithm of the marginal likelihood (for more information, see
  \code{\link{LML}}).

  Hypothesis testing with Bayes factors is more robust than frequentist
  hypothesis testing, since the Bayesian form avoids model selection
  bias, evaluates evidence in favor of the null hypothesis, includes
  model uncertainty, and allows non-nested models to be compared (though
  of course the model must have the same dependent variable). Also,
  frequentist significance tests become biased in favor of rejecting the
  null hypothesis with sufficiently large sample size.

  The Bayes factor for comparing two models may be approximated as the
  ratio of the marginal likelihood of the data in model 1 and model
  2. Formally, the Bayes factor in this case is

\deqn{B =
  \frac{p(\textbf{y}|\mathcal{M}_1)}{p(\textbf{y}|\mathcal{M}_2)} =
  \frac{\int
    p(\textbf{y}|\Theta_1,\mathcal{M}_1)p(\Theta_1|\mathcal{M}_1)d\Theta_1}{\int
    p(\textbf{y}|\Theta_2,\mathcal{M}_2)p(\Theta_2|\mathcal{M}_2)d\Theta_2}}{B
  = (p(y|M_1)) / p(y|M_2) = (integral of p(y|Theta_1,M_1)p(Theta_1|M_1)dTheta_1)
  / (integral of p(y|Theta_2,M_2)p(Theta_2|M_2)dTheta_2)}

  where \eqn{p(\textbf{y}|\mathcal{M}_1)}{p(y|M_1)} is the marginal
  likelihood of the data in model 1.

  The \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, \code{\link{LaplacesDemon}},
  \code{\link{PMC}}, and \code{\link{VariationalBayes}} functions each
  return the \code{\link{LML}}, the approximate logarithm of the
  marginal likelihood of the data, in each fitted object of class
  \code{iterquad}, \code{laplace}, \code{demonoid}, \code{pmc}, or
  \code{vb}. The \code{BayesFactor} function calculates matrix \code{B},
  a matrix of Bayes factors, where each element of matrix \code{B} is a
  comparison of two models. Each Bayes factor is calculated as the
  exponentiated difference of \code{LML} of model 1
  (\eqn{\mathcal{M}_1}{M[1]}) and \code{LML} of model 2
  (\eqn{\mathcal{M}_2}{M[2]}), and the hypothesis for each element of
  matrix \code{B} is that the model associated with the row is greater
  than the model associated with the column. For example, element
  \code{B[3,2]} is the Bayes factor that model 3 is greater than model
  2. The 'Strength of Evidence' aids in the interpretation (Jeffreys,
  1961).

  A table for the interpretation of the strength of evidence for Bayes
  factors is available at
  \url{https://web.archive.org/web/20150214194051/http://www.bayesian-inference.com/bayesfactors}.

  Each Bayes factor, \code{B}, is the posterior odds in favor of the
  hypothesis divided by the prior odds in favor of the hypothesis, where
  the hypothesis is usually
  \eqn{\mathcal{M}_1 > \mathcal{M}_2}{M[1] > M[2]}. For example, when
  \code{B[3,2]=2}, the data favor \eqn{\mathcal{M}_3}{M[3]} over
  \eqn{\mathcal{M}_2}{M[2]} with 2:1 odds.

  It is also popular to consider the natural logarithm of the Bayes
  factor. The scale of the logged Bayes factor is the same above and
  below one, which is more appropriate for visual comparisons. For
  example, when comparing two Bayes factors at 0.5 and 2, the logarithm
  of these Bayes factors is -0.69 and 0.69.

  Gelman finds Bayes factors generally to be irrelevant, because they
  compute the relative probabilities of the models conditional on one
  of them being true. Gelman prefers approaches that measure the
  distance of the data to each of the approximate models (Gelman et al.,
  2004, p. 180), such as with posterior predictive checks (see the
  \code{\link{predict.iterquad}} function regarding iterative
  quadrature, \code{\link{predict.laplace}} function in the context of
  Laplace Approximation, \code{\link{predict.demonoid}} function in the
  context of MCMC, \code{\link{predict.pmc}} function in the context
  of PMC, or \code{\link{predict.vb}} function in the context of
  Variational Bayes). Kass et al. (1995) asserts this can be done
  without assuming one model is the true model.
}
\value{
  \code{BayesFactor} returns an object of class \code{bayesfactor} that
  is a list with the following components:
  \item{B}{This is a matrix of Bayes factors.}
  \item{Hypothesis}{
    This is the hypothesis, and is stated as 'row > column', indicating
    that the model associated with the row of an element in matrix
    \code{B} is greater than the model associated with the column of
    that element.}
  \item{Strength.of.Evidence}{
    This is the strength of evidence in favor of the hypothesis.}
  \item{Posterior.Probability}{
    This is a vector of the posterior probability of each model, given
    flat priors.}
}
\references{
  Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
  Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
  Hall, London.

  Jeffreys, H. (1961). "Theory of Probability, Third Edition". Oxford
  University Press: Oxford, England.

  Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". \emph{Journal
  of the American Statistical Association}, 90(430), p. 773--795.

  Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
  Posterior Simulation with the Laplace-Metropolis Estimator".
  \emph{Journal of the American Statistical Association}, 92,
  p. 648--655.
}
\author{Statisticat, LLC.}
\seealso{
  \code{\link{is.bayesfactor}},
  \code{\link{is.proper}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LML}},
  \code{\link{PMC}},
  \code{\link{predict.demonoid}},
  \code{\link{predict.iterquad}},
  \code{\link{predict.laplace}},
  \code{\link{predict.pmc}},
  \code{\link{predict.vb}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
# The following example fits a model as Fit1, then adds a predictor, and
# fits another model, Fit2. The two models are compared with Bayes
# factors.

library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
J <- 2
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,10]+1)))
X[,2] <- CenterScale(X[,2])

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

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

########################  Laplace Approximation  ##########################
Fit1 <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
     Iterations=10000)
Fit1

##############################  Demon Data  ###############################
data(demonsnacks)
J <- 3
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(7,8)]))
X[,2] <- CenterScale(X[,2])
X[,3] <- CenterScale(X[,3])

#########################  Data List Preparation  #########################
mon.names <- c("sigma","mu[1]")
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) return(c(rnormv(Data$J,0,10), rhalfcauchy(1,5)))
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

########################  Laplace Approximation  ##########################
Fit2 <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
     Iterations=10000)
Fit2

#############################  Bayes Factor  ##############################
Model.list <- list(M1=Fit1, M2=Fit2)
BayesFactor(Model.list)
}
\keyword{Bayes Factor}
\keyword{Hypothesis Testing}
\keyword{Model Selection}
