\name{VariationalBayes}
\alias{VariationalBayes}
\title{Variational Bayes}
\description{
  The \code{VariationalBayes} function is a numerical approximation
  method for deterministically estimating the marginal posterior
  distributions, target distributions, in a Bayesian model with
  approximated distributions by minimizing the Kullback-Leibler
  Divergence (\code{\link{KLD}}) between the target and its
  approximation.
}
\usage{
VariationalBayes(Model, parm, Data, Covar=NULL, Interval=1.0E-6,
     Iterations=1000, Method="Salimans2", Samples=1000, sir=TRUE,
     Stop.Tolerance=1.0E-5, CPUs=1, Type="PSOCK")
}
\arguments{
  \item{Model}{This required argument receives the model from a
    user-defined function. The user-defined function is where the model
    is specified. \code{VariationalBayes} passes two arguments to
    the model function, \code{parms} and \code{Data}. For more
    information, see the \code{\link{LaplacesDemon}} function and
    ``LaplacesDemon Tutorial'' vignette.}
  \item{parm}{This argument requires a vector of initial values equal in
    length to the number of parameters. \code{VariationalBayes} will
    attempt to optimize these initial values for the parameters, where
    the optimized values are the posterior means, for later use with the
    \code{\link{IterativeQuadrature}}, \code{\link{LaplacesDemon}}, or
    \code{\link{PMC}} function. The \code{\link{GIV}} function may be
    used to randomly generate initial values. Parameters must be
    continuous.}
  \item{Data}{This required argument accepts a list of data. The list of
    data must include \code{mon.names} which contains monitored variable
    names, and \code{parm.names} which contains parameter
    names. \code{VariationalBayes} must be able to determine the
    sample size of the data, and will look for a scalar sample size
    variable \code{n} or \code{N}. If not found, it will look for
    variable \code{y} or \code{Y}, and attempt to take its number of
    rows as sample size. \code{VariationalBayes} needs to determine
    sample size due to the asymptotic nature of this method. Sample size
    should be at least \eqn{\sqrt{J}}{sqrt(J)} with \eqn{J} exchangeable
    parameters.}
  \item{Covar}{This argument defaults to \code{NULL}, but may otherwise
    accept a \eqn{K \times K}{K x K} covariance matrix (where \eqn{K}
    is the number of dimensions or parameters) of the parameters. When
    the model is updated for the first time and prior variance or
    covariance is unknown, then \code{Covar=NULL} should be used. Once
    \code{VariationalBayes} has finished updating, it may be desired
    to continue updating where it left off, in which case the covariance
    matrix from the last run can be input into the next run.}
  \item{Interval}{This argument receives an interval for estimating
    approximate gradients. The logarithm of the unnormalized joint
    posterior density of the Bayesian model is evaluated at the current
    parameter value, and again at the current parameter value plus this
    interval.}
  \item{Iterations}{This argument accepts an integer that determines the
    number of iterations that \code{VariationalBayes} will attempt
    to maximize the logarithm of the unnormalized joint posterior
    density. \code{Iterations} defaults to 1000.
    \code{VariationalBayes} will stop before this number of
    iterations if the tolerance is less than or equal to the
    \code{Stop.Tolerance} criterion. The required amount of computer
    memory increases with \code{Iterations}. If computer memory is
    exceeded, then all will be lost.}
  \item{Method}{This optional argument currently accepts only
    \code{Salimans2}, which is the second algorithm in Salimans and
    Knowles (2013).}
  \item{Samples}{This argument indicates the number of posterior samples
    to be taken with sampling importance resampling via the
    \code{\link{SIR}} function, which occurs only when
    \code{sir=TRUE}. Note that the number of samples should increase
    with the number and intercorrelations of the parameters.}
  \item{sir}{This logical argument indicates whether or not Sampling
    Importance Resampling (SIR) is conducted via the \code{\link{SIR}}
    function to draw independent posterior samples. This argument
    defaults to \code{TRUE}. Even when \code{TRUE}, posterior samples
    are drawn only when \code{VariationalBayes} has
    converged. Posterior samples are required for many other functions,
    including \code{plot.vb} and \code{predict.vb}. The only
    time that it is advantageous for \code{sir=FALSE} is when
    \code{VariationalBayes} is used to help the initial values for
    \code{\link{IterativeQuadrature}}, \code{\link{LaplacesDemon}}, or
    \code{\link{PMC}}, and it is unnecessary for time to be spent on
    sampling. Less time can be spent on sampling by increasing
    \code{CPUs}, which parallelizes the sampling.}
  \item{Stop.Tolerance}{This argument accepts any positive number and
    defaults to 1.0E-3. Tolerance is calculated each iteration, and the
    criteria varies by algorithm. The algorithm is considered to have
    converged to the user-specified \code{Stop.Tolerance} when the
    tolerance is less than or equal to the value of
    \code{Stop.Tolerance}, and the algorithm terminates at the end of
    the current iteration. Often, multiple criteria are used, in
    which case the maximum of all criteria becomes the tolerance. For
    example, when partial derivatives are taken, it is commonly required
    that the Euclidean norm of the partial derivatives is a criterion,
    and another common criterion is the Euclidean norm of the
    differences between the current and previous parameter
    values. Several algorithms have other, specific tolerances.}
  \item{CPUs}{This argument accepts an integer that specifies the number
    of central processing units (CPUs) of the multicore computer or
    computer cluster. This argument defaults to \code{CPUs=1}, in which
    parallel processing does not occur. Parallelization occurs only for
    sampling with \code{\link{SIR}} when \code{sir=TRUE}.}
  \item{Type}{This argument specifies the type of parallel processing to
    perform, accepting either \code{Type="PSOCK"} or
    \code{Type="MPI"}.}
}
\details{
  Variational Bayes (VB) is a family of numerical approximation
  algorithms that is a subset of variational inference algorithms, or
  variational methods. Some examples of variational methods include the
  mean-field approximation, loopy belief propagation, tree-reweighted
  belief propagation, and expectation propagation (EP).

  Variational inference for probabilistic models was introduced in the
  field of machine learning, influenced by statistical physics
  literature (Saul et al., 1996; Saul and Jordan, 1996; Jaakkola, 1997).
  The mean-field methods in Neal and Hinton (1999) led to variational
  algorithms.
  
  Variational inference algorithms were later generalized for conjugate
  exponential-family models (Attias, 1999, 2000; Wiegerinck, 2000;
  Ghahramani and Beal, 2001; Xing et al., 2003). These algorithms still
  require different designs for different model forms. Salimans and
  Knowles (2013) introduced general-purpose VB algorithms for Gaussian
  posteriors.
  
  A VB algorithm deterministically estimates the marginal posterior
  distributions (target distributions) in a Bayesian model with
  approximated distributions by minimizing the Kullback-Leibler
  Divergence (\code{\link{KLD}}) between the target and its
  approximation. The complicated posterior distribution is approximated
  with a simpler distribution. The simpler, approximated distribution is
  called the variational approximation, or approximation distribution,
  of the posterior. The term variational is derived from the calculus of
  variations, and regards optimization algorithms that select the best
  function (which is a distribution in VB), rather than merely selecting
  the best parameters.

  VB algorithms often use Gaussian distributions as approximating
  distributions. In this case, both the mean and variance of the
  parameters are estimated.
  
  Usually, a VB algorithm is slower to convergence than a Laplace
  Approximation algorithm, and faster to convergence than a Monte Carlo
  algorithm such as Markov chain Monte Carlo (MCMC). VB often provides
  solutions with comparable accuracy to MCMC in less time. Though Monte
  Carlo algorithms provide a numerical approximation to the exact
  posterior using a set of samples, VB provides a locally-optimal,
  exact analytical solution to an approximation of the posterior. VB is
  often more applicable than MCMC to big data or large-dimensional
  models.

  Since VB is deterministic, it is asymptotic and subject to the same
  limitations with respect to sample size as Laplace Approximation.
  However, VB estimates more parameters than Laplace Approximation,
  such as when Laplace Approximation optimizes the posterior mode of a
  Gaussian distribution, while VB optimizes both the Gaussian mean and
  variance.

  Traditionally, VB algorithms required customized equations. The
  \code{VariationalBayes} function uses general-purpose algorithms. A
  general-purpose VB algorithm is less efficient than an algorithm
  custom designed for the model form. However, a general-purpose
  algorithm is applied consistently and easily to numerous model forms.

  When \code{Method="Salimans2"}, the second algorithm of Salimans and
  Knowles (2013) is used. This requires the gradient and Hessian, which
  is more efficient with a small number of parameters as long as the
  posterior is twice differentiable. The step size is constant. This
  algorithm is suitable for marginal posterior distributions that are
  Gaussian and unimodal. A stochastic approximation algorithm is used
  in the context of fixed-form VB, inspired by considering fixed-form VB
  to be equivalent to performing a linear regression with the sufficient
  statistics of the approximation as independent variables and the
  unnormalized logarithm of the joint posterior density as the dependent
  variable. The number of requested iterations should be large, since the
  step-size decreases for larger requested iterations, and a small
  step-size will eventually converge. A large number of requested
  iterations results in a smaller step-size and better convergence
  properties, so hope for early convergence. However convergence is
  checked only in the last half of the iterations after the algorithm
  begins to average the mean and variance from the samples of the
  stochastic approximation. The history of stochastic samples is
  returned.
}
\value{
  \code{VariationalBayes} returns an object of class \code{vb}
  that is a list with the following components:
  \item{Call}{This is the matched call of \code{VariationalBayes}.}
  \item{Converged}{This is a logical indicator of whether or not
    \code{VariationalBayes} converged within the specified
    \code{Iterations} according to the supplied \code{Stop.Tolerance}
    criterion. Convergence does not indicate that the global maximum has
    been found, but only that the tolerance was less than or equal to
    the \code{Stop.Tolerance} criterion.}
  \item{Covar}{This is the estimated covariance matrix. The
    \code{Covar} matrix may be scaled and input into the \code{Covar}
    argument of the \code{\link{LaplacesDemon}} or \code{\link{PMC}}
    function for further estimation, or the diagonal of this matrix may
    be used to represent the posterior variance of the parameters,
    provided the algorithm converged and matrix inversion was
    successful. To scale this matrix for use with Laplace's Demon or
    PMC, multiply it by \eqn{2.38^2/d}, where \eqn{d} is the number
    of initial values.}
  \item{Deviance}{This is a vector of the iterative history of the
    deviance in the \code{VariationalBayes} function, as it sought
    convergence.}
  \item{History}{This is an array of the iterative history of the
    parameters in the \code{VariationalBayes} function, as it sought
    convergence. The first matrix is for means and the second matrix is
    for variances.}
  \item{Initial.Values}{This is the vector of initial values that was
    originally given to \code{VariationalBayes} in the \code{parm}
    argument.}
  \item{LML}{This is an approximation of the logarithm of the marginal
    likelihood of the data (see the \code{\link{LML}} function for more
    information). When the model has converged and \code{sir=TRUE}, the
    NSIS method is used. When the model has converged and
    \code{sir=FALSE}, the LME method is used. This is the
    logarithmic form of equation 4 in Lewis and Raftery (1997). As a
    rough estimate of Kass and Raftery (1995), the LME-based LML is
    worrisome when the sample size of the data is less than five times
    the number of parameters, and \code{LML} should be adequate in most
    problems when the sample size of the data exceeds twenty times the
    number of parameters (p. 778). The LME is inappropriate with
    hierarchical models. However \code{LML} is estimated, it is useful
    for comparing multiple models with the \code{BayesFactor} function.}
  \item{LP.Final}{This reports the final scalar value for the logarithm
    of the unnormalized joint posterior density.}
  \item{LP.Initial}{This reports the initial scalar value for the
    logarithm of the unnormalized joint posterior density.}
  \item{Minutes}{This is the number of minutes that
    \code{VariationalBayes} was running, and this includes the
    initial checks as well as drawing posterior samples and creating
    summaries.}
  \item{Monitor}{When \code{sir=TRUE}, a number of independent
    posterior samples equal to \code{Samples} is taken, and the draws
    are stored here as a matrix. The rows of the matrix are the samples,
    and the columns are the monitored variables.}
  \item{Posterior}{When \code{sir=TRUE}, a number of independent
    posterior samples equal to \code{Samples} is taken, and the draws
    are stored here as a matrix. The rows of the matrix are the samples,
    and the columns are the parameters.}
  \item{Step.Size.Final}{This is the final, scalar \code{Step.Size}
    value at the end of the \code{VariationalBayes} algorithm.}
  \item{Step.Size.Initial}{This is the initial, scalar \code{Step.Size}.}
  \item{Summary1}{This is a summary matrix that summarizes the
    point-estimated posterior means and variances. Uncertainty around
    the posterior means is estimated from the estimated covariance
    matrix. Rows are parameters. The following columns are included:
    Mean, SD (Standard Deviation), LB (Lower Bound), and UB (Upper
    Bound). The bounds constitute a 95\% probability interval.}
  \item{Summary2}{This is a summary matrix that summarizes the
    posterior samples drawn with sampling importance resampling
    (\code{\link{SIR}}) when \code{sir=TRUE}, given the point-estimated
    posterior means and covariance matrix. Rows are parameters. The
    following columns are included: Mean, SD (Standard Deviation),
    LB (Lower Bound), and UB (Upper Bound). The bounds constitute a 95\%
    probability interval.}
  \item{Tolerance.Final}{This is the last \code{Tolerance} of the
    \code{VariationalBayes} algorithm.}
  \item{Tolerance.Stop}{This is the \code{Stop.Tolerance} criterion.}
}
\references{
  Attias, H. (1999). "Inferring Parameters and Structure of Latent
  Variable Models by Variational Bayes". In \emph{Uncertainty in
  Artificial Intelligence}.
  
  Attias, H. (2000). "A Variational Bayesian Framework for Graphical
  Models". In \emph{Neural Information Processing Systems}.

  Ghahramani, Z. and Beal, M. (2001). "Propagation Algorithms for
  Variational Bayesian Learning". In \emph{Neural Information Processing
  Systems}, p. 507--513.

  Jaakkola, T. (1997). "Variational Methods for Inference and Estimation
  in Graphical Models". PhD thesis, Massachusetts Institute of
  Technology.

  Salimans, T. and Knowles, D.A. (2013). "Fixed-Form Variational
  Posterior Approximation through Stochastic Linear Regression".
  \emph{Bayesian Analysis}, 8(4), p. 837--882.

  Neal, R. and Hinton, G. (1999). "A View of the EM Algorithm that
  Justifies Incremental, Sparse, and Other Variants". In Learning in
  Graphical Models, p. 355--368. MIT Press, 1999.

  Saul, L. and Jordan, M. (1996). "Exploiting Tractable Substructures in
  Intractable Networks". \emph{Neural Information Processing Systems}.

  Saul, L., Jaakkola, T., and Jordan, M. (1996). "Mean Field Theory for
  Sigmoid Belief Networks". \emph{Journal of Artificial Intelligence
  Research}, 4, p. 61--76.

  Wiegerinck, W. (2000). "Variational Approximations Between Mean Field
  Theory and the Junction Tree Algorithm". In \emph{Uncertainty in
  Artificial Intelligence}.

  Xing, E., Jordan, M., and Russell, S. (2003). "A Generalized Mean
  Field Algorithm for Variational Inference in Exponential Families". In
  \emph{Uncertainty in Artificial Intelligence}.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BayesFactor}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{GIV}},
  \code{\link{LML}},
  \code{\link{PMC}}, and
  \code{\link{SIR}}.
}
\examples{
# The accompanying Examples vignette is a compendium of examples.
####################  Load the LaplacesDemon Library  #####################
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,10]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

#########################  Data List Preparation  #########################
mon.names <- "mu[1]"
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
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=mu[1],
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     return(Modelout)
     }

############################  Initial Values  #############################
#Initial.Values <- GIV(Model, MyData, PGF=TRUE)
Initial.Values <- rep(0,J+1)

#Fit <- VariationalBayes(Model, Initial.Values, Data=MyData, Covar=NULL,
#     Iterations=1000, Method="Salimans2", Stop.Tolerance=1e-3, CPUs=1)
#Fit
#print(Fit)
#PosteriorChecks(Fit)
#caterpillar.plot(Fit, Parms="beta")
#plot(Fit, MyData, PDF=FALSE)
#Pred <- predict(Fit, Model, MyData, CPUs=1)
#summary(Pred, Discrep="Chi-Square")
#plot(Pred, Style="Covariates", Data=MyData)
#plot(Pred, Style="Density", Rows=1:9)
#plot(Pred, Style="Fitted")
#plot(Pred, Style="Jarque-Bera")
#plot(Pred, Style="Predictive Quantiles")
#plot(Pred, Style="Residual Density")
#plot(Pred, Style="Residuals")
#Levene.Test(Pred)
#Importance(Fit, Model, MyData, Discrep="Chi-Square")

#Fit$Covar is scaled (2.38^2/d) and submitted to LaplacesDemon as Covar.
#Fit$Summary[,1] is submitted to LaplacesDemon as Initial.Values.
#End
}
\keyword{Adaptive}
\keyword{Optimization}
