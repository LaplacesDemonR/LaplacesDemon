\name{Model.Specification.Time}
\alias{Model.Spec.Time}
\title{Model Specification Time}
\description{
  The \code{Model.Spec.Time} function returns the time in minutes to
  evaluate a model specification a number of times, as well as
  the evaluations per minute, and componentwise iterations per minute.
}
\usage{
Model.Spec.Time(Model, Initial.Values, Data, n=1000)
}
\arguments{
  \item{Model}{This requried argument is a model specification
    function. For more information, see \code{\link{LaplacesDemon}}.}
  \item{Initial.Values}{This required argument is a vector of initial
    values for the parameters.}
  \item{Data}{This required argument is a list of data. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{n}{This is the number of evaluations of the model specification,
    and accuracy increases with \code{n}.}
}
\details{
  The largest single factor to affect the run-time of an algorithm --
  whether it is with \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, \code{\link{LaplacesDemon}},
  \code{\link{PMC}}, or \code{\link{VariationalBayes}} -- is the time
  that it takes to evaluate the model specification. This has also been
  observed in Rosenthal (2007). It is highly recommended that a user of
  the \code{LaplacesDemon} package should attempt to reduce the run-time
  of the model specification, usually by testing alternate forms of code
  for speed. This is especially true with big data, such as with the
  \code{\link{BigData}} function.

  Every function in the LaplacesDemon package is byte-compiled, which is
  a recent option in R. This reduces run-time, thanks to Tierney's
  compiler package in base R. The model specification, however, is
  specified by the user, and should be byte-compiled. The reduction in
  run-time may range from mild to dramatic, depending on the model. It
  is highly recommended that users concerned with run-time activate the
  compiler package and use the \code{cmpfun} function, as per the
  example below.
  
  A model specification function that is optimized for speed and
  involves many records may result in a model update run-time that is
  considerably less than other popular MCMC-based software algorithms
  that loop through records, even when those algorithms are coded in
  \code{C} or other fast languages. For a comparison, see the
  ``Laplace's Demon Tutorial'' vignette.

  However, if a model specification function in the LaplacesDemon
  package is not fully vectorized (contains \code{for} loops and
  \code{apply} functions), then run-time will typically be slower than
  other popular MCMC-based software algorithms.

  The speed of calculating the model specification function is
  affected by parameter constraints, such as with the
  \code{\link{interval}} function. Parameter constraints may add
  considerable variability to the speed of this calculation, and usually
  more variation occurs with initial values that are far from the target
  distributions.

  Distributions in the \code{LaplacesDemon} package usually have logical
  checks to ensure correctness. These checks may slow the 
  calculation of the density, for example. If the user is confident
  these checks are unnecessary for their model, then the user may
  copy the code to a new function name and comment-out the checks to
  improve speed.

  When speed is of paramount importance, a user is encouraged to
  experiment with writing the model specification function in another
  language such as in \code{C++} with the \code{Rcpp} package, and
  calling drop-in replacement functions from within the \code{Model}
  function, or re-writing the model function entirely in \code{C++}.
  For an introduction to including \code{C++} in \pkg{LaplacesDemon},
  see \url{https://web.archive.org/web/20150227225556/http://www.bayesian-inference.com/softwarearticlescppsugar}.

  When a model specification function is computationally expensive,
  another approach to reduce run-time may be for the user to write
  parallelized code within the model, splitting up difficult tasks and
  assigning each to a separate CPU.

  Another use for \code{Model.Spec.Time} is to allow the user to make an
  informed decision about which MCMC algorithm to select, given the
  speed of their model specification. For example, the Adaptive
  Metropolis-within-Gibbs (AMWG) of Roberts and Rosenthal (2009) is
  currently the adaptive MCMC algorithm of choice in general in many
  cases, but this choice is conditional on run-time. While other
  MCMC algorithms in \code{LaplacesDemon} evaluate the model
  specification function once per iteration, componentwise algorithms
  such as in the MWG family evaluate it once per parameter per
  iteration, significantly increasing run-time per iteration in large
  models. The \code{Model.Spec.Time} function may forewarn the user of
  the associated run-time, and if it should be decided to go with an
  alternate MCMC algorithm, such as AMM, in which each element of its
  covariance matrix must stabilize for the chains to become
  stationary. AMM, for example, will require many more iterations of 
  burn-in (for more information, see the \code{\link{burnin}} function),
  but with numerous iterations, allows more thinning. A general
  recommendation may be to use AMWG when
  \code{Componentwise.Iters.per.Minute} >= 1000, but this is subjective
  and best determined by each user for each model. A better decision may
  be made by comparing MCMC algorithms with the \code{\link{Juxtapose}}
  function for a particular model.
  
  Following are a few common suggestions for increasing the speed of
  \code{R} code in the model specification function. There are often
  exceptions to these suggestions, so case-by-case experimentation is
  also suggested.

  \itemize{
    \item Replace exponents with multiplied terms, such as \code{x^2}
    with \code{x*x}.
    \item Replace \code{mean(x)} with \code{sum(x)/length(x)}.
    \item Replace parentheses (when possible) with curly brackets, such
    as \code{x*(a+b)} with \code{x*\{a+b\}}.
    \item Replace \code{tcrossprod(Data$X, t(beta))} with
    \code{Data$X \%*\% beta} when there are few predictors, and avoid
    \code{tcrossprod(beta, Data$X)}, which is always slowest.
    \item Vectorize functions and eliminate \code{apply} and \code{for}
    functions. There are often specialized functions. For example,
    \code{max.col(X)} is faster than \code{apply(X, 1, which.max)}.
  }

  When seeking speed, things to consider beyond the LaplacesDemon
  package are the Basic Linear Algebra System (BLAS) and hardware. The
  ATLAS (Automatically Tuned Linear Algebra System) should be worth
  installing (and there are several alternatives), but this discussion
  is beyond the scope of this documentation. Lastly, the speed at which
  the computer can process iterations is limited by its hardware, and
  more should be considered than merely the CPU (Central Processing
  Unit). Again, though, this is beyond the scope of this documentation.
}
\value{
  The \code{Model.Spec.Time} function returns a list with three
  components:

  \item{Time}{This is the time in minutes to evaluate the model
    specification \code{n} times.}
  \item{Evals.per.Minute}{This is the number of evaluations of the
    model specification per minute: \code{n} / \code{Time}. This is also
    the number of iterations per minute in an algorithm that is not
    componentwise, where one evaluation occurs per iteration.}
  \item{Componentwise.Iters.per.Minute}{This is the number of iterations
    per minute in a componentwise MCMC algorithm, such as AMWG or
    MWG. It is the evaluations per minute divided by the number of
    parameters, since an evaluation must occur for each parameter, for
    each iteration.}
}
\references{
  Roberts, G.O. and Rosenthal, J.S. (2009). "Examples of Adaptive
  MCMC". \emph{Computational Statistics and Data Analysis}, 18,
  p. 349--367.
}
\author{Statisticat, LLC.}
\seealso{
  \code{\link{.C}},
  \code{\link{.Fortran}},

  \code{\link{apply}},
  \code{\link{BigData}},
  \code{\link{interval}},
  \code{\link{IterativeQuadrature}},
  \code{\link{Juxtapose}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{max.col}},
  \code{\link{PMC}},
  \code{\link{system.time}}, and
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
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}

#########################  Data List Preparation  #########################
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) return(c(rnormv(Data$J,0,10), rhalfcauchy(1,5)))
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

##########################  Model Specification  ##########################
Model <- function(parm, Data)
     {
     ### Parameters
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
     parm[Data$pos.sigma] <- sigma
     ### Log of Prior Densities
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

############################  Model.Spec.Time  ############################
### Not byte-compiled
Model.Spec.Time(Model, Initial.Values, MyData)
### Byte-compiled
library(compiler)
Model <- cmpfun(Model)
Model.Spec.Time(Model, Initial.Values, MyData)
}
\keyword{Utility}
