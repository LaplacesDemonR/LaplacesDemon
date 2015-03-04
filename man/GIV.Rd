\name{GIV}
\alias{GIV}
\title{Generate Initial Values}
\description{
  The \code{GIV} function generates initial values for use with the
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, \code{\link{PMC}}, and
  \code{\link{VariationalBayes}} functions.
}
\usage{
GIV(Model, Data, n=1000, PGF=FALSE)
}
\arguments{
  \item{Model}{This required argument is a model specification
    function. For more information, see \code{\link{LaplacesDemon}}.}
  \item{Data}{This required argument is a list of data. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{n}{This is the number of attempts to generate acceptable
    initial values.}
  \item{PGF}{Logical. When \code{TRUE}, a Parameter-Generating Function
    (PGF) is required to be in \code{Data}, and \code{GIV} will generate
    initial values according to the user-specified PGF. This argument
    defaults to \code{FALSE}, in which case initial values are generated
    randomly without respect to a user-specified function.}
}
\details{
  Initial values are required for optimization or sampling algorithms. A
  user may specify initial values, or use the \code{GIV} function for
  random generation. Initial values determined by the user may fail to
  produce a finite posterior in complicated models, and the \code{GIV}
  function is here to help.

  \code{GIV} has several uses. First, the
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, and \code{\link{VariationalBayes}}
  functions use \code{GIV} internally if unacceptable initial values are
  discovered. Second, the user may use \code{GIV} when developing their
  model specification function, \code{Model}, to check for potential
  problems. Third, the user may prefer to randomly generate acceptable
  initial values. Lastly, \code{GIV} is recommended when running
  multiple or parallel chains with the \code{\link{LaplacesDemon.hpc}}
  function (such as for later use with the \code{Gelman.Diagnostic}) for
  dispersed starting locations. For dispersed starting locations,
  \code{GIV} should be run once for each parallel chain, and the results
  should be stored per row in a matrix of initial values. For more
  information, see the \code{LaplacesDemon.hpc} documentation for
  initial values.

  It is strongly recommended that the user specifies a
  Parameter-Generating Function (PGF), and includes this function in the
  list of data. Although the PGF may be specified according to the prior
  distributions (possibly considered as a Prior-Generating Function), it
  is often specified with a more restricted range. For example, if a
  user has a model with the following prior distributions

  \deqn{\beta_j \sim \mathcal{N}(0, 1000), j=1,\dots,5}{beta_j ~ N(0,
    1000), j=1,\dots,5}
  \deqn{\sigma \sim \mathcal{HC}(25)}{sigma ~ HC(25)}

  then the PGF, given the prior distributions, is

  \code{PGF <- function(Data) return(c(rnormv(5,0,1000),rhalfcauchy(1,25)))}

  However, the user may not want to begin with initial values that could
  be so far from zero (as determined by the variance of 1000), and may
  instead prefer

  \code{PGF <- function(Data) return(c(rnormv(5,0,10),rhalfcauchy(1,5)))}

  When \code{PGF=FALSE}, initial values are attempted to be constrained
  to the interval \eqn{[-100,100]}. This is done to prevent numeric
  overflows with parameters that are exponentiated within the model
  specification function. First, \code{GIV} passes the upper and lower
  bounds of this interval to the model, and any changed parameters are
  noted.

  At this point, it is hoped that a non-finite posterior is not
  found. If found, then the remainder of the process is random and
  without the previous bounds. This can be particularly problematic in
  the case of, say, initial values that are the elements of a matrix
  that must be positive-definite, especially with large matrices. If a
  random solution is not found, then \code{GIV} will fail.
  
  If the posterior is finite and \code{PGF=FALSE}, then initial values
  are randomly generated with a normal proposal and a small variance at
  the center of the returned range of each parameter. As \code{GIV}
  fails to find acceptable initial values, the algorithm iterates toward
  its maximum number of iterations, \code{n}. In each iteration, the
  variance increases for the proposal.
  
  Initial values are considered acceptable only when the first two
  returned components of \code{Model} (which are \code{LP} and
  \code{Dev}) are finite, and when initial values do not change through
  constraints, as returned in the fifth component of the list:
  \code{parm}.

  If \code{GIV} fails to return acceptable initial values, then it is
  best to study the model specification function. When the model is
  complicated, here is a suggestion. Remove the log-likelihood,
  \code{LL}, from the equation that calculates the logarithm of the
  unnormalized joint posterior density, \code{LP}. For example, convert
  \code{LP <- LL + beta.prior} to \code{LP <- beta.prior}. Now, maximize
  \code{LP}, which is merely the set of prior densities, with any
  optimization algorithm. Replace \code{LL}, and run the model with
  initial values that are in regions of high prior density (preferably
  with \code{PGF=TRUE}. If this fails, then the model specification
  should be studied closely, because a non-finite posterior should
  (especially) never be associated with regions of high prior density.
}
\value{
  The \code{GIV} function returns a vector equal in length to the
  number of parameters, and each element is an initial value for the
  associated parameter in \code{Data$parm.names}. When \code{GIV} fails
  to find acceptable initial values, each returned element is \code{NA}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{as.initial.values}},
  \code{\link{Gelman.Diagnostic}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

#########################  Data List Preparation  #########################
mon.names <- c("LP","sigma")
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

########################  Generate Initial Values  ########################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)
}
\keyword{Initial Values}
\keyword{Utility}
