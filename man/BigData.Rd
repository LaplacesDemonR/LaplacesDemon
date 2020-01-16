\name{BigData}
\alias{BigData}
\title{Big Data}
\description{
  This function enables Bayesian inference with data that is too large
  for computer memory (RAM) with the simplest method: reading in batches
  of data (where each batch is a section of rows), applying a function
  to the batch, and combining the results.
}
\usage{
BigData(file, nrow, ncol, size=1, Method="add", CPUs=1, Type="PSOCK",
FUN, ...)
}
\arguments{
  \item{file}{This required argument accepts a path and filename that
    must refer to a .csv file, and that must contain only a numeric
    matrix without a header, row names, or column names.}
  \item{nrow}{This required argument accepts a scalar integer that
    indicates the number of rows in the big data matrix.}
  \item{ncol}{This required argument accepts a scalar integer that
    indicates the number of columns in the big data matrix.}
  \item{size}{This argument accepts a scalar integer that specifies the
    number of rows of each batch. The last batch is not required to have
    the same number of rows as the other batches. The largest possible
    size, and therefore the fewest number of batches, should be
    preferred.}
  \item{Method}{This argument accepts a scalar string, defaults to
    "add", and alternatively accepts "rbind". When
    \code{Method="rbind"}, the user-specified function \code{FUN} is
    applied to each batch, and results are combined together by rows.
    For example, if calculating \eqn{\mu = \textbf{X}\beta}{mu = XB} in,
    say, 10 batches, then the output column vector \eqn{\mu}{mu} is equal
    to the number of rows of the big data set.}
  \item{CPUs}{This argument accepts an integer that specifies the number
    of central processing units (CPUs) of the multicore computer or
    computer cluster. This argument defaults to \code{CPUs=1}, in which
    parallel processing does not occur.}
  \item{Type}{This argument specifies the type of parallel processing to
    perform, accepting either \code{Type="PSOCK"} or
    \code{Type="MPI"}.}
  \item{FUN}{This required argument accepts a user-specified function
    that will be performed on each batch. The first argument in the
    function must be the data.}
  \item{...}{Additional arguments are used within the user-specified
    function. Additional arguments often refer to parameters.}
}
\details{
  Big data is defined loosely here as data that is too large for
  computer memory (RAM). The \code{BigData} function uses the
  split-apply-combine strategy with a big data set. The unmanageable
  big data set is split into smaller, manageable pieces (batches),
  a function is applied to each batch, and results are combined.

  Each iteration, the \code{BigData} function opens a connection to a
  big data set and keeps the connection open while the \code{scan}
  function reads in each batch of data (elsewhere, batches are often
  referred to chunks). A user-specified function is applied to each
  batch of data, the results are combined together, the connection is
  closed, and the results are returned.

  As an introductory example, suppose a statistician updates a linear
  regression model, but the design matrix \eqn{\textbf{X}}{X} is too
  large for computer memory. Suppose the design matrix has 100 million
  rows, and the statistician specifies \code{size=1e6}. The statistician
  combines dependent variable \eqn{\textbf{y}} with design matrix
  \eqn{\textbf{X}}. Each iteration in \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, \code{\link{LaplacesDemon}},
  \code{\link{PMC}}, or \code{\link{VariationalBayes}}, the
  \code{BigData} function sequentially reads in one million rows of the
  combined data \eqn{\textbf{X}}, calculates expectation vector
  \eqn{\mu}, and finally returns the sum of the log-likelihood. The sum
  of the log-likelihood is added together for all batches, and returned.
  
  There are many limitations with this function.

  This function is not fast, in the sense that the entire big data set
  is processed in batches, each iteration. With iterative methods, this
  may perform well, albeit slowly.

  There are many functions that cannot be performed on batches, though
  most models in the Examples vignette may easily be updated with big
  data.

  Large matrices of samples are unaddressed, only the data.

  Although many (but not all) models may be estimated, many additional
  functions in this package will not work when applied after the model
  has updated. Instead, a batch or random sample of data (see the
  \code{\link{read.matrix}} function for sampling from big data) should
  be used in the usual way, in the \code{Data} argument, and the
  \code{Model} function coded in the usual way without the
  \code{BigData} function.

  Parallel processing may be performed when the user specifies
  \code{CPUs} to be greater than one, implying that the specified number
  of CPUs exists and is available. Parallelization may be performed on a
  multicore computer or a computer cluster. Either a Simple Network of
  Workstations (SNOW) or Message Passing Interface (MPI) is used. Each
  call to \code{BigData} establishes and closes the parallelization,
  which is costly, and unfortunately results in copious output to the
  console. With small data sets, parallel processing may be slower, due
  to computer network communication. With larger data sets, the user
  should experience a faster run-time.

  There have been several alternative approaches suggested for big data.

  Huang and Gelman (2005) propose that the user creates batches by
  sampling from big data, updating a separate Bayesian model on each
  batch, and combining the results into a consensus posterior. This
  many-mini-model approach may be faster when feasible, because multiple
  models may be updated in parallel, say one per CPU. Such results will
  work with all functions in this package. With the many-mini-model
  approach, several methods are proposed for combining posterior samples
  from batch-level models, such as by using a normal approximation,
  updating from prior to posterior sequentially (the posterior from the
  last batch becomes the prior of the next batch), sample from the full
  posterior via importance sampling from the batched posteriors, and
  more.

  Scott et al. (2013) propose a method that they call Consensus Monte
  Carlo, which consists of breaking the data down into chunks, calling
  each chunk a shard, and use a many-mini-model approach as well, but
  propose their own method of weighting the posteriors back together.
  
  Balakrishnan and Madigan (2006) introduced a Sequential Monte Carlo
  (SMC) sampler, a refinement of an earlier proposal, that was designed
  for big data. It makes one pass through the massive data set, after an
  initial MCMC estimation on a small sample. Each particle is updated
  for each record, resulting in numerous evaluations per record.

  Welling and Teh (2011) proposed a new class of MCMC sampler in which
  only a random sample of big data is used each iteration. The
  stochastic gradient Langevin dynamics (SGLD) algorithm is available
  in the \code{\link{LaplacesDemon}} function.

  An important alternative to consider is using the \code{ff} package,
  where "ff" stands for fast access file. The \code{ff} package has been
  tested successfully with updating a model in \code{LaplacesDemon}.
  Once the big data set, say \eqn{\textbf{X}}{X}, is an object of
  class \code{ff_matrix}, simply include it in the list of data as
  usual, and modify the \code{Model} specification function
  appropriately. For example, change \code{mu <- tcrossprod(X, t(beta))}
  to \code{mu <- tcrossprod(X[], t(beta))}. The \code{ff} package is
  not included as a dependency in the \code{LaplacesDemon} package, so
  it must be installed and activated.
}
\value{
  The \code{BigData} function returns output that is the result of
  performing a user-specified function on batches of big data. Output is
  a matrix, and may have one or more column vectors.
}
\references{
  Balakrishnan, S. and Madigan, D. (2006). "A One-Pass Sequential Monte
  Carlo Method for Bayesian Analysis of Massive Datasets".
  \emph{Bayesian Analysis}, 1(2), p. 345--362.

  Huang, Z. and Gelman, A. (2005) "Sampling for Bayesian Computation
  with Large Datasets". \emph{SSRN eLibrary}.

  Scott, S.L., Blocker, A.W. and Bonassi, F.V. (2013). "Bayes and Big
  Data: The Consensus Monte Carlo Algorithm". In \emph{Bayes 250}.

  Welling, M. and Teh, Y.W. (2011). "Bayesian Learning via Stochastic
  Gradient Langevin Dynamics". \emph{Proceedings of the 28th
  International Conference on Machine Learning (ICML)}, p. 681--688.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.RAM}},
  \code{\link{PMC}},
  \code{\link{PMC.RAM}},
  \code{\link{read.matrix}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
### Below is an example of a linear regression model specification
### function in which BigData reads in a batch of 1,000 records of
### Data$N records from a data set that is too large to fully open
### in memory. The example simulates on 10,000 records, which is
### not big data; it's just a toy example. The data set is file X.csv,
### and the first column of matrix X is the dependent variable y. The
### user supplies a function to BigData along with parameters beta and
### sigma. When each batch of 1,000 records is read in,
### mu = XB is calculated, and then the LL is calculated as
### y ~ N(mu, sigma^2). These results are added together from all
### batches, and returned as LL.

library(LaplacesDemon)
N <- 10000
J <- 10 #Number of predictors, including the intercept
X <- matrix(1,N,J)
for (j in 2:J) {X[,j] <- rnorm(N,runif(1,-3,3),runif(1,0.1,1))}
beta.orig <- runif(J,-3,3)
e <- rnorm(N,0,0.1)
y <- as.vector(tcrossprod(beta.orig, X) + e)
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), log.sigma=0))
PGF <- function(Data) return(c(rnormv(Data$J,0,0.01),
     log(rhalfcauchy(1,1))))
MyData <- list(J=J, PGF=PGF, N=N, mon.names=mon.names,
     parm.names=parm.names) #Notice that X and y are not included here
filename <- tempfile("X.csv")  
write.table(cbind(y,X), filename, sep=",", row.names=FALSE,
  col.names=FALSE)

Model <- function(parm, Data)
     {
     ### Parameters
     beta <- parm[1:Data$J]
     sigma <- exp(parm[Data$J+1])
     ### Log(Prior Densities)
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     ### Log-Likelihood
     LL <- BigData(file=filename, nrow=Data$N, ncol=Data$J+1, size=1000,
          Method="add", CPUs=1, Type="PSOCK",
          FUN=function(x, beta, sigma) sum(dnorm(x[,1], tcrossprod(x[,-1],
               t(beta)), sigma, log=TRUE)), beta, sigma)
     ### Log-Posterior
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma),
          yhat=0,#rnorm(length(mu), mu, sigma),
          parm=parm)
     return(Modelout)
     }

### From here, the user may update the model as usual.
}
\keyword{High Performance Computing}
\keyword{Utility}
