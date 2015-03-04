\name{Blocks}
\alias{Blocks}
\title{Blocks}
\description{
  The \code{Blocks} function returns a list of \eqn{N} blocks of
  parameters, for use with some MCMC algorithms in the
  \code{\link{LaplacesDemon}} function. Blocks may be created either
  sequentially, or from a hierarchical clustering of the posterior
  correlation matrix.
}
\usage{
Blocks(Initial.Values, N, PostCor=NULL)
}
\arguments{
  \item{Initial.Values}{This required argument is a vector of initial
    values.}
  \item{N}{This optional argument indicates the desired number of
    blocks. If omitted, then the truncated square root of the number of
    initial values is used. If a posterior correlation matrix is
    supplied to \code{PostCor}, then \code{N} may be a scalar, or
    have length two. If \code{N} has length two, then the first element
    indicates the minimum number of blocks, and the second element
    indicates the maximum number of blocks, and the number of blocks
    is the maximum of the mean silhouette width for each hierarchical
    cluster solution.}
  \item{PostCor}{This optional argument defaults to \code{NULL}, in
    which case sequential blocking is performed. If a posterior
    correlation matrix is supplied, then blocks are created based on
    hierarchical clustering.}
}
\details{
  Usually, there is more than one target distribution in MCMC, in which
  case it must be determined whether it is best to sample from target
  distributions individually, in groups, or all at once. Blockwise
  sampling (also called block updating) refers to splitting a
  multivariate vector into groups called blocks, and each block is
  sampled separately. A block may contain one or more parameters.

  Parameters are usually grouped into blocks such that parameters within
  a block are as correlated as possible, and parameters between blocks
  are as independent as possible. This strategy retains as much of the
  parameter correlation as possible for blockwise sampling, as opposed
  to componentwise sampling where parameter correlation is ignored.
  The \code{\link{PosteriorChecks}} function can be used on the output
  of previous runs to find highly correlated parameters. See examples
  below.

  Advantages of blockwise sampling are that a different MCMC algorithm
  may be used for each block (or parameter, for that matter), creating a
  more specialized approach (though different algorithms by block are
  not supported here), the acceptance of a newly proposed state is
  likely to be higher than sampling from all target distributions at
  once in high dimensions, and large proposal covariance matrices can
  be reduced in size, which is most helpful again in high dimensions.

  Disadvantages of blockwise sampling are that correlations probably
  exist between parameters between blocks, and each block is updated
  while holding the other blocks constant, ignoring these correlations
  of parameters between blocks. Without simultaneously taking
  everything into account, the algorithm may converge slowly or never
  arrive at the proper solution. However, there are instances when it
  may be best when everything is not taken into account at once, such
  as in state-space models. Also, as the number of blocks increases,
  more computation is required, which slows the algorithm. In general,
  blockwise sampling allows a more specialized approach at the expense
  of accuracy, generalization, and speed. Blockwise sampling is offered
  in the following algorithms: Adaptive-Mixture Metropolis (AMM),
  Adaptive Metropolis-within-Gibbs (AMWG), Automated Factor Slice
  Sampler (AFSS), Elliptical Slice Sampler (ESS), Hit-And-Run Metropolis
  (HARM), Metropolis-within-Gibbs (MWG), Random-Walk Metropolis (RWM),
  Robust Adaptive Metropolis (RAM), Slice Sampler (Slice), and the
  Univariate Eigenvector Slice Sampler (UESS).

  Large-dimensional models often require blockwise sampling. For
  example, with thousands of parameters, a componentwise algorithm
  must evaluate the model specification function once per parameter per
  iteration, resulting in an algorithm that may take longer than is
  acceptable to produce samples. Algorithms that require derivatives,
  such as the family of Hamiltonian Monte Carlo (HMC), require even more
  evaluations of the model specification function per iteration, and
  quickly become too costly in large dimensions. Finally, algorithms
  with multivariate proposals often have difficulty producing an
  accepted proposal in large-dimensional models. The most practical
  solution is to group parameters into \eqn{N} blocks, and each
  iteration the algorithm evaluates the model specification function
  \eqn{N} times, each with a reduced set of parameters.

  The \code{Blocks} function performs either a sequential assignment of
  parameters to blocks when posterior correlation is not supplied, or
  uses hierarchical clustering to create blocks based on posterior
  correlation. If posterior correlation is supplied, then the user may
  specify a range of the number of blocks to consider, and the optimal
  number of blocks is considered to be the maximum of the mean
  silhouette width of each hierarchical clustering. Silhouette width
  is calculated as per the \code{cluster} package. Hierarchical
  clustering is performed on the distance matrix calculated from the
  dissimilarity matrix (1 - abs(PostCor)) of the posterior correlation
  matrix. With sequential assignment, the number of parameters per
  block is approximately equal. With hierarchical clustering, the
  number of parameters per block may vary widely. Creating blocks from
  hierarchical clustering performs well in practice, though there are
  many alternative methods the user may consider outside of this
  function, such as using factor analysis, model-based clustering, or
  other methods.

  Aside from sequentially-assigned blocks, or blocks based on posterior
  correlation, it is also common to group parameters with similar uses,
  such as putting regression effects parameters into one block, and
  autocorrelation parameters into another block. Another popular way to
  group parameters into blocks is by time-period for some time-series
  models. These alternative blocking strategies are unsupported in the
  \code{Blocks} function, and best left to user discretion.

  Some MCMC algorithms that accept blocked parameters also require
  blocked variance-covariance matrices. The \code{Blocks} function
  does not return these matrices, because it may not be necessary,
  or when it is, the user may prefer identity matrices, scaled identity
  matrices, or matrices with explicitly-defined elements.

  If the user is looking for a place to begin with blockwise sampling,
  then the recommended, default approach (when blocked parameters by
  time-period are not desired in a time-series) is to begin with a
  trial run of the adaptive, unblocked HARM algorithm (since covariance
  matrices are not required) for the purposes of obtaining a posterior
  correlation matrix. Next, create blocks with the \code{Blocks}
  function based on the posterior correlation matrix obtained from the
  trial run. Finally, run the desired, blocked algorithm with the newly
  created blocks (and possibly user-specified covariance matrices),
  beginning where the trial run ended.

  If hierarchical clustering is used, then it is important to note that
  hierarchical clustering has no idea that the user intends to perform
  blockwise sampling in MCMC. If hierarchical clustering returns
  numerous small blocks, then the user may consider combining some or
  all of those blocks. For example, if several 1-parameter blocks are
  returned, then blockwise sampling will equal componentwise sampling
  for those blocks, which will iterate slower. Conversely, if
  hierarchical clustering returns one or more big blocks, each with
  enough parameters that multivariate sampling will have difficulty
  getting an accepted proposal, or an accepted proposal that moves
  more than a small amount, then the user may consider subdividing
  these big blocks into smaller, more manageable blocks, though with
  the understanding that more posterior correlation is unaccounted for.
}
\value{
  The \code{Blocks} function returns an object of class \code{blocks},
  which is a list. Each component of the list is a block of parameters,
  and parameters are indicated by their position in the initial values
  vector.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}} and
  \code{\link{PosteriorChecks}}.
}
\examples{
library(LaplacesDemon)

### Create the default number of sequentially assigned blocks:
Initial.Values <- rep(0,1000)
MyBlocks <- Blocks(Initial.Values)
MyBlocks

### Or, a pre-specified number of sequentially assigned blocks:
#Initial.Values <- rep(0,1000)
#MyBlocks <- Blocks(Initial.Values, N=20)

### If scaled diagonal covariance matrices are desired:
#VarCov <- list()
#for (i in 1:length(MyBlocks))
#  VarCov[[i]] <- diag(length(MyBlocks[[i]]))*2.38^2/length(MyBlocks[[i]])

### Or, determine the number of blocks in the range of 2 to 50 from
### hierarchical clustering on the posterior correlation matrix of an
### object, say called Fit, output from LaplacesDemon:
#MyBlocks <- Blocks(Initial.Values, N=c(2,50),
#  PostCor=cor(Fit$Posterior1))
#lapply(MyBlocks, length) #See the number of parameters per block

### Or, create a pre-specified number of blocks from hierarchical
### clustering on the posterior correlation matrix of an object,
### say called Fit, output from LaplacesDemon:
#MyBlocks <- Blocks(Initial.Values, N=20, PostCor=cor(Fit$Posterior1))

### Posterior correlation from a previous trial run could be obtained
### with either method below (though cor() will be fastest because
### additional checks are not calculated for the parameters):
#rho <- cor(Fit$Posterior1)
#rho <- PosteriorChecks(Fit)$Posterior.Correlation
}
\keyword{MCMC}
\keyword{Multicollinearity}
\keyword{Posterior Correlation}
\keyword{Utility}
