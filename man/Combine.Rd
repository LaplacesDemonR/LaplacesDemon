\name{Combine}
\alias{Combine}
\title{Combine Demonoid Objects}
\description{
  This function combines objects of class \code{demonoid}.
}
\usage{
Combine(x, Data, Thinning=1)
}
\arguments{
  \item{x}{This is a list of objects of class \code{demonoid}, and this
    list may be an object of class \code{demonoid.hpc}.}
  \item{Data}{This is the data, and must be identical to the data used
    to create the \code{demonoid} objects with
    \code{\link{LaplacesDemon}}.}
  \item{Thinning}{This is the amount of thinning to apply to the
    posterior samples after appending them together. \code{Thinning}
    defaults to 1, in which case all samples are retained. For example,
    in the case of, say, \code{Thinning=10}, then only every 10th sample
    would be retained. When combining parallel chains, \code{Thinning}
    is often left to its default. When combining consecutive updates,
    \code{Thinning} is usually applied, with the value equal to the
    number of objects of class \code{demonoid}. For more information on
    thinning, see the \code{Thin} function.}
}
\details{
  The purpose of the \code{Combine} function is to enable a user to
  combine objects of class \code{demonoid} for one of three
  reasons. First, parallel chains from \code{\link{LaplacesDemon.hpc}}
  may be combined after convergence is assessed with
  \code{\link{Gelman.Diagnostic}}. Second, consecutive updates of single
  chains from \code{\link{LaplacesDemon}} or parallel chains from
  \code{\link{LaplacesDemon.hpc}} may be combined when the computer has
  insufficient random-access memory (RAM) for the user to update once 
  with enough iterations. Third, consecutive single-chain or
  parallel-chain updates may be combined when it seems that the
  logarithm of the joint posterior distribution, \code{LP}, seems to be
  oscillating up and down, which is described in more detail below.

  The most common use regards the combination of parallel chains output
  from \code{\link{LaplacesDemon.hpc}}. Typically, a user with parallel
  chains examines them graphically with the
  \code{\link{caterpillar.plot}} and \code{plot} (actually,
  \code{\link{plot.demonoid}}) functions, and assesses convergence
  with the \code{\link{Gelman.Diagnostic}} function. Thereafter, the
  parallel chain output in the object of class \code{demonoid.hpc}
  should be combined into a single object of class \code{demonoid},
  before doing posterior predictive checks and making inferences. In
  this case, the \code{Thinning} argument usually is recommended to
  remain at its default.

  It is also common with a high-dimensional model (a model with a large
  number of parameters) to need more posterior samples than allowed by
  the random-access memory (RAM) of the computer. In this case, it is
  best to use the \code{\link{LaplacesDemon.RAM}} function to estimate
  the amount of RAM that a given model will require with a given number
  of iterations, and then update \code{\link{LaplacesDemon}} almost as
  much as RAM allows, and save the output object of class
  \code{demonoid}. Then, the user is advised to continue onward with a
  consecutive update (after using \code{\link{as.initial.values}} and
  anything else appropriate to prepare for the consecutive
  update). Suppose a user desires to update a gigantic model with
  thousands of parameters, and with the aid of
  \code{\link{LaplacesDemon.RAM}}, estimates that they can safely update
  only 100,000 iterations, and that 150,000 iterations would exceed RAM
  and crash the computer. The patient user can update several
  consecutive models, each with retaining only 1,000 thinned posterior
  samples, and combine them later with the \code{Combine} function, by
  placing multiple objects into a list, as described below. In this way,
  it is possible for a user to update models that otherwise far exceed
  computer RAM.

  Less commonly, multiple updates of single-chain objects should be
  combined into a single object of class \code{demonoid}. This is most
  useful in complicated models that are run for large numbers of
  iterations, where it may be suspected that stationarity has been
  achieved, but that thinning is insufficient, and the samples may be
  combined and thinned. If followed, then these suggestions may continue
  seemingly to infinity, and the unnormalized logarithm of the joint
  posterior density, \code{LP}, may seem to oscillate, sometimes
  improving and getting higher, and getting lower during other updates.
  For this purpose, the prior covariance matrix of the last model is
  retained (rather than combining them). This may be an unpleasant
  surprise for combining parallel updates, so be aware of it.

  In these cases, which usually involve complicated models with high
  autocorrelation in the chains, the user may opt to use parallel
  processing with the \code{\link{LaplacesDemon.hpc}} function, or may
  use the \code{\link{LaplacesDemon}} function as follows. The user
  should save (meaning, not overwrite) each object of class
  \code{demonoid}, place multiple objects into a list, and use the
  \code{Combine} function to combine these objects.

  For example, suppose a user names the object Fit, as in the
  \code{\link{LaplacesDemon}} example. Now, rather than overwriting
  object Fit, object Fit is renamed, after updating a million
  iterations, to Fit1. As suggested by \code{\link{Consort}}, another
  million iterations are used, but now to create object Fit2. Further
  suppose this user specified \code{Thinning=1000} in
  \code{\link{LaplacesDemon}}, meaning that the million iterations are
  thinned by 1,000, so only 1,000 iterations are retained in each
  object, Fit1 and Fit2. In this case, \code{Combine} combines the
  information in Fit1 and Fit2, and returns an object the user names
  Fit3. Fit3 has only 1,000 iterations, which is the result of appending
  the iterations in Fit1 and Fit2, and thinning by 2. If 2,000,000
  iterations were updated from the beginning, and were thinned by 2,000,
  then the same information exists now in Fit3. The \code{\link{Consort}}
  function can now be applied to Fit3, to see if stationarity is found.
  If not, then more objects of class \code{demonoid} can be collected and
  combined.
}
\value{
  This function returns an object of class \code{demonoid}. For more
  information on an object of class \code{demonoid}, see the
  \code{LaplacesDemon} function.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{\code{\link{caterpillar.plot}},
  \code{\link{Gelman.Diagnostic}},
  \code{\link{LaplacesDemon}}, 
  \code{\link{LaplacesDemon.hpc}}, and
  \code{\link{Thin}}.}
\keyword{MCMC}
\keyword{Multiple Chains}
\keyword{Parallel Chains}
\keyword{Utility}