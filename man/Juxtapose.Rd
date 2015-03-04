\name{Juxtapose}
\alias{Juxtapose}
\title{Juxtapose MCMC Algorithm Inefficiency}
\description{
  This function gives a side-by-side comparison of (or juxtaposes) the
  inefficiency of MCMC algorithms in \code{\link{LaplacesDemon}} for
  applied use, and is a valuable tool for selecting what is likely to be
  the least inefficient algorithm for the user's current model, prior to
  updating the final, intended model.
}
\usage{
Juxtapose(x)
}
\arguments{
  \item{x}{This is a list of multiple components. Each component must
    be an object of class \code{demonoid}.}
}
\details{
  Laplace's Demon recommends using the \code{Juxtapose} function on the
  user's model (or most likely a simplified version of it) with a
  smaller, simulated data set to select the least inefficient MCMC
  algorithm before using real data and updating the model for numerous
  iterations. The least inefficient MCMC algorithm differs for different
  models and data sets. Using \code{Juxtapose} in this way does not
  guarantee that the selected algorithm will remain the best choice with
  real data, but it should be better than otherwise selecting an
  algorithm.

  The user must make a decision regarding their model and data. The more
  similar the model and data is to the final, intended model and data,
  the more appropriate will be the results of the \code{Juxtapose}
  function. However, if the full model and data are used, then the user
  may as well instead skip using \code{Juxtapose} and proceed directly
  to \code{\link{LaplacesDemon}}. Replacing the actual data set with a
  smaller, simulated set is fairly straightforward, but the
  decision-making will most likely focus on what is the best way to
  reduce the full model specification. A simple approach may be to
  merely reduce the number of predictors. However, complicated models
  may have several components that slow down estimation time, and extend
  the amount of time until global stationarity is estimated. Laplace's
  Demon offers no guidance here, and leaves it in the realm of user
  discretion.

  First, the user should simulate a smaller data set, and if best,
  reduce the model specification. Next, the user must select candidate
  algorithms. Then, the user must update each algorithm with
  \code{\link{LaplacesDemon}} for numerous iterations, with the goal of
  achieving stationarity for all parameters early in the
  iterations. Each update should begin with the same model specification
  function, vector of initial values, and data. Each output object of
  class \code{demonoid} should be renamed. An example follows.

  Suppose a user considers three candidate algorithms for their model:
  AMWG, NUTS, and twalk. The user updates each model, saving the model
  that used the AMWG algorithm as, say, \code{Fit1}, the NUTS model as
  \code{Fit2}, and the twalk model as \code{Fit3}.

  Next, the output model objects are put in a list and passed to the
  \code{Juxtapose} function. See the example below.

  The \code{Juxtapose} function uses an internal version of the
  \code{\link{IAT}}, which is a slightly modified version of that found
  in the \code{SamplerCompare} package. The \code{Juxtapose} function
  returns an object of class \code{juxtapose}. It is a matrix in which
  each row is a result and each column is an algorithm.

  The rows are:
  \itemize{
    \item \code{iter.min}: This is the iterations per minute.
    \item \code{t.iter.min}: This is the thinned iterations per minute.
    \item \code{prop.stat}: This is the proportion of iterations that
    were stationary.
    \item \code{IAT.025}: This is the 2.5\% quantile of the integrated
    autocorrelation time of the worst parameter, estimated only on
    samples when all parameters are estimated to be globally stationary.
    \item \code{IAT.500}: This is the median integrated autocorrelation
    time of the worst parameter, estimated only on samples when all
    parameters are estimated to be globally stationary.
    \item \code{IAT.975}: This is the 97.5\% quantile of the integrated
    autocorrelation time of the worst parameter, estimated only on
    samples when all parameters are estimated to be globally stationary.
    \item \code{ISM.025}: This is the 2.5\% quantile of the number of
    independent samples per minute.
    \item \code{ISM.500}: This is the median of the number of the
    independent samples per minute. The least inefficient MCMC algorithm
    has the highest \code{ISM.500}.
    \item \code{ISM.975}: This is the 97.5\% quantile of the number of
    the independent samples per minute.
  }

  As for calculating \eqn{ISM}, let \eqn{TIM} be the observed number of
  thinned iterations per minute, \eqn{PS} be the percent of iterations
  in which all parameters were estimated to be globally stationary, and
  \eqn{IAT_q}{IAT[q]} be a quantile from a simulated distribution of the
  integrated autocorrelation time among the parameters.

  \deqn{ISM = \frac{PS \times TIM}{IAT_q}}{ISM = (PS x TIM) / IAT[q]}
  
  There are various ways to measure the inefficiency of MCMC
  samplers. \code{\link{IAT}} is used perhaps most often. As with the
  \code{SamplerCompare} package, Laplace's Demon uses the worst
  parameter, in terms of \code{\link{IAT}}. Often, the number of
  evaluations or number of parameters is considered. The
  \code{Juxtapose} function, instead considers the final criterion of
  MCMC efficiency, in an applied context, to be \code{ISM}, or the
  number of Independent (thinned) Samples per Minute. The algorithm with
  the highest \code{ISM.500} is the best, or least inefficient,
  algorithm with respect to its worst \code{\link{IAT}}, the proportion
  of iterations required to seem to have global stationarity, and the
  number of (thinned) iterations per minute.

  A disadvantage of using time is that it will differ by computer, and
  is less likely to be reported in a journal. The advantage, though, is
  that it is more meaningful to a user. Increases in the number of
  evaluations, parameters, and time should all correlate well, but time
  may enlighten a user as to expected run-time given the model just
  studied, even though the real data set will most likely be larger than
  the simulated data used initially. NUTS is an example of a sampler in
  which the number of evaluations varies per iteration. For an
  alternative approach, see Thompson (2010).

  The \code{Juxtapose} function also adjusts \code{ISM} by
  \code{prop.stat}, the proportion of the iterations in which all chains
  were estimated to be stationary. This adjustment is weighted by
  burn-in iterations, penalizing an algorithm that took longer to
  achieve global stationarity. The goal, again, is to assist the user in
  selecting the least inefficient MCMC algorithm in an applied setting.

  The \code{Juxtapose} function has many other potential uses than those
  described above. One additional use of the \code{Juxtapose} function is
  to compare inefficiencies within a single algorithm in which
  algorithmic specifications varied with different model
  updates. Another use is to investigate parallel chains in an object of
  class \code{demonoid.hpc}, as returned from the
  \code{\link{LaplacesDemon.hpc}} function. Yet another use is to
  compare the effects of small changes to a model specification
  function, such as with priors, or due to an increase in the amount of
  simulated data.

  An object of class \code{juxtapose} may be plotted with the
  \code{\link{plot.juxtapose}} function, which displays \code{ISM} by
  default, or optionally \code{IAT}. For more information, see the
  \code{\link{plot.juxtapose}} function.

  Independent samples per minute, calculated as \code{\link{ESS}}
  divided by minutes of run-time, are also available by parameter in the
  \code{\link{PosteriorChecks}} function.
}
\value{
  This function returns an object of class \code{juxtapose}. It is a
  \eqn{9 \times J}{9 x J} matrix with nine results for \eqn{J} MCMC
  algorithms.
}
\references{
  Thompson, M. (2010). "Graphical Comparison of MCMC Performance". ArXiv
  e-prints, eprint 1011.4458.
}
\seealso{
  \code{\link{IAT}},
  \code{\link{is.juxtapose}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{plot.juxtapose}}, and
  \code{\link{PosteriorChecks}}.
}
\examples{
### Update three demonoid objects, each from different MCMC algorithms.
### Suppose Fit1 was updated with AFSS, Fit2 with AMWG, and
### Fit3 with NUTS. Then, compare the inefficiencies:
#Juxt <- Juxtapose(list(Fit1=Fit1, Fit2=Fit2, Fit3=Fit3)); Juxt
#plot(Juxt, Style="ISM")
}
\keyword{Diagnostic}
\keyword{Utilities}
