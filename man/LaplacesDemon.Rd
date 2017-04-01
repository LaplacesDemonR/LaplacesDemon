\name{LaplacesDemon}
\alias{LaplacesDemon}
\alias{LaplacesDemon.hpc}
\title{Laplace's Demon}
\description{
  The \code{LaplacesDemon} function is the main function of Laplace's
  Demon. Given data, a model specification, and initial values,
  \code{LaplacesDemon} maximizes the logarithm of the unnormalized joint
  posterior density with MCMC and provides samples of the marginal
  posterior distributions, deviance, and other monitored variables.

  The \code{LaplacesDemon.hpc} function extends \code{LaplacesDemon} to
  parallel chains for multicore or cluster high performance computing.
}
\usage{
LaplacesDemon(Model, Data, Initial.Values, Covar=NULL, Iterations=10000,
     Status=100, Thinning=10, Algorithm="MWG", Specs=list(B=NULL),
     Debug=list(DB.chol=FALSE, DB.eigen=FALSE, DB.MCSE=FALSE,
     DB.Model=TRUE), LogFile="", ...)
LaplacesDemon.hpc(Model, Data, Initial.Values, Covar=NULL,
     Iterations=10000, Status=100, Thinning=10, Algorithm="MWG",
     Specs=list(B=NULL), Debug=list(DB.chol=FALSE, DB.eigen=FALSE,
     DB.MCSE=FALSE, DB.Model=TRUE), LogFile="", Chains=2, CPUs=2,
     Type="PSOCK", Packages=NULL, Dyn.libs=NULL)
}
\arguments{
  \item{Model}{This required argument receives the model from a
    user-defined function that must be named Model. The user-defined
    function is where the model is specified. \code{LaplacesDemon}
    passes two arguments to the model function, \code{parms} and
    \code{Data}, and receives five arguments from the model function:
    \code{LP} (the logarithm of the unnormalized joint posterior),
    \code{Dev} (the deviance), \code{Monitor} (the monitored variables),
    \code{yhat} (the variables for posterior predictive checks), and
    \code{parm}, the vector of parameters, which may be constrained in
    the model function. More information on the Model specification
    function may be found in the "LaplacesDemon Tutorial" vignette, and
    the \code{\link{is.model}} function. Many examples of model
    specification functions may be found in the "Examples" vignette.}
  \item{Data}{This required argument accepts a list of data. The list of
    data must contain \code{mon.names} which contains monitored variable
    names, and must contain \code{parm.names} which contains parameter
    names. The \code{\link{as.parm.names}} function may be helpful for
    preparing the data, and the \code{\link{is.data}} function may be
    helpful for checking data.}
  \item{Initial.Values}{For \code{LaplacesDemon}, this argument requires
    a vector of initial values equal in length to the number of
    parameters. For \code{LaplacesDemon.hpc}, this argument also accepts
    a vector, in which case the same initial values will be applied to
    all parallel chains, or the argument accepts a matrix in which each
    row is a parallel chain and the number of columns is equal in length
    to the number of parameters. When a matrix is supplied for
    \code{LaplacesDemon.hpc}, each parallel chain begins with its own
    initial values that are preferably dispersed. For both
    \code{LaplacesDemon} and \code{LaplacesDemon.hpc}, each initial
    value will be the starting point for an adaptive chain or a
    non-adaptive Markov chain of a parameter. Parameters are assumed to
    be continuous, unless specified to be discrete (see \code{dparm}
    below), which is not accepted by all algorithms (see
    \code{\link{dcrmrf}} for an alternative). If all initial values are
    set to zero, then Laplace's Demon will attempt to optimize the
    initial values with the \code{\link{LaplaceApproximation}}
    function. After Laplace's Demon finishes updating, it may be desired
    to continue updating from where it left off. To continue, this
    argument should receive the last iteration of the previous update.
    For example, if the output object is called Fit, then
    \code{Initial.Values=as.initial.values(Fit)}. Initial values may be
    generated randomly with the \code{\link{GIV}} function.}
  \item{Covar}{This argument defaults to \code{NULL}, but may otherwise
    accept a \eqn{K \times K}{K x K} proposal covariance matrix
    (where \eqn{K} is the number of dimensions or parameters), a
    variance vector, or a list of covariance matrices (for blockwise
    sampling in some algorithms). When the model is updated for the
    first time and prior variance or covariance is unknown, then
    \code{Covar=NULL} should be used. Some algorithms require
    covariance, some only require variance, and some require neither.
    Laplace's Demon automatically converts the user input to the
    required form. Once Laplace's Demon has finished updating, it may
    be desired to continue updating where it left off, in which case
    the proposal covariance matrix from the last run can be input into
    the next run. The covariance matrix may also be input from the
    \code{\link{LaplaceApproximation}} function, if used.}
  \item{Iterations}{This required argument accepts integers larger than
    10, and determines the number of iterations that Laplace's Demon
    will update the parameters while searching for target
    distributions. The required amount of computer memory will increase
    with \code{Iterations}. If computer memory is exceeded, then all
    will be lost. The \code{\link{Combine}} function can be used later
    to combine multiple updates.}
  \item{Status}{This argument accepts an integer between 1 and the
    number of iterations, and indicates how often, in iterations, the
    user would like the status printed to the screen or log
    file. Usually, the following is reported: the number of iterations,
    the proposal type (for example, multivariate or componentwise, or
    mixture, or subset), and LP. For example, if a model is updated for
    1,000 iterations and \code{Status=200}, then a status message will
    be printed at the following iterations: 200, 400, 600, 800, and
    1,000.}
  \item{Thinning}{This argument accepts integers between 1 and the
    number of iterations, and indicates that every nth iteration will be
    retained, while the other iterations are discarded. If
    \code{Thinning=5}, then every 5th iteration will be
    retained. Thinning is performed to reduce autocorrelation and the
    number of marginal posterior samples.}
  \item{Algorithm}{This argument accepts the abbreviated name of the
    MCMC algorithm, which must appear in quotes. A list of MCMC
    algorithms appears below in the Details section, and the
    abbreviated name is in parenthesis.}
  \item{Specs}{This argument defaults to \code{NULL}, and accepts a list
    of specifications for the MCMC algorithm declared in the
    \code{Algorithm} argument. The specifications associated with each
    algorithm may be seen below in the examples, must appear in the
    order shown, and are described in the details section below.}
  \item{Debug}{This argument accepts a list of logical scalars that
    control whether or not errors or warnings are reported due to a
    \code{try} function or non-finite values. List components include
    \code{DB.chol} regarding \code{chol}, \code{DB.eigen} regarding
    \code{eigen}, \code{DB.MCSE} regarding \code{\link{MCSE}}, and
    \code{DB.Model} regarding the Model specification function. Errors
    and warnings should be investigated, but do not necessarily indicate
    a faulty Model specification function or a bug in the software. For
    example, a sampler may make a proposal that would result in a matrix
    that is not positive definite, when it should be. This kind of error
    or warning is acceptable, provided the sampler handles it correctly
    by rejecting the proposal, and provided the Model specification
    function is not causing the issue. Oftentimes, blockwise sampling
    with carefully chosen blocks will mostly or completely eliminate
    errors or warnings that occur otherwise in larger, multivariate
    proposals. Similarly, debugged componentwise algorithms tend to
    provide more information than multivariate algorithms, since
    usually the parameter and both its current and proposed values may
    be reported. If confident in the Model specification function, and
    errors or warnings are produced frequently that are acceptable,
    then consider setting \code{DB.Model=FALSE} for cleaner output and
    faster sampling. If the Model specification function is not faulty
    and there is a bug in \code{LaplacesDemon}, then please report it
    with a bug description and reproducible code on 
    \url{https://github.com/LaplacesDemonR/LaplacesDemon/issues}.}
  \item{LogFile}{This argument is used to specify a log file name in
    quotes in the working directory as a destination, rather than the
    console, for the output messages of \code{cat} and \code{stop}
    commands. It is helpful to assign a log file name when
    using multiple cores, such as with \code{LaplacesDemon.hpc}. Doing
    so allows the user to check the progress in the log. A number of log
    files are created, one for each chain, and one for the overall
    process.}
  \item{Chains}{This argument is required only for
    \code{LaplacesDemon.hpc}, and indicates the number of parallel
    chains.}
  \item{CPUs}{This argument is required for parallel independent or
    interactive chains in \code{LaplacesDemon} or
    \code{LaplacesDemon.hpc}, and indicates the number of central
    processing units (CPUs) of the computer or cluster. For example,
    when a user has a quad-core computer, \code{CPUs=4}.}
  \item{Type}{This argument defaults to \code{"PSOCK"} and uses the
    Simple Network of Workstations (SNOW) for parallelization.
    Alternatively, \code{Type="MPI"} may be specified to use Message
    Passing Interface (MPI) for parallelization.}
  \item{Packages}{This optional argument is for use with parallel
    independent or interacting chains, and defaults to \code{NULL}. This
    argument accepts a vector of package names to load into each
    parallel chain. If the \code{Model} specification depends on any
    packages, then these package names need to be in this vector.}
  \item{Dyn.libs}{This optional argument is for use with parallel
    independent or interacting chain, and defaults to \code{NULL}. This
    argument accepts a vector of the names of dynamic link libraries
    (shared objects) to load into each parallel chain. The libraries
    must be located in the working directory.}
  \item{...}{Additional arguments are unused.}
}
\details{
  \code{LaplacesDemon} offers numerous MCMC algorithms for numerical
  approximation in Bayesian inference. The algorithms are

  \itemize{
    \item Adaptive Directional Metropolis-within-Gibbs (ADMG)
    \item Adaptive Griddy-Gibbs (AGG)
    \item Adaptive Hamiltonian Monte Carlo (AHMC)
    \item Adaptive Metropolis (AM)
    \item Adaptive Metropolis-within-Gibbs (AMWG)
    \item Adaptive-Mixture Metropolis (AMM)
    \item Affine-Invariant Ensemble Sampler (AIES)
    \item Componentwise Hit-And-Run Metropolis (CHARM)
    \item Delayed Rejection Adaptive Metropolis (DRAM)
    \item Delayed Rejection Metropolis (DRM)
    \item Differential Evolution Markov Chain (DEMC)
    \item Elliptical Slice Sampler (ESS)
    \item Gibbs Sampler (Gibbs)
    \item Griddy-Gibbs (GG)
    \item Hamiltonian Monte Carlo (HMC)
    \item Hamiltonian Monte Carlo with Dual-Averaging (HMCDA)
    \item Hit-And-Run Metropolis (HARM)
    \item Independence Metropolis (IM)
    \item Interchain Adaptation (INCA)
    \item Metropolis-Adjusted Langevin Algorithm (MALA)
    \item Metropolis-Coupled Markov Chain Monte Carlo (MCMCMC)
    \item Metropolis-within-Gibbs (MWG)
    \item Multiple-Try Metropolis (MTM)
    \item No-U-Turn Sampler (NUTS)
    \item Oblique Hyperrectangle Slice Sampler (OHSS)
    \item Preconditioned Crank-Nicolson (pCN)
    \item Random Dive Metropolis-Hastings (RDMH)
    \item Random-Walk Metropolis (RWM)
    \item Reflective Slice Sampler (RSS)
    \item Refractive Sampler (Refractive)
    \item Reversible-Jump (RJ)
    \item Robust Adaptive Metropolis (RAM)
    \item Sequential Adaptive Metropolis-within-Gibbs (SAMWG)
    \item Sequential Metropolis-within-Gibbs (SMWG)
    \item Slice Sampler (Slice)
    \item Stochastic Gradient Langevin Dynamics (SGLD)
    \item Tempered Hamiltonian Monte Carlo (THMC)
    \item t-walk (twalk)
    \item Univariate Eigenvector Slice Sampler (UESS)
    \item Updating Sequential Adaptive Metropolis-within-Gibbs (USAMWG)
    \item Updating Sequential Metropolis-within-Gibbs (USMWG)
  }

  It is a goal for the documentation in the \pkg{LaplacesDemon} to be
  extensive. However, details of MCMC algorithms are best explored
  online at \url{https://web.archive.org/web/20150206014000/http://www.bayesian-inference.com/mcmc}, as well
  as in the "LaplacesDemon Tutorial" vignette, and the "Bayesian
  Inference" vignette. Algorithm specifications (\code{Specs}) are
  listed below:

  \itemize{
    \item \code{A} is used in AFSS, HMCDA, MALA, NUTS, OHSS, and UESS. In
    MALA, it is the maximum acceptable value of the Euclidean norm of
    the adaptive parameters mu and sigma, and the Frobenius norm of the
    covariance matrix. In AFSS, HMCDA, NUTS, OHSS, and UESS, it is the
    number of initial, adaptive iterations to be discarded as burn-in.
    \item \code{Adaptive} is the iteration in which adaptation begins,
    and is used in AM, AMM, DRAM, INCA, and Refractive. Most of these
    algorithms adapt according to an observed covariance matrix, and
    should sample before beginning to adapt.
    \item \code{alpha.star} is the target acceptance rate in MALA and
    RAM, and is optional in CHARM and HARM. The recommended value for
    multivariate proposals is \code{alpha.star=0.234}, for componentwise
    proposals is \code{alpha.star=0.44}, and for MALA is
    \code{alpha.star=0.574}.
    \item \code{at} affects the traverse move in twalk. \code{at=6} is
    recommended. It helps when some parameters are highly correlated,
    and the correlation structure may change through the
    state-space. The traverse move is associated with an acceptance rate
    that decreases as the number of parameters increases, and is the
    reason that \code{n1} is used to select a subset of parameters each
    iteration. If adjusted, it is recommended to stay in the interval
    [2,10].
    \item \code{aw} affects the walk move in twalk, and \code{aw=1.5} is
    recommended. If adjusted, it is recommended to stay in the
    interval [0.3,2].
    \item \code{beta} is a scale parameter for AIES, and defaults to 2,
    or an autoregressive parameter for pCN.
    \item \code{bin.n} is the scalar size parameter for a binomial prior
    distribution of model size for the RJ algorithm.
    \item \code{bin.p} is the scalar probability parameter for a
    binomial prior distribution of model size for the RJ algorithm.
    \item \code{B} is a list of blocked parameters. Each component of
    the list represents a block of parameters, and contains a vector in
    which each element is the position of the associated parameter in
    parm.names. This function is optional in the AFSS, AMM, AMWG, ESS,
    HARM, MWG, RAM, RWM, Slice, and UESS algorithms. For more
    information on blockwise sampling, see the \code{\link{Blocks}}
    function.
    \item \code{Begin} indicates the time-period in which to begin
    updating (filtering or predicting) in the USAMWG and USMWG
    algorithms.
    \item \code{Bounds} is used in the Slice algorithm. It is a vector
    of length two with the lower and upper boundary of the slice. For
    continuous parameters, it is often set to negative and positive
    infinity, while for discrete parameters it is set to the minimum
    and maximum discrete values to be sampled. When blocks are used,
    this must be supplied as a list with the same number of list
    components as the number of blocks.
    \item \code{delta} is used in HMCDA, MALA, and NUTS. In HMCDA and
    NUTS, it is the target acceptance rate, and the recommended value is
    0.65 in HMCDA and 0.6 in NUTS. In MALA, it is a constant in the
    bounded drift function, may be in the interval [1e-10,1000], and 1
    is the default.
    \item \code{Dist} is the proposal distribution in RAM, and may
    either be \code{Dist="t"} for t-distributed or \code{Dist="N"} for
    normally-distributed.
    \item \code{dparm} accepts a vector of integers that indicate
    discrete parameters. This argument is for use with the AGG or GG
    algorithm.
    \item \code{Dyn} is a \eqn{T \times K}{T x K} matrix of dynamic
    parameters, where \eqn{T} is the number of time-periods and \eqn{K}
    is the number of dynamic parameters. \code{Dyn} is used by SAMWG,
    SMWG, USAMWG, and USMWG. Non-dynamic parameters are updated first in
    each sampler iteration, then dynamic parameters are updated in a
    random order in each time-period, and sequentially by time-period.
    \item \code{epsilon} is used in AHMC, HMC, HMCDA, MALA, NUTS, SGLD,
    and THMC. It is the step-size in all algorithms except MALA. It is
    a vector equal in length to the number of parameters in AHMC, HMC,
    and THMC. It is a scalar in HMCDA and NUTS. It is either a scalar
    or a vector equal in length to the number of iterations in SGLD.
    When \code{epsilon=NULL} in HMCDA or NUTS (only), a reasonable
    initial value is found. In MALA, it is a vector of length two. The
    first element is the acceptable minimum of adaptive scale sigma, and
    the second element is added to the diagonal of the covariance matrix
    for regularization.
    \item \code{FC} is used in Gibbs and accepts a function that
    receives two arguments: the vector of all parameters and the list of
    data (similar to the Model specification function). FC must return
    the updated vector of all parameters. The user specifies FC to
    calculate the full conditional distribution of one or more
    parameters.
    \item \code{file} is the quoted name of a numeric matrix of data,
    without headers, for SGLD. The big data set must be a .csv
    file. This matrix has \code{Nr} rows and \code{Nc} columns. Each
    iteration, SGLD will randomly select a block of rows, where the
    number of rows is specified by the \code{size} argument.
    \item \code{Fit} is an object of class \code{demonoid} in the USAMWG
    and USMWG algorithms. Posterior samples before the time-period
    specified in the \code{Begin} argument are not updated, and are used
    instead from \code{Fit}.
    \item \code{gamma} controls the step size in DEMC or the decay of
    adaptation in MALA and RAM. In DEMC, it is positive and defaults to
    \eqn{2.38 / \sqrt{2J}}{2.38/sqrt(2J)} when \code{NULL}, where
    \eqn{J} is the length of initial values. For RAM, it is in the
    interval (0.5,1], and 0.66 is recommended. For MALA, it is in the
    interval (1,\code{Iterations}), and defaults to 1.
    \item \code{Grid} accepts either a vector or a list of vectors of
    evenly-spaced points on a grid for the AGG or GG algorithm. When the
    argument is a vector, the same grid is applied to all
    parameters. When the argument is a list, each component in the list
    has a grid that is applied to the corresponding parameter. The
    algorithm will evaluate each continuous parameter at the latest
    value plus each point in the grid, or each discrete parameter (see
    \code{dparm}) at each grid point (which should be each discrete
    value).
    \item \code{K} is a scalar number of proposals in MTM.
    \item \code{L} is a scalar number of leapfrog steps in AHMC, HMC, and
    THMC. When \code{L=1}, the algorithm reduces to Langevin Monte Carlo
    (LMC).
    \item \code{lambda} is used in HMCDA and MCMCMC. In HMCDA, it is a
    scalar trajectory length. In MCMCMC, it is either a scalar that
    controls temperature spacing, or a vector of temperature spacings.
    \item \code{Lmax} is a scalar maximum for \code{L} (see above) in
    HMCDA and NUTS.
    \item \code{m} is used in the AFSS, AHMC, HMC, Refractive, RSS, Slice,
    THMC, and UESS algorithms. In AHMC, HMC, and THMC, it is a
    \eqn{J \times J}{J x J} mass matrix for \eqn{J} initial values. In
    AFSS and UESS, it is a scalar, and is the maximum number of steps
    for creating the slice interval. In Refractive and RSS, it is a
    scalar, and is the number of steps to take per iteration. In Slice,
    it is either a scalar or a list with as many list components as
    blocks. It must be an integer in [1,Inf], and indicates the maximum
    number of steps for creating the slice interval.
    \item \code{mu} is a vector that is equal in length to the initial
    values. This vector will be used as the mean of the proposal
    distribution, and is usually the posterior mode of a
    previously-updated \code{\link{LaplaceApproximation}}.
    \item \code{MWG} is used in Gibbs to specify a vector of parameters
    that are to receive Metropolis-within-Gibbs updates. Each element is
    an integer that indicates the parameter.
    \item \code{Nc} is either the number of (un-parallelized) parallel
    chains in DEMC (and must be at least 3) or the number of columns of
    big data in SGLD.
    \item \code{Nr} is the number of rows of big data in SGLD.
    \item \code{n} is the number of previous iterations in ADMG, AFSS,
    AMM, AMWG, OHSS, RAM, and UESS.
    \item \code{n1} affects the size of the subset of each set of points
    to adjust, and is used in twalk. It relates to the number of
    parameters, and \code{n1=4} is recommended. If adjusted, it is
    recommended to stay in the interval [2,20].
    \item \code{parm.p} is a vector of probabilities for parameter
    selection in the RJ algorithm, and must be equal in length to
    the number of initial values.
    \item \code{r} is a scalar used in the Refractive algorithm to
    indicate the ratio between r1 and r2.
    \item \code{Periodicity} specifies how often in iterations the
    adaptive algorithm should adapt, and is used by AHMC, AM, AMM, AMWG,
    DRAM, INCA, SAMWG, and USAMWG. If \code{Periodicity=10}, then the
    algorithm adapts every 10th iteration. A higher \code{Periodicity}
    is associated with an algorithm that runs faster, because it does
    not have to calculate adaptation as often, though the algorithm
    adapts less often to the target distributions, so it is a
    trade-off. It is recommended to use the lowest value that runs
    fast enough to suit the user, or provide sufficient adaptation.
    \item \code{selectable} is a vector of indicators of whether or not
    a parameter is selectable for variable selection in the RJ
    algorithm. Non-selectable parameters are assigned a zero, and are
    always in the model. Selectable parameters are assigned a one. This
    vector must be equal in length to the number of initial values.
    \item \code{selected} is a vector of indicators of whether or not
    each parameter is selected when the RJ algorithm begins, and
    must be equal in length to the number of initial values.
    \item \code{SIV} stands for secondary initial values and is used by
    twalk. \code{SIV} must be the same length as \code{Initial.Values},
    and each element of these two vectors must be unique from each
    other, both before and after being passed to the \code{Model}
    function. \code{SIV} defaults to \code{NULL}, in which case values
    are generated with \code{\link{GIV}}.
    \item \code{size} is the number of rows of big data to be read into
    SGLD each iteration.
    \item \code{smax} is the maximum allowable tuning parameter sigma,
    the standard deviation of the conditional distribution, in the AGG
    algorithm.
    \item \code{Temperature} is used in the THMC algorithm to heat up
    the momentum in the first half of the leapfrog steps, and then cool
    down the momentum in the last half. \code{Temperature} must be
    positive. When greater than 1, THMC should explore more diffuse
    distributions, and may be helpful with multimodal distributions.
    \item \code{Type} is used in the Slice algorithm. It is either a
    scalar or a list with the same number of list components as blocks.
    This accepts \code{"Continuous"} for continuous parameters,
    \code{"Nominal"} for discrete parameters that are unordered, and
    \code{"Ordinal"} for discrete parameters that are ordered.
    \item \code{w} is used in AFSS, AMM, DEMC, Refractive, RSS, and
    Slice. It is a mixture weight for both the AMM and DEMC algorithms,
    and in these algorithms it is in the interval (0,1]. For AMM, it is
    recommended to use \code{w=0.05}, as per Roberts and Rosenthal
    (2009). The two mixture components in AMM are adaptive multivariate
    and static/symmetric univariate proposals. The mixture is determined
    at each iteration with mixture weight \code{w}. In the AMM
    algorithm, a higher value of \code{w} is associated with more
    static/symmetric univariate proposals, and a lower \code{w} is
    associated with more adaptive multivariate proposals. AMM will be
    unable to include the multivariate mixture component until it has
    accumulated some history, and models with more parameters will take
    longer to be able to use adaptive multivariate proposals. In DEMC,
    it indicates the probability that each iteration uses a snooker
    update, rather than a projection update, and the recommended default
    is \code{w=0.1}. In the Refractive algorithm, \code{w} is a scalar
    step size parameter. In AFSS, RSS, and the Slice algorithms, this is
    a step size interval for creating the slice interval. In AFSS and
    RSS, a scalar or vector equal in length the number of initial values
    is accepted. In Slice, a scalar or a list with a number of list
    components equal to the number of blocks is accepted.
    \item \code{Z} accepts a \eqn{T \times J}{T x J} matrix or \eqn{T
    \times J \times Nc}{T x J x Nc} array of thinned samples for \eqn{T}
    thinned iterations, \eqn{J} parameters, and \eqn{Nc} chains for
    DEMC. \code{Z} defaults to \code{NULL}. The matrix of thinned
    posterior samples from a previous run may be used, in which case the
    samples are copied across the chains.
  }
}
\value{
  \code{LaplacesDemon} returns an object of class \code{demonoid}, and
  \code{LaplacesDemon.hpc} returns an object of class
  \code{demonoid.hpc} that is a list of objects of class
  \code{demonoid}, where the number of components in the list
  is the number of parallel chains. Each object of class \code{demonoid} 
  is a list with the following components:
  \item{Acceptance.Rate}{This is the acceptance rate of the MCMC
    algorithm, indicating the percentage of iterations in which the
    proposals were accepted. For more information on acceptance rates,
    see the \code{\link{AcceptanceRate}} function.}
  \item{Algorithm}{This reports the specific algorithm used.}
  \item{Call}{This is the matched call of \code{LaplacesDemon}.}
  \item{Covar}{This stores the \eqn{K \times K}{K x K} proposal
    covariance matrix (where \eqn{K} is the dimension or number of
    parameters), variance vector, or list of covariance matrices.
    If variance or covariance is used for adaptation, then this
    covariance is returned. Otherwise, the variance of the samples of
    each parameter is returned. If the model is updated in the future,
    then this vector, matrix, or list can be used to start the next
    update where the last update left off. Only the diagonal of this
    matrix is reported in the associated \code{print} function.}
  \item{CovarDHis}{This \eqn{N \times K}{N x K} matrix stores the
    diagonal of the proposal covariance matrix of each adaptation in
    each of \eqn{N} rows for \eqn{K} dimensions, where the dimension is
    the number of parameters or length of the initial values vector. The
    proposal covariance matrix should change less over time. An
    exception is that the AHMC algorithm stores an algorithm
    specification here, which is not the diagonal of the proposal
    covariance matrix.}
  \item{Deviance}{This is a vector of the deviance of the model, with a
    length equal to the number of thinned samples that were retained.
    Deviance is useful for considering model fit, and is equal to the
    sum of the log-likelihood for all rows in the data set, which is
    then multiplied by negative two.}
  \item{DIC1}{This is a vector of three values: Dbar, pD, and DIC. Dbar
    is the mean deviance, pD is a measure of model complexity indicating
    the effective number of parameters, and DIC is the Deviance
    Information Criterion, which is a model fit statistic that is the
    sum of Dbar and pD. \code{DIC1} is calculated over all retained
    samples. Note that pD is calculated as \code{var(Deviance)/2} as in
    Gelman et al. (2004).}
  \item{DIC2}{This is identical to \code{DIC1} above, except that it is
    calculated over only the samples that were considered by the
    \code{BMK.Diagnostic} to be stationary for all parameters. If
    stationarity (or a lack of trend) is not estimated for all
    parameters, then \code{DIC2} is set to missing values.}
  \item{Initial.Values}{This is the vector of \code{Initial.Values},
    which may have been optimized with the
    \code{\link{IterativeQuadrature}} or
    \code{\link{LaplaceApproximation}} function.}
  \item{Iterations}{This reports the number of \code{Iterations} for
    updating.}
  \item{LML}{This is an approximation of the logarithm of the marginal
    likelihood of the data (see the \code{\link{LML}} function for more
    information). \code{LML} is estimated only with stationary samples,
    and only with a non-adaptive algorithm, including Adaptive
    Griddy-Gibbs (AGG), Affine-Invariant Ensemble Sampler (AIES),
    Componentwise Hit-And-Run (CHARM), Delayed Rejection Metropolis
    (DRM), Elliptical Slice Sampling (ESS), Gibbs Sampler (Gibbs),
    Griddy-Gibbs (GG), Hamiltonian Monte Carlo (HMC), Hit-And-Run
    Metropolis (HARM), Independence Metropolis (IM), Metropolis-Coupled
    Markov Chain Monte Carlo (MCMCMC), Metropolis-within-Gibbs (MWG),
    Multiple-Try Metropolis, No-U-Turn Sampler (NUTS), Random Dive
    Metropolis-Hastings (RDMH), Random-Walk Metropolis (RWM), Reflective
    Slice Sampler (RSS), Refractive Sampler (Refractive),
    Reversible-Jump (RJ), Sequential Metropolis-within-Gibbs (SMWG),
    Slice Sampler (Slice), Stochastic Gradient Langevin Dynamics (SGLD),
    Tempered Hamiltonian Monte Carlo (THMC), or t-walk (twalk).
    \code{LML} is estimated with nonparametric self-normalized
    importance sampling (NSIS), given LL and the marginal posterior
    samples of the parameters. \code{LML} is useful for comparing
    multiple models with the \code{\link{BayesFactor}} function.}
  \item{Minutes}{This indicates the number of minutes that
    \code{LaplacesDemon} was running, and includes the initial checks as
    well as time it took the \code{\link{LaplaceApproximation}}
    function, assessing stationarity, effective sample size (ESS), and
    creating summaries.}
  \item{Model}{This contains the model specification \code{Model}.}
  \item{Monitor}{This is a vector or matrix of one or more monitored
    variables, which are variables that were specified in the
    \code{Model} function to be observed as chains (or Markov chains,
    if \code{Adaptive=0}), but that were not deviance or parameters.}
  \item{Parameters}{This reports the number of parameters.}
  \item{Posterior1}{This is a matrix of marginal posterior distributions
    composed of thinned samples, with a number of rows equal to the
    number of thinned samples and a number of columns equal to the
    number of parameters. This matrix includes all thinned samples.}
  \item{Posterior2}{This is a matrix equal to \code{Posterior1}, except
    that rows are included only if stationarity (a lack of trend) is
    indicated by the \code{\link{BMK.Diagnostic}} for all parameters.
    If stationarity did not occur, then this matrix is missing.}
  \item{Rec.BurnIn.Thinned}{This is the recommended burn-in for the
    thinned samples, where the value indicates the first row that was
    stationary across all parameters, and previous rows are discarded
    as burn-in. Samples considered as burn-in are discarded because they
    do not represent the target distribution and have not adequately
    forgotten the initial value of the chain (or Markov chain, if
    \code{Adaptive=0}).}
  \item{Rec.BurnIn.UnThinned}{This is the recommended burn-in for all
    samples, in case thinning will not be necessary.}
  \item{Rec.Thinning}{This is the recommended value for the
    \code{Thinning} argument according to the autocorrelation in the
    thinned samples, and it is limited to the interval [1,1000].}
  \item{Specs}{This is an optional list of algorithm specifications.}
  \item{Status}{This is the value in the \code{Status} argument.}
  \item{Summary1}{This is a matrix that summarizes the marginal
    posterior distributions of the parameters, deviance, and monitored
    variables over all samples in \code{Posterior1}. The following
    summary statistics are included: mean, standard deviation, MCSE
    (Monte Carlo Standard Error), ESS is the effective sample size due
    to autocorrelation, and finally the 2.5\%, 50\%, and 97.5\%
    quantiles are reported. MCSE is essentially a standard deviation
    around the marginal posterior mean that is due to uncertainty
    associated with using MCMC. The acceptable size of the MCSE
    depends on the acceptable uncertainty associated around the
    marginal posterior mean. Laplace's Demon prefers to continue
    updating until each MCSE is less than 6.27\% of each marginal
    posterior standard deviation (see the \code{\link{MCSE}} and
    \code{\link{Consort}} functions). The default \code{IMPS} method
    is used. Next, the desired precision of ESS depends on the user's
    goal, and Laplace's Demon prefers to continue until each ESS is at
    least 100, which should be enough to describe 95\% boundaries of an
    approximately Gaussian distribution (see the \code{\link{ESS}} for
    more information).}
  \item{Summary2}{This matrix is identical to the matrix in
    \code{Summary1}, except that it is calculated only on the
    stationary samples found in \code{Posterior2}. If universal
    stationarity was not estimated for the parameters, then this matrix
    is set to missing values.}
  \item{Thinned.Samples}{This is the number of thinned samples that
    were retained.}
  \item{Thinning}{This is the value of the \code{Thinning} argument.}
}
\references{
  Atchade, Y.F. (2006). "An Adaptive Version for the Metropolis Adjusted
  Langevin Algorithm with a Truncated Drift". \emph{Methodology and
  Computing in Applied Probability}, 8, p. 235--254.

  Bai, Y. (2009). "An Adaptive Directional Metropolis-within-Gibbs
  Algorithm". Technical Report in Department of Statistics at the
  University of Toronto.

  Beskos, A., Roberts, G.O., Stuart, A.M., and Voss, J. (2008). "MCMC
  Methods for Diffusion Bridges". Stoch. Dyn., 8, p. 319--350.

  Boyles, L.B. and Welling, M. (2012). "Refractive Sampling".

  Craiu, R.V., Rosenthal, J., and Yang, C. (2009). "Learn From Thy
  Neighbor: Parallel-Chain and Regional Adaptive MCMC". \emph{Journal
  of the American Statistical Assocation}, 104(488), p. 1454--1466.
  
  Christen, J.A. and Fox, C. (2010). "A General Purpose Sampling
  Algorithm for Continuous Distributions (the t-walk)". \emph{Bayesian
  Analysis}, 5(2), p. 263--282.

  Dutta, S. (2012). "Multiplicative Random Walk Metropolis-Hastings
  on the Real Line". \emph{Sankhya B}, 74(2), p. 315--342.

  Duane, S., Kennedy, A.D., Pendleton, B.J., and Roweth, D. (1987).
  "Hybrid Monte Carlo". \emph{Physics Letters}, B, 195, p. 216--222.
  
  Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
  Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
  Hall, London.

  Geman, S. and Geman, D. (1984). "Stochastic Relaxation, Gibbs
  Distributions, and the Bayesian Restoration of Images". \emph{IEEE
  Transactions on Pattern Analysis and Machine Intelligence}, 6(6),
  p. 721--741.

  Geyer, C.J. (1991). "Markov Chain Monte Carlo Maximum Likelihood". In
  Keramidas, E.M. Computing Science and Statistics: Proceedings of the
  23rd Symposium of the Interface. Fairfax Station VA: Interface
  Foundation. p. 156--163.

  Goodman J, and Weare, J. (2010). "Ensemble Samplers with Affine
  Invariance". \emph{Communications in Applied Mathematics and
  Computational Science}, 5(1), p. 65--80.

  Green, P.J. (1995). "Reversible Jump Markov Chain Monte Carlo
  Computation and Bayesian Model Determination". \emph{Biometrika},
  82, p. 711--732.

  Haario, H., Laine, M., Mira, A., and Saksman, E. (2006). "DRAM:
  Efficient Adaptive MCMC". \emph{Statistical Computing}, 16,
  p. 339--354.
     
  Haario, H., Saksman, E., and Tamminen, J. (2001). "An Adaptive
  Metropolis Algorithm". \emph{Bernoulli}, 7, p. 223--242.

  Hoffman, M.D. and Gelman. A. (2012). "The No-U-Turn Sampler:
  Adaptively Setting Path Lengths in Hamiltonian Monte
  Carlo". \emph{Journal of Machine Learning Research}, p. 1--30.

  Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". \emph{Journal
  of the American Statistical Association}, 90(430), p. 773--795.

  Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
  Posterior Simulation with the Laplace-Metropolis Estimator".
  \emph{Journal of the American Statistical Association}, 92,
  p. 648--655.

  Liu, J., Liang, F., and Wong, W. (2000). "The Multiple-Try Method and
  Local Optimization in Metropolis Sampling". \emph{Journal of the
  American Statistical Association}, 95, p. 121--134.

  Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., and Teller,
  E. (1953). "Equation of State Calculations by Fast Computing
  Machines". \emph{Journal of Chemical Physics}, 21, p. 1087--1092.

  Mira, A. (2001). "On Metropolis-Hastings Algorithms with Delayed
  Rejection". \emph{Metron}, Vol. LIX, n. 3-4, p. 231--241.

  Murray, I., Adams, R.P., and MacKay, D.J. (2010). "Elliptical Slice
  Sampling". \emph{Journal of Machine Learning Research}, 9,
  p. 541--548.
  
  Neal, R.M. (2003). "Slice Sampling" (with discussion). \emph{Annals of
  Statistics}, 31(3), p. 705--767.

  Ritter, C. and Tanner, M. (1992), "Facilitating the Gibbs Sampler: the
  Gibbs Stopper and the Griddy-Gibbs Sampler", \emph{Journal of the
  American Statistical Association}, 87, p. 861--868.

  Roberts, G.O. and Rosenthal, J.S. (2009). "Examples of Adaptive
  MCMC". \emph{Computational Statistics and Data Analysis}, 18,
  p. 349--367.

  Roberts, G.O. and Tweedie, R.L. (1996). "Exponential Convergence of
  Langevin Distributions and Their Discrete Approximations".
  \emph{Bernoulli}, 2(4), p. 341--363.
  
  Rosenthal, J.S. (2007). "AMCMC: An R interface for adaptive MCMC".
  \emph{Computational Statistics and Data Analysis}, 51, p. 5467--5470.

  Smith, R.L. (1984). "Efficient Monte Carlo Procedures for Generating
  Points Uniformly Distributed Over Bounded Region". \emph{Operations
  Research}, 32, p. 1296--1308.

  Ter Braak, C.J.F. and Vrugt, J.A. (2008). "Differential Evolution
  Markov Chain with Snooker Updater and Fewer Chains", \emph{Statistics
  and Computing}, 18(4), p. 435--446.

  Tibbits, M., Groendyke, C., Haran, M., Liechty, J. (2014). "Automated
  Factor Slice Sampling". \emph{Journal of Computational and Graphical
  Statistics}, 23(2), p. 543--563.

  Thompson, M.D. (2011). "Slice Sampling with Multivariate Steps".
  \url{http://hdl.handle.net/1807/31955}

  Vihola, M. (2011). "Robust Adaptive Metropolis Algorithm with Coerced
  Acceptance Rate". \emph{Statistics and Computing}. Springer,
  Netherlands.

  Welling, M. and Teh, Y.W. (2011). "Bayesian Learning via Stochastic
  Gradient Langevin Dynamics". \emph{Proceedings of the 28th
  International Conference on Machine Learning (ICML)}, p. 681--688.

}
\author{Statisticat, LLC.,
     Silvere Vialet-Chabrand \email{silvere@vialet-chabrand.com}}
\seealso{
  \code{\link{AcceptanceRate}},
  \code{\link{as.initial.values}},
  \code{\link{as.parm.names}},
  \code{\link{BayesFactor}},
  \code{\link{Blocks}},
  \code{\link{BMK.Diagnostic}},
  \code{\link{Combine}},
  \code{\link{Consort}},
  \code{\link{dcrmrf}},
  \code{\link{ESS}},
  \code{\link{GIV}},
  \code{\link{is.data}},
  \code{\link{is.model}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon.RAM}},
  \code{\link{LML}}, and
  \code{\link{MCSE}}.
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
#library(compiler)
#Model <- cmpfun(Model) #Consider byte-compiling for more speed

set.seed(666)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

####################  Automated Factor Slice Sampler  #####################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))
Fit
print(Fit)
#Consort(Fit)
#plot(BMK.Diagnostic(Fit))
#PosteriorChecks(Fit)
#caterpillar.plot(Fit, Parms="beta")
#BurnIn <- Fit$Rec.BurnIn.Thinned
#plot(Fit, BurnIn, MyData, PDF=FALSE)
#Pred <- predict(Fit, Model, MyData, CPUs=1)
#summary(Pred, Discrep="Chi-Square")
#plot(Pred, Style="Covariates", Data=MyData)
#plot(Pred, Style="Density", Rows=1:9)
#plot(Pred, Style="ECDF")
#plot(Pred, Style="Fitted")
#plot(Pred, Style="Jarque-Bera")
#plot(Pred, Style="Predictive Quantiles")
#plot(Pred, Style="Residual Density")
#plot(Pred, Style="Residuals")
#Levene.Test(Pred)
#Importance(Fit, Model, MyData, Discrep="Chi-Square")

#############  Adaptive Directional Metropolis-within-Gibbs  ##############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="ADMG", Specs=list(n=0, Periodicity=50))

########################  Adaptive Griddy-Gibbs  ##########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AGG", Specs=list(Grid=GaussHermiteQuadRule(3)$nodes,
#     dparm=NULL, smax=Inf, CPUs=1, Packages=NULL, Dyn.libs=NULL))

##################  Adaptive Hamiltonian Monte Carlo  #####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AHMC", Specs=list(epsilon=0.02, L=2, m=NULL,
#     Periodicity=10))

##########################  Adaptive Metropolis  ##########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

###################  Adaptive Metropolis-within-Gibbs  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AMWG", Specs=list(B=NULL, n=0, Periodicity=50))

######################  Adaptive-Mixture Metropolis  ######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AMM", Specs=list(Adaptive=500, B=NULL, n=0,
#     Periodicity=10, w=0.05))

###################  Affine-Invariant Ensemble Sampler  ###################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AIES", Specs=list(Nc=2*length(Initial.Values), Z=NULL,
#     beta=2, CPUs=1, Packages=NULL, Dyn.libs=NULL))

#################  Componentwise Hit-And-Run Metropolis  ##################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="CHARM", Specs=NULL)

###########  Componentwise Hit-And-Run (Adaptive) Metropolis  #############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="CHARM", Specs=list(alpha.star=0.44))

#################  Delayed Rejection Adaptive Metropolis  #################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DRAM", Specs=list(Adaptive=500, Periodicity=10))

#####################  Delayed Rejection Metropolis  ######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DRM", Specs=NULL)

##################  Differential Evolution Markov Chain  ##################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DEMC", Specs=list(Nc=3, Z=NULL, gamma=NULL, w=0.1))

#######################  Elliptical Slice Sampler  ########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="ESS", Specs=list(B=NULL))

#############################  Gibbs Sampler  #############################
### NOTE: Unlike the other samplers, Gibbs requires specifying a
### function (FC) that draws from full conditionals.
#FC <- function(parm, Data)
#     {
#     ### Parameters
#     beta <- parm[Data$pos.beta]
#     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
#     sigma2 <- sigma*sigma
#     ### Hyperparameters
#     betamu <- rep(0,length(beta))
#     betaprec <- diag(length(beta))/1000
#     ### Update beta
#     XX <- crossprod(Data$X)
#     Xy <- crossprod(Data$X, Data$y)
#     IR <- backsolve(chol(XX/sigma2 + betaprec), diag(length(beta)))
#     btilde <- crossprod(t(IR)) %*% (Xy/sigma2 + betaprec %*% betamu)
#     beta <- btilde + IR %*% rnorm(length(beta))
#     return(c(beta,sigma))
#     }
##library(compiler)
##FC <- cmpfun(FC) #Consider byte-compiling for more speed
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="Gibbs", Specs=list(FC=FC, MWG=pos.sigma))


#############################  Griddy-Gibbs  ##############################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="GG", Specs=list(Grid=seq(from=-0.1, to=0.1, len=5),
#     dparm=NULL, CPUs=1, Packages=NULL, Dyn.libs=NULL))

#######################  Hamiltonian Monte Carlo  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HMC", Specs=list(epsilon=0.001, L=2, m=NULL))

#############  Hamiltonian Monte Carlo with Dual-Averaging  ###############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=1, Thinning=1,
#     Algorithm="HMCDA", Specs=list(A=500, delta=0.65, epsilon=NULL,
#     Lmax=1000, lambda=0.1))

#######################  Hit-And-Run Metropolis  ##########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HARM", Specs=NULL)

##################  Hit-And-Run (Adaptive) Metropolis  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HARM", Specs=list(alpha.star=0.234, B=NULL))

########################  Independence Metropolis  ########################
### Note: the mu and Covar arguments are populated from a previous Laplace
### Approximation.
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=Fit$Covar, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="IM",
#     Specs=list(mu=Fit$Summary1[1:length(Initial.Values),1]))

#########################  Interchain Adaptation  #########################
#Initial.Values <- rbind(Initial.Values, GIV(Model, MyData, PGF=TRUE))
#Fit <- LaplacesDemon.hpc(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="INCA", Specs=list(Adaptive=500, Periodicity=10),
#     LogFile="MyLog", Chains=2, CPUs=2, Type="PSOCK", Packages=NULL,
#     Dyn.libs=NULL)

################  Metropolis-Adjusted Langevin Algorithm  #################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MALA", Specs=list(A=1e7, alpha.star=0.574, gamma=1,
#          delta=1, epsilon=c(1e-6,1e-7)))

#############  Metropolis-Coupled Markov Chain Monte Carlo  ###############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MCMCMC", Specs=list(lambda=1, CPUs=2, Packages=NULL,
#     Dyn.libs=NULL))

#######################  Metropolis-within-Gibbs  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MWG", Specs=list(B=NULL))

########################  Multiple-Try Metropolis  ########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MTM", Specs=list(K=4, CPUs=1, Packages=NULL, Dyn.libs=NULL))

##########################  No-U-Turn Sampler  ############################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=1, Thinning=1,
#     Algorithm="NUTS", Specs=list(A=500, delta=0.6, epsilon=NULL,
#     Lmax=Inf))

#################  Oblique Hyperrectangle Slice Sampler  ##################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="OHSS", Specs=list(A=Inf, n=0))

#####################  Preconditioned Crank-Nicolson  #####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="pCN", Specs=list(beta=0.1))

######################  Robust Adaptive Metropolis  #######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RAM", Specs=list(alpha.star=0.234, B=NULL, Dist="N",
#     gamma=0.66, n=0))

###################  Random Dive Metropolis-Hastings  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RDMH", Specs=NULL)

##########################  Refractive Sampler  ###########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="Refractive", Specs=list(Adaptive=1, m=2, w=0.1, r=1.3))

###########################  Reversible-Jump  #############################
#bin.n <- J-1
#bin.p <- 0.2
#parm.p <- c(1, rep(1/(J-1),(J-1)), 1)
#selectable <- c(0, rep(1,J-1), 0)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RJ", Specs=list(bin.n=bin.n, bin.p=bin.p,
#          parm.p=parm.p, selectable=selectable,
#          selected=c(0,rep(1,J-1),0)))

########################  Random-Walk Metropolis  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RWM", Specs=NULL)

########################  Reflective Slice Sampler  #######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RSS", Specs=list(m=5, w=1e-5))

##############  Sequential Adaptive Metropolis-within-Gibbs  ##############
#NOTE: The SAMWG algorithm is only for state-space models (SSMs)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="SAMWG", Specs=list(Dyn=Dyn, Periodicity=50))

##################  Sequential Metropolis-within-Gibbs  ###################
#NOTE: The SMWG algorithm is only for state-space models (SSMs)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="SMWG", Specs=list(Dyn=Dyn))

#############################  Slice Sampler  #############################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=1, Thinning=1,
#     Algorithm="Slice", Specs=list(B=NULL, Bounds=c(-Inf,Inf), m=100,
#     Type="Continuous", w=1))

#################  Stochastic Gradient Langevin Dynamics  #################
#NOTE: The Data and Model functions must be coded differently for SGLD.
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=10, Thinning=10,
#     Algorithm="SGLD", Specs=list(epsilon=1e-4, file="X.csv", Nr=1e4,
#     Nc=6, size=10))

###################  Tempered Hamiltonian Monte Carlo  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="THMC", Specs=list(epsilon=0.001, L=2, m=NULL,
#     Temperature=2))

###############################  t-walk  #################################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="twalk", Specs=list(SIV=NULL, n1=4, at=6, aw=1.5))

#################  Univariate Eigenvector Slice Sampler  #################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="UESS", Specs=list(A=Inf, B=NULL, m=100, n=0))

##########  Updating Sequential Adaptive Metropolis-within-Gibbs  #########
#NOTE: The USAMWG algorithm is only for state-space model updating
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
#     Covar=NULL, Iterations=100000, Status=100, Thinning=100,
#     Algorithm="USAMWG", Specs=list(Dyn=Dyn, Periodicity=50, Fit=Fit,
#     Begin=T.m))

##############  Updating Sequential Metropolis-within-Gibbs  ##############
#NOTE: The USMWG algorithm is only for state-space model updating
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
#     Covar=NULL, Iterations=100000, Status=100, Thinning=100,
#     Algorithm="USMWG", Specs=list(Dyn=Dyn, Fit=Fit, Begin=T.m))

#End
}
\keyword{Adaptive}
\keyword{Adaptive Directional Metropolis-within-Gibbs}
\keyword{Adaptive Griddy-Gibbs}
\keyword{Adaptive Hamiltonian Monte Carlo}
\keyword{Adaptive Metropolis}
\keyword{Adaptive-Mixture Metropolis}
\keyword{Adaptive Metropolis-within-Gibbs}
\keyword{Affine-Invariant Ensemble Sampler}
\keyword{Bayesian Inference}
\keyword{Delayed Rejection Adaptive Metropolis}
\keyword{Delayed Rejection Metropolis}
\keyword{Differential Evolution Markov Chain}
\keyword{Elliptical Slice Sampler}
\keyword{Gibbs Sampler}
\keyword{Griddy-Gibbs}
\keyword{Hamiltonian Monte Carlo}
\keyword{High Performance Computing}
\keyword{Hit-And-Run}
\keyword{Independence Metropolis}
\keyword{Interchain Adaptation}
\keyword{MCMC}
\keyword{Metropolis-Adjusted Langevin Algorithm}
\keyword{Metropolis-Coupled Markov Chain Monte Carlo}
\keyword{Metropolis-within-Gibbs}
\keyword{Multiple-Try Metropolis}
\keyword{Multiple Chains}
\keyword{No-U-Turn Sampler}
\keyword{Oblique Hyperrectangle Slice Sampler}
\keyword{Optimization}
\keyword{Parallel Chains}
\keyword{Preconditioned Crank-Nicolson}
\keyword{Reflective Slice Sampler}
\keyword{Refractive Sampler}
\keyword{Reversible-Jump}
\keyword{Robust Adaptive Metropolis}
\keyword{Random-Walk Metropolis}
\keyword{Sequential Adaptive Metropolis-within-Gibbs}
\keyword{Sequential Metropolis-within-Gibbs}
\keyword{Slice Sampler}
\keyword{Stochastic Gradient Langevin Dynamics}
\keyword{Tempered Hamiltonian Monte Carlo}
\keyword{t-walk}
\keyword{Univariate Eigenvector Slice Sampler}
\keyword{Updating Sequential Adaptive Metropolis-within-Gibbs}
\keyword{Updating Sequential Metropolis-within-Gibbs}
