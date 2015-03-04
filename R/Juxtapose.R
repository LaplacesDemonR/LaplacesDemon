###########################################################################
# Juxtapose                                                               #
#                                                                         #
# The purpose of the Juxtapose function is to compare the inefficiency of #
# multiple updates of LaplacesDemon, each with a different algorithm, but #
# the same model, initial values, and data. The internal ar.act and       #
# ar.act1 functions are slightly modified versions of those in the        #
# SamplerCompare package (modified to use the rmvn function), and the     #
# method of assessing inefficiency differs.                               #
###########################################################################

Juxtapose <- function(x)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(!is.list(x)) stop("The x argument must be a list.")
     ### IAT with probability intervals
     ar.act1  <- function(y)
          {
          stopifnot(NCOL(y) == 1)
          if(length(unique(y)) < 5)
               return(list(act=NA, act.025=NA, act.975=NA, se=NA,
                    order=NA))
          order.max <- NULL
          repeat {
               A <- ar.yw(y, demean=FALSE, order.max=order.max)
               if(A$order == 0)
                    A <- ar.yw(y, demean=FALSE, order.max=1, aic=FALSE)
               pi <- A$ar
               pi.var <- A$asy.var.coef
               if(kappa(pi.var) < 1 / sqrt(.Machine$double.eps) ||
                    isTRUE(order.max == 1))
                    break
               order.max <- floor(A$order * 0.7)}
          acf <- matrix(ARMAacf(ar=pi)[2:(A$order+1)])
          act <- (1-sum(pi*acf))/(1-sum(pi))^2
          simulation.length <- min(max(40, length(y)), 5000)
          pi.var2 <- as.symmetric.matrix(pi.var)
          AX <- rmvn(simulation.length, pi, pi.var2)
          act.sim <- numeric(simulation.length)
          for (i in 1:simulation.length) {
               pi.sim <- AX[i,]
               acf.sim <- ARMAacf(ar=pi.sim)[2:(A$order+1)]
               if(any(abs(polyroot(c(-1,pi.sim)))<1)) act.sim[i] <- Inf
               else act.sim[i] <- (1-sum(pi.sim*acf.sim)) /
                    (1-sum(pi.sim))^2}
          act.sim[is.na(act.sim)] <- Inf
          act.025 <- as.numeric(quantile(act.sim, 0.025))
          act.975 <- as.numeric(quantile(act.sim, 0.975))
          se <- (act.975-act.025)/(2*1.96)
          return(list(act=act, se=se, act.025=act.025, act.975=act.975,
               order=A$order))
          }
     ar.act <- function(Y, true.mean=NULL)
          {
          Y <- as.matrix(Y)
          stopifnot(is.null(true.mean) || ncol(Y)==length(true.mean))
          if(is.null(true.mean)) mu <- colMeans(Y)
          else mu <- true.mean
          acts <- sapply(1:ncol(Y), function(i) ar.act1(Y[,i]-mu[i]))
          max.i <- which.max(unlist(acts['act',]))
          if(length(max.i)!=1) {
               return(list(act=NA, act.025=NA, act.975=NA, se=NA,
                    order=NA))}
          else return(acts[,max.i])
          }
     ### Back to Juxtapose
     lenx <- length(x)
     ### Set up Output
     out <- matrix(NA, 9, lenx)
     algs <- rep(NA, lenx)
     for (i in 1:lenx) {
          if(!identical(class(x[[i]]), "demonoid"))
               stop("A component of x was found not be of class demonoid.")
          ### Use the abbreviated name of the algorithm
          algs[i] <- switch(x[[i]][["Algorithm"]],
               "Adaptive Directional Metropolis-within-Gibbs"="ADMG",
               "Adaptive Griddy-Gibbs"="AGG",
               "Adaptive Hamiltonian Monte Carlo"="AHMC",
               "Adaptive Metropolis"="AM",
               "Adaptive Metropolis-within-Gibbs"="AMWG",
               "Adaptive-Mixture Metropolis"="AMM",
               "Affine-Invariant Ensemble Sampler"="AIES",
               "Automated Factor Slice Sampler"="AFSS",
               "Componentwise Hit-And-Run Metropolis"="CHARM",
               "Delayed Rejection Adaptive Metropolis"="DRAM",
               "Delayed Rejection Metropolis"="DRM",
               "Differential Evolution Markov Chain"="DEMC",
               "Elliptical Slice Sampler"="ESS",
               "Experimental"="Exper",
               "Gibbs Sampler"="Gibbs",
               "Griddy-Gibbs"="GG",
               "Hamiltonian Monte Carlo"="HMC",
               "Hamiltonian Monte Carlo with Dual-Averaging"="HMCDA",
               "Hit-And-Run Metropolis"="HARM",
               "Independence Metropolis"="IM",
               "Interchain Adaptation"="INCA",
               "Metropolis-Adjusted Langevin Algorithm"="MALA",
               "Metropolis-Coupled Markov Chain Monte Carlo"="MCMCMC",
               "Metropolis-within-Gibbs"="MWG",
               "Multiple-Try Metropolis"="MTM",
               "No-U-Turn Sampler"="NUTS",
               "Oblique Hyperrectangle Slice Sampler"="OHSS",
               "Preconditioned Crank-Nicolson"="pCN",
               "Random Dive Metropolis-Hastings"="RDMH",
               "Random-Walk Metropolis"="RWM",
               "Reflective Slice Sampler"="RSS",
               "Refractive Sampler"="Refractive",
               "Reversible-Jump"="RJ",
               "Robust Adaptive Metropolis"="RAM",
               "Sequential Adaptive Metropolis-within-Gibbs"="SAMWG",
               "Sequential Metropolis-within-Gibbs"="SMWG",
               "Slice Sampler"="Slice",
               "Stochastic Gradient Langevin Dynamics"="SGLD",
               "Tempered Hamiltonian Monte Carlo"="THMC",
               "t-walk"="t-walk",
               "Univariate Eigenvector Slice Sampler"="UESS",
               "Updating Sequential Adaptive Metropolis-within-Gibbs"="USAMWG",
               "Updating Sequential Metropolis-within-Gibbs"="USMWG")
          }
     colnames(out) <- algs
     rownames(out) <- c("iter.min","t.iter.min","prop.stat","IAT.025",
          "IAT.500","IAT.975","ISM.025","ISM.500","ISM.975")
     class(out) <- "juxtapose"
     ### Juxtapose Algorithms
     for (i in 1:lenx) {
          iter.min <- x[[i]][["Iterations"]] / x[[i]][["Minutes"]]
          t.iter.min <- x[[i]][["Iterations"]] / x[[i]][["Thinning"]] /
               x[[i]][["Minutes"]]
          if(x[[i]][["Rec.BurnIn.Thinned"]] >= x[[i]][["Thinned.Samples"]])
               prop.stat <- 0
          else prop.stat <- 1 - (x[[i]][["Rec.BurnIn.Thinned"]] /
               x[[i]][["Thinned.Samples"]])
          if(all(is.na(x[[i]][["Summary2"]]))) {
               iat.500 <- iat.025 <- iat.975 <- Inf}
          else {
               iat.temp <- ar.act(x[[i]][["Posterior2"]])
               iat.500 <- iat.temp$act
               iat.025 <- iat.temp$act.025
               iat.975 <- iat.temp$act.975}
          ism.025 <- prop.stat * t.iter.min / iat.975
          ism.500 <- prop.stat * t.iter.min / iat.500
          ism.975 <- prop.stat * t.iter.min / iat.025
          out[1,i] <- round(iter.min,2)
          out[2,i] <- round(t.iter.min,2)
          out[3,i] <- round(prop.stat,2)
          out[4,i] <- round(iat.025,2)
          out[5,i] <- round(iat.500,2)
          out[6,i] <- round(iat.975,2)
          out[7,i] <- round(ism.025,2)
          out[8,i] <- round(ism.500,2)
          out[9,i] <- round(ism.975,2)
          }
     return(out)
     }

#End
