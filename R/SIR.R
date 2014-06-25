###########################################################################
# Sampling Importance Resampling (SIR)                                    #
#                                                                         #
# The purpose of the SIR function is to perform sampling importance       #
# re-sampling, usually to draw samples from the posterior as output from  #
# LaplaceApproxmation function. This function is similar to the sir       #
# function in the LearnBayes package.                                     #
###########################################################################

SIR <- function(Model, Data, mu, Sigma, n=1000, CPUs=1, Type="PSOCK") 
     {
     if(missing(Model)) stop("The Model function is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(missing(mu)) stop("The mu argument is required.")
     if(!is.vector(mu)) mu <- as.vector(mu)
     if(missing(Sigma)) stop("The Sigma argument is required.")
     if(!is.symmetric.matrix(Sigma)) Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma)) Sigma <- as.positive.definite(Sigma)
     if(length(mu) != nrow(Sigma)) stop("mu and Sigma are incompatible.")
     ### Sampling
     k <- length(mu)
     theta <- rmvn(n, mu, Sigma)
     theta[which(!is.finite(theta))] <- 0
     colnames(theta) <- Data[["parm.names"]]
     ### Importance
     lf <- matrix(0, n, 1)
     ### Non-Parallel Processing
     if(CPUs == 1) {
          for (i in 1:n) {
               mod <- Model(theta[i,], Data)
               lf[i] <- mod[["LP"]]
               theta[i,] <- mod[["parm"]]}
          }
     else { ### Parallel Processing
          detectedCores <- max(detectCores(),
               as.integer(Sys.getenv("NSLOTS")), na.rm=TRUE)
          cat("\n\nCPUs Detected:", detectedCores, "\n")
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n")
               CPUs <- detectedCores}
          cl <- makeCluster(CPUs, Type)
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
               ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          mod <- parLapply(cl, 1:nrow(theta),
               function(x) Model(theta[x,], Data))
          stopCluster(cl)
          lf <- unlist(lapply(mod,
               function(x) x[["LP"]]))[1:nrow(theta)]
          theta <- matrix(unlist(lapply(mod,
               function(x) x[["parm"]])), nrow(theta), ncol(theta))
          rm(mod)}
     lp <- dmvn(theta, mu, Sigma, log=TRUE)
     md <- max(lf - lp)
     lw <- lf - lp - md
     if(any(!is.finite(lw))) 
          lw[!is.finite(lw)] <- min(lw[is.finite(lw)])
     probs <- exp(lw - logadd(lw))
     ### Resampling
     options(warn=-1)
     indices <- try(sample.int(n, size=n, replace=TRUE, prob=probs),
          silent=TRUE)
     options(warn=0)
     if(inherits(indices, "try-error")) indices <- 1:n
     if(k > 1) theta <- theta[indices,]
     else theta <- theta[indices]
     return(theta)
     }

#End
