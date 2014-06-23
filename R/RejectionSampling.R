###########################################################################
# RejectionSampling                                                       #
#                                                                         #
# The purpose of the RejectionSampling function is to perform rejection   #
# sampling.                                                               #
###########################################################################

RejectionSampling <- function(Model, Data, mu, S, df=Inf, logc, n=1000,
     CPUs=1, Type="PSOCK")
     {
     ### Initial Checks
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(missing(mu)) stop("The mu argument is required.")
     if(missing(S)) stop("The S argument is required.")
     if(missing(df)) stop("The df argument is required.")
     if(missing(logc)) stop("The logc argument is required.")
     df <- abs(df)
     ### Rejection Sampling
     k <- length(mu)
     if(df == Inf) theta <- rmvn(n, mu, S)
     else theta <- rmvt(n, mu, S, df)
     lf <- rep(0, nrow(theta))
     ### Non-Parallel Processing
     if(CPUs == 1) {
          for (i in 1:nrow(theta)) {
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
     if(df == Inf) lg <- dmvn(theta, mu, S, log=TRUE)
     else lg <- dmvt(theta, mu, S, df, log=TRUE)
     prob <- exp(lf - lg - logc)
     if(k == 1) theta <- theta[runif(n) < prob]
     else theta <- theta[runif(n) < prob,]
     theta <- class("rejection")
     return(theta)
     }

#End
