###########################################################################
# predict.vb                                                              #
#                                                                         #
# The purpose of the predict.vb function is to predict y[new] or y[rep],  #
# and later provide posterior predictive checks for objects of class vb.  #
###########################################################################

predict.vb <- function(object, Model, Data, CPUs=1, Type="PSOCK", ...)
     {
     ### Initial Checks
     if(missing(object)) stop("The object argument is required.")
     if(object$Converged == FALSE)
          stop("VariationalBayes did not converge.")
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(is.null(Data[["y"]]) & is.null(Data[["Y"]]))
          stop("Data must have y or Y.")
     if(!is.null(Data[["y"]])) y <- as.vector(Data[["y"]])
     if(!is.null(Data[["Y"]])) y <- as.vector(Data[["Y"]])
     CPUs <- abs(round(CPUs))
     ### p(y[rep] | y), Deviance, and Monitors
     Dev <- rep(NA, nrow(object$Posterior))
     monitor <- matrix(NA, length(Data[["mon.names"]]),
          nrow(object$Posterior))
     lengthcomp <- as.vector(Model(object$Posterior[1,], Data)[["yhat"]])
     if(!identical(length(lengthcomp), length(y)))
          stop("y and yhat differ in length.")
     yhat <- matrix(NA, length(y), nrow(object$Posterior))
     ### Non-Parallel Processing
     if(CPUs == 1) {
          for (i in 1:nrow(object$Posterior)) {
               mod <- Model(object$Posterior[i,], Data)
               Dev[i] <- as.vector(mod[["Dev"]])
               monitor[,i] <- as.vector(mod[["Monitor"]])
               yhat[,i] <- as.vector(mod[["yhat"]])}
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
          mod <- parLapply(cl, 1:nrow(object$Posterior),
               function(x) Model(object$Posterior[x,], Data))
          stopCluster(cl)
          Dev <- unlist(lapply(mod,
               function(x) x[["Dev"]]))[1:nrow(object$Posterior)]
          monitor <- matrix(unlist(lapply(mod,
               function(x) x[["Monitor"]])), length(Data[["mon.names"]]),
               nrow(object$Posterior))
          yhat <- matrix(unlist(lapply(mod,
               function(x) x[["yhat"]])), length(y),
               nrow(object$Posterior))
          rm(mod)}
     rownames(monitor) <- Data[["mon.names"]]
     ### Warnings
     if(any(is.na(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.na(yhat)), " missing values.")
     if(any(is.nan(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.nan(yhat)), " non-numeric (NaN) values.")
     if(any(is.infinite(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.infinite(yhat)), " infinite values.")
     if(any(!is.finite(Dev)))
          cat("\nWARNING: Deviance has non-finite values.")
     ### Create Output
     predicted <- list(y=y, yhat=yhat, Deviance=Dev,
          monitor=monitor)
     class(predicted) <- "vb.ppc"
     return(predicted)
     }

#End
