###########################################################################
# predict.pmc                                                             #
#                                                                         #
# The purpose of the predict.pmc function is to predict y[new] or y[rep], #
# and later provide posterior predictive checks for objects of class pmc. #
###########################################################################

predict.pmc <- function(object, Model, Data, CPUs=1, Type="PSOCK", ...)
     {
     ### Initial Checks
     if(missing(object)) stop("The object argument is required.")
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(is.null(Data[["y"]]) & is.null(Data[["Y"]]))
          stop("Data must have y or Y.")
     if(!is.null(Data[["y"]])) y <- as.vector(Data[["y"]])
     if(!is.null(Data[["Y"]])) y <- as.vector(Data[["Y"]])
     CPUs <- abs(round(CPUs))
     ### p(y[rep] | y)
     post <- as.matrix(object$Posterior2)
     Dev <- rep(NA, nrow(post))
     yhat <- matrix(NA, length(y), nrow(post))
     lengthcomp <- as.vector(Model(post[1,], Data)[["yhat"]])
     if(!identical(length(lengthcomp), length(y)))
          stop("y and yhat differ in length.")
     ### Non-Parallel Processing
     if(CPUs == 1) {
          for (i in 1:nrow(post)) {
               mod <- Model(post[i,], Data)
               Dev[i] <- as.vector(mod[["Dev"]])
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
          mod <- parLapply(cl, 1:nrow(post),
               function(x) Model(post[x,], Data))
          stopCluster(cl)
          Dev <- unlist(lapply(mod,
               function(x) x[["Dev"]]))[1:nrow(post)]
          yhat <- matrix(unlist(lapply(mod,
               function(x) x[["yhat"]])), length(y), nrow(post))
          rm(mod)}
     ### Warnings
     if(any(is.na(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.na(yhat)), " missing values.\n")
     if(any(is.nan(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.nan(yhat)), " non-numeric (NaN) values.\n")
     if(any(is.infinite(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.infinite(yhat)), " infinite values.\n")
     if(any(!is.finite(Dev)))
          cat("\nWARNING: Deviance has non-finite values.")
     ### Create Output
     predicted <- list(y=y, yhat=yhat, Deviance=Dev)
     class(predicted) <- "pmc.ppc"
     return(predicted)
     }

#End
