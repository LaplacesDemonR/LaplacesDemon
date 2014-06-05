###########################################################################
# plot.laplace                                                            #
#                                                                         #
# The purpose of the plot.laplace function is to plot an object of class  #
# laplace.                                                                #
###########################################################################

plot.laplace <- function(x, Data=NULL, PDF=FALSE, Parms=NULL, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(class(x) != "laplace")
          stop("x must be of class laplace.")
     if(is.null(Data)) stop("The Data argument is NULL.")
     if(any(is.na(x$History))) stop("There is no history to plot.")
     ### Selecting Parms
     if(is.null(Parms)) {
          History <- x$History
          Posterior <- x$Posterior}
     else {
          Parms <- sub("\\[","\\\\[",Parms)
          Parms <- sub("\\]","\\\\]",Parms)
          Parms <- sub("\\.","\\\\.",Parms)
          if(length(grep(Parms[1], colnames(x$History))) == 0)
               stop("Parameter in Parms does not exist.")
          keepcols <- grep(Parms[1], colnames(x$History))
          if(length(Parms) > 1) {
               for (i in 2:length(Parms)) {
                    if(length(grep(Parms[i],
                         colnames(x$History))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keepcols <- c(keepcols, grep(Parms[i],
                         colnames(x$History)))}}
          History <- as.matrix(x$History[,keepcols])
          colnames(History) <- colnames(x$History)[keepcols]
          if(all(!is.na(x$Posterior))) {
               Posterior <- as.matrix(x$Posterior[,keepcols])
               colnames(Posterior) <- colnames(History)}
          else Posterior <- x$Posterior
          }
     if(PDF == TRUE)
          {
          pdf("LaplaceApproximation.Plots.pdf")
          par(mfrow=c(2,2))
          }
     else {par(mfrow=c(2,2), ask=TRUE)}
     ### Plot Parameter
     for (j in 1:ncol(History))
          {
          plot(1:nrow(History), History[,j],
               type="l", xlab="Iterations", ylab="Value",
               main=colnames(History)[j])
          if({x$Converged == TRUE} & !any(is.na(Posterior))) {
               plot(density(Posterior[,j]),
                    xlab="Value", main=colnames(Posterior)[j])
               polygon(density(Posterior[,j]),
                    col="black", border="black")
               abline(v=0, col="red", lty=2)}
          }
     ### Plot Deviance History
     plot(1:length(x$Deviance), x$Deviance, type="l", xlab="Iterations",
          ylab="Value", main="Deviance")
     ### Plot Monitor
     if({x$Converged == TRUE} & !any(is.na(x$Monitor))) {
          for (j in 1:ncol(x$Monitor)) {
               plot(density(x$Monitor[,j]),
                    xlab="Value", main=Data[["mon.names"]][j])
               polygon(density(x$Monitor[,j]),
                    col="black", border="black")
               abline(v=0, col="red", lty=2)}
          }
     if(PDF == TRUE) dev.off()
     }

#End
