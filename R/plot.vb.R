###########################################################################
# plot.vb                                                                 #
#                                                                         #
# The purpose of the plot.vb function is to plot an object of class vb.   #
###########################################################################

plot.vb <- function(x, Data=NULL, PDF=FALSE, Parms=NULL, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(class(x) != "vb") stop("x must be of class vb.")
     if(is.null(Data)) stop("The Data argument is NULL.")
     if(any(is.na(x$History))) stop("There is no history to plot.")
     ### Selecting Parms
     if(is.null(Parms)) {
          History1 <- x$History[,,1]
          History2 <- x$History[,,2]
          Posterior <- x$Posterior}
     else {
          Parms <- sub("\\[","\\\\[",Parms)
          Parms <- sub("\\]","\\\\]",Parms)
          Parms <- sub("\\.","\\\\.",Parms)
          if(length(grep(Parms[1], dimnames(x$History)[[2]])) == 0)
               stop("Parameter in Parms does not exist.")
          keepcols <- grep(Parms[1], dimnames(x$History[[2]]))
          if(length(Parms) > 1) {
               for (i in 2:length(Parms)) {
                    if(length(grep(Parms[i],
                         dimnames(x$History)[[2]])) == 0)
                         stop("Parameter in Parms does not exist.")
                    keepcols <- c(keepcols, grep(Parms[i],
                         dimnames(x$History)[[2]]))}}
          History1 <- as.matrix(x$History[,keepcols,1])
          History2 <- as.matrix(x$History[,keepcols,2])
          colnames(History1) <- dimnames(x$History)[[2]][keepcols]
          colnames(History2) <- dimnames(x$History)[[2]][keepcols]
          if(all(!is.na(x$Posterior))) {
               Posterior <- as.matrix(x$Posterior[,keepcols])
               colnames(Posterior) <- colnames(History1)}
          else Posterior <- x$Posterior
          }
     if(PDF == TRUE)
          {
          pdf("VariationalBayes.Plots.pdf")
          par(mfrow=c(3,3))
          }
     else {par(mfrow=c(3,3), ask=TRUE)}
     ### Plot Parameter
     for (j in 1:ncol(History1))
          {
          plot(1:nrow(History1), History1[,j],
               type="l", xlab="Iterations", ylab="Value (Mean)",
               main=colnames(History1)[j])
          plot(1:nrow(History2), History2[,j],
               type="l", xlab="Iterations", ylab="Value (Variance)",
               main=colnames(History2)[j])
          if({x$Converged == TRUE} & !any(is.na(Posterior))) {
               plot(density(Posterior[,j]),
                    xlab="Value", main=colnames(Posterior)[j])
               polygon(density(Posterior[,j]),
                    col="black", border="black")
               abline(v=0, col="red", lty=2)
               }
          else {
               draws <- rnorm(1000, x$Summary1[j,1], x$Summary1[j,2])
               plot(density(draws), xlab="Value",
                    main=colnames(History1)[j])
               polygon(density(draws), col="black", border="black")
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
