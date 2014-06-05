###########################################################################
# plot.iterquad                                                           #
#                                                                         #
# The purpose of the plot.iterquad function is to plot an object of class #
# iterquad.                                                               #
###########################################################################

plot.iterquad <- function(x, Data=NULL, PDF=FALSE, Parms=NULL, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(class(x) != "iterquad")
          stop("x must be of class iterquad.")
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
          pdf("IterativeQuadrature.Plots.pdf")
          par(mfrow=c(2,2))
          }
     else {par(mfrow=c(2,2), ask=TRUE)}
     ### Plot Parameter
     for (j in 1:ncol(History))
          {
          plot(1:nrow(History), History[,j],
               type="l", xlab="Iterations", ylab="Value",
               main=colnames(History)[j])
          LPw <- x$LPw[,j]
          Mw <- x$M[,j]
          if(sum(Mw) > 0) Mw <- Mw / sum(Mw)
          Z <- x$Z[,j]
          LPw <- as.vector(by(LPw, Z, max))
          Mw <- as.vector(by(Mw, Z, max))
          Z <- unique(Z)
          o <- order(Z)
          LPw <- LPw[o]
          Mw <- Mw[o]
          Z <- Z[o]
          covdens <- dnorm(Z, History[nrow(History),j],
               sqrt(diag(x$Covar))[j])
          if(sum(covdens) > 0) covdens <- covdens / sum(covdens)
          if({x$Converged == TRUE} & !any(is.na(Posterior))) {
               dens <- density(Posterior[,j])
               if(sum(dens$y) > 0) dens$y <- dens$y / sum(dens$y)
               x.lim <- range(c(Z,dens$x))
               y.lim <- c(0, max(Mw, LPw, covdens, dens$y))
               }
          else {
               x.lim <- range(Z)
               y.lim <- c(0, max(Mw, LPw, covdens))
               }
          plot(Z, Mw, type="h", xlim=x.lim, ylim=y.lim, col="white",
               main=colnames(History)[j], xlab="Value",
               ylab="Normalized Density")
          polygon(c(Z, rev(Z)), c(rep(0,length(LPw)), rev(covdens)),
               col=rgb(0,255,0,50,maxColorValue=255), border=NA)
          if({x$Converged == TRUE} & !any(is.na(Posterior)))
               polygon(dens,
                    col=rgb(0,0,255,50,maxColorValue=255), border=NA)
          polygon(c(Z, rev(Z)),
               c(rep(0,length(Mw)),rev(Mw)),
               col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          polygon(c(Z, rev(Z)),
               c(rep(0,length(LPw)),rev(LPw)),
               col=rgb(255,0,0,50,maxColorValue=255), border=NA)
          abline(v=0, col="red", lty=2)
          lines(Z, Mw, type="h")
          lines(Z, LPw, type="h", col="red")
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
