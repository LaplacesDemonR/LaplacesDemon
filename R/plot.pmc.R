###########################################################################
# plot.pmc                                                                #
#                                                                         #
# The purpose of the plot.pmc function is to plot an object of class pmc. #
###########################################################################

plot.pmc <- function(x, BurnIn=0, Data=NULL, PDF=FALSE, Parms=NULL, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(class(x) != "pmc")
          stop("x must be of class pmc.")
     if(is.null(Data)) stop("The Data argument is NULL.")
     if(BurnIn >= nrow(x$Posterior2)) BurnIn <- 0
     Stat.at <- BurnIn + 1
     ### Selecting Parms
     if(is.null(Parms)) {
          Posterior <- x$Posterior2
          keepcols <- 1:ncol(Posterior)}
     else {
          Parms <- sub("\\[","\\\\[",Parms)
          Parms <- sub("\\]","\\\\]",Parms)
          Parms <- sub("\\.","\\\\.",Parms)
          if(length(grep(Parms[1], colnames(x$Posterior2))) == 0)
               stop("Parameter in Parms does not exist.")
          keepcols <- grep(Parms[1], colnames(x$Posterior2))
          if(length(Parms) > 1) {
               for (i in 2:length(Parms)) {
                    if(length(grep(Parms[i],
                         colnames(x$Posterior2))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keepcols <- c(keepcols,
                         grep(Parms[i], colnames(x$Posterior2)))}}
          Posterior <- as.matrix(x$Posterior2[,keepcols])
          colnames(Posterior) <- colnames(x$Posterior2)[keepcols]}
     if(PDF == TRUE) {
          pdf("PMC.Plots.pdf")
          par(mfrow=c(2,2))
          }
     else {par(mfrow=c(2,2), ask=TRUE)}
     ### Plot Parameters
     for (j in 1:ncol(Posterior)) {
          ### Plot Parameter Trace Plots
          k <- keepcols[j]
          LL <- Me <- UL <- matrix(0, x$M, x$Iterations)
          for (m in 1:x$M) {for (i in Stat.at:x$Iterations) {
               LL[m,i] <- as.vector(quantile(x$Posterior1[,k,i,m],
                    probs=0.025))
               Me[m,i] <- as.vector(quantile(x$Posterior1[,k,i,m],
                    probs=0.500))
               UL[m,i] <- as.vector(quantile(x$Posterior1[,k,i,m],
                    probs=0.975))}}
          plot(Stat.at:x$Iterations, Me[1,Stat.at:x$Iterations],
               ylim=c(min(LL[,Stat.at:x$Iterations]),
                    max(UL[,Stat.at:x$Iterations])), pch=20,
               xlab="Iterations", ylab="Value",
               main=colnames(Posterior)[j])
          for (i in Stat.at:x$Iterations) lines(c(i,i), c(LL[1,i], UL[1,i]))
          if(x$M > 1) {
               for (m in 2:x$M) {
                    points(Stat.at:x$Iterations+(m-1)*0.1,
                         Me[m,Stat.at:x$Iterations], col=m, pch=20)
                    for (i in Stat.at:x$Iterations) {
                         lines(c(i+(m-1)*0.1,i+(m-1)*0.1),
                              c(LL[m,i], UL[m,i]), col=m)}}}
          ### Plot Parameter Densities
          plot(density(Posterior[Stat.at:x$Thinned.Samples,j]),
               xlab="Value", main=colnames(Posterior)[j])
          polygon(density(Posterior[Stat.at:x$Thinned.Samples,j]),
               col="black", border="black")
          abline(v=0, col="red", lty=2)
          }
     ### Plot Deviance Density
     plot(density(x$Deviance[Stat.at:length(x$Deviance)]),
          xlab="Value", main="Deviance")
     polygon(density(x$Deviance[Stat.at:length(x$Deviance)]), col="black",
               border="black")
     abline(v=0, col="red", lty=2)
     ### Plot Monitored Variable Densities
     if(is.vector(x$Monitor)) {J <- 1; nn <- length(x$Monitor)}
     else if(is.matrix(x$Monitor)) {
          J <- ncol(x$Monitor); nn <- nrow(x$Monitor)}
     for (j in 1:J) {
          plot(density(x$Monitor[Stat.at:nn,j]),
               xlab="Value", main=Data[["mon.names"]][j])
          polygon(density(x$Monitor[Stat.at:nn,j]), col="black",
               border="black")
          abline(v=0, col="red", lty=2)}
     ### Plot Convergence Diagnostics
     plot(x$Perplexity, ylim=c(0,1), type="l", xlab="Iterations", ylab="",
          sub="Perplexity=black; ESSN=red", main="Convergence")
     lines(x$ESSN, col="red")
     W <- Thin(x$W, By=x$Thinning)
     boxplot(W, outline=FALSE, col="red", xlab="Iterations",
          ylab="Importance Weights")
     if(x$M > 1) {
          plot(x$alpha[1,], col=1, ylim=c(0,1), type="l",
               main="Mixture Probabilities", xlab="Iterations",
               ylab="alpha")
          for (m in 2:x$M) {
               lines(x$alpha[m,], col=m)}}
     if(PDF == TRUE) dev.off()
     }

#End
