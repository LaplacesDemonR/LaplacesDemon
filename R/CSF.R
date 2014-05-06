###########################################################################
# CSF                                                                     #
#                                                                         #
# The purpose of the CSF function is to provide a visual MCMC diagnostic  #
# based on the cumulative sample function (CSF).                          #
###########################################################################

CSF <- function(x, name, method="Quantiles", quantiles=c(.025,.5,.975),
     output=FALSE)
     {
     if(missing(x)) stop("The x argument is required.")
     if(missing(name)) name <- "x"
     if(is.constant(x)) stop("x must not be constant.")
     if(!is.vector(x)) x <- as.vector(x)
     if(method == "ESS") {
          y <- rep(0, length(x))
          for (i in 1:length(x)) {
               test <- try(ESS(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          plot(y, type="l", xlab="Cumulative Sample", ylab="ESS")
          if(output == TRUE) return(y)}
     if(method == "Geweke.Diagnostic") {
          y <- rep(0, length(x))
          for (i in 1:length(x)) {
               test <- try(Geweke.Diagnostic(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample",
               ylab="Geweke Diagnostic")
          abline(h=2, lty=2, col="red"); abline(h=-2, lty=2, col="red")
          if(output == TRUE) return(y)}
     if(method == "HPD") {
          Y <- matrix(0, length(x), 2)
          for (i in 1:length(x)) {
               test <- try(as.vector(p.interval(x[1:i], HPD=TRUE,
                    MM=FALSE)[1,]), silent=TRUE)
               if(!inherits(test, "try-error")) Y[i,] <- test}
          plot(x, type="l", col="gray", xlab="Sample Size",
               ylab="HPD (95%)")
          for (i in 1:2) lines(Y[,i], col="black")
          if(output == TRUE) return(Y)}
     if(method == "is.stationary") {
          y <- rep(FALSE, length(x))
          for (i in 1:length(x)) {
               test <- try(is.stationary(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample",
               ylab="Stationary Indicator")
          if(output == TRUE) return(y)}
     if(method == "Kurtosis") {
          kurtosis <- function(x) {  
               m4 <- mean((x-mean(x))^4) 
               kurt <- m4/(sd(x)^4)-3  
               return(kurt)}
          y <- rep(0, length(x))
          for (i in 1:length(x)) {
               test <- try(kurtosis(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample", ylab="Kurtosis")
          if(output == TRUE) return(y)}
     if(method == "MCSE") {
          y <- rep(1, length(x))
          for (i in 1:length(x)) {
               test <- try(MCSE(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          y[which(!is.finite(y))] <- 1
          y[which(y > 1)] <- 1
          plot(y, type="l", xlab="Cumulative Sample", ylab="MCSE")
          if(output == TRUE) return(y)}
     if(method == "MCSE.bm") {
          y <- rep(1, length(x))
          for (i in 1:length(x)) {
               test <- try(MCSE(x[1:i], method="batch.means"), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          y[which(!is.finite(y))] <- 1
          y[which(y > 1)] <- 1
          plot(y, type="l", xlab="Cumulative Sample", ylab="MCSE")
          if(output == TRUE) return(y)}
     if(method == "MCSE.sv") {
          y <- rep(1, length(x))
          for (i in 1:length(x)) {
               test <- try(MCSE(x[1:i], method="sample.variance"),
                    silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          y[which(!is.finite(y))] <- 1
          y[which(y > 1)] <- 1
          plot(y, type="l", xlab="Cumulative Sample", ylab="MCSE")
          if(output == TRUE) return(y)}
     if(method == "Mean") {
          y <- rep(0, length(x))
          for (i in 1:length(x)) {
               test <- try(mean(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample", ylab="Mean")
          if(output == TRUE) return(y)}
     if(method == "Mode") {
          y <- rep(0, length(x))
          for (i in 1:length(x)) {
               test <- try(Mode(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test[1]}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample", ylab="Mode")
          if(output == TRUE) return(y)}
     if(method == "N.Modes") {
          y <- rep(1, length(x))
          for (i in 1:length(x)) {
               test <- try(Modes(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- length(test$modes)}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample", ylab="Number of Modes")
          if(output == TRUE) return(y)}
     if(method == "Precision") {
          y <- rep(0, length(x))
          for (i in 1:length(x)) {
               test <- try(var(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- 1 / test}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample", ylab="Precision")
          if(output == TRUE) return(y)}
     if(method == "Quantiles") {
          Y <- matrix(0, length(x), length(quantiles))
          for (i in 1:length(x)) {
               test <- try(quantile(x[1:i], probs=quantiles), silent=TRUE)
               if(!inherits(test, "try-error")) Y[i,] <- test}
          plot(x, type="l", col="gray", xlab="Sample Size", ylab="Quantiles")
          for (i in 1:ncol(Y)) lines(Y[,i], col="black")
          if(output == TRUE) return(Y)}
     if(method == "Skewness") {
          skewness <-  function(x) {
               m3 <- mean((x-mean(x))^3)
               skew <- m3/(sd(x)^3)
               return(skew)}
          y <- rep(0, length(x))
          for (i in 1:length(x)) {
               test <- try(skewness(x[1:i]), silent=TRUE)
               if(!inherits(test, "try-error")) y[i] <- test}
          par(mfrow=c(2,1))
          plot(1:length(x), x, type="l", xlab="Iterations", ylab=name)
          panel.smooth(1:length(x), x, pch="")
          plot(y, type="l", xlab="Cumulative Sample", ylab="Skewness")
          if(output == TRUE) return(y)}
     }

#End
