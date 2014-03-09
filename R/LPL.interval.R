###########################################################################
# LPL.interval                                                            #
#                                                                         #
# The purpose of the LPL.interval function is to estimate the lowest      #
# posterior loss (LPL) interval.                                          #
###########################################################################

LPL.interval <- function(Prior, Posterior, prob=0.95, plot=FALSE,
     PDF=FALSE)
     {
     ### Initial Checks
     if(missing(Prior)) stop("Prior is required.")
     if(missing(Posterior)) stop("Posterior is required.")
     if(!is.vector(Prior)) Prior <- as.vector(Prior)
     if(!is.vector(Posterior)) Posterior <- as.vector(Posterior)
     if(length(Prior) != length(Posterior))
          stop("Length mismatch between Prior and Posterior.")
     if(any(!is.finite(Prior) | !is.finite(Posterior)))
          stop("Non-finite values found in Prior or Posterior.")
     name <- names(Posterior)
     if(is.null(name)) name <- "Value"
     ### Expected Posterior Loss
     ord <- order(Posterior)
     Prior <- Prior[ord]
     Posterior <- Posterior[ord]
     loss <- KLD(Prior, Posterior)[[3]]
     ### Plot Expected Posterior Loss
     if(plot == TRUE) {
          if(PDF == TRUE) pdf("LPL.Plot.pdf")
          par(mfrow=c(2,1))
          plot(Posterior, loss, type="l", main="Posterior Loss",
               xlab=name, ylab="E(Posterior Loss)")
          polygon(c(min(Posterior), Posterior, max(Posterior)),
               c(min(loss), loss, min(loss)),
               col="gray", border="gray")}
     ### Find LPL Interval
     n <- length(loss)
     gap <- max(1, min(n - 1, round(n * prob)))
     loss.sum <- init <- 1:(n - gap)
     for (i in 1:length(init)) {
          loss.sum[i] <- sum(loss[init[i]:(init[i]+gap)])}
     min.init <- init[which.min(loss.sum)]
     ans <- cbind(Posterior[min.init], Posterior[min.init+gap])
     colnames(ans) <- c("Lower","Upper")
     attr(ans, "LPL.Interval") <- prob
     ### Shade LPL Area
     if(plot == TRUE) {
          polygon(c(min(Posterior[min.init]),
               Posterior[min.init:(min.init+gap)],
               max(Posterior[min.init+gap])),
               c(min(loss[min.init:(min.init+gap)]),
               loss[min.init:(min.init+gap)],
               min(loss[min.init:(min.init+gap)])),
               col="black", border="black")
          abline(v=0, col="red", lty=2)
          kde <- kde.low <- kde.high <- density(Posterior)
          kde.low$x <- kde$x[kde$x < ans[1,1]]
          kde.low$y <- kde$y[which(kde$x < ans[1,1])]
          kde.high$x <- kde$x[kde$x > ans[1,2]]
          kde.high$y <- kde$y[which(kde$x > ans[1,2])]
          plot(kde, xlab=name, ylab="Density",
               main="LPL Probability Interval")
          polygon(kde, col="black", border="black")
          polygon(c(min(kde.low$x), kde.low$x, max(kde.low$x)),
               c(min(kde.low$y), kde.low$y, min(kde.low$y)),
               col="gray", border="gray")
          polygon(c(min(kde.high$x), kde.high$x, max(kde.high$x)),
               c(min(kde.high$y), kde.high$y, min(kde.high$y)),
               col="gray", border="gray")
          abline(v=0, col="red", lty=2)
          if(PDF == TRUE) dev.off()}
     return(ans)
     }

#End
