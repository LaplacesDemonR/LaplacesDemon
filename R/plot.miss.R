###########################################################################
# plot.miss                                                               #
#                                                                         #
# The purpose of the plot.miss function is to plot an object of class     #
# miss.                                                                   #
###########################################################################

plot.miss <- function(x, PDF=FALSE, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(PDF == TRUE) {
          pdf("MISS.Plots.pdf")
          par(mfrow=c(3,3))
          }
     else par(mfrow=c(3,3), ask=TRUE)
     ### Plot Imputations
     for (i in 1:nrow(x$Imp)) {
          plot(1:ncol(x$Imp), x$Imp[i,], type="l", xlab="Iterations",
               ylab="Value", main=paste("Imp[", i, ",]", sep=""))
          panel.smooth(1:ncol(x$Imp), x$Imp[i,], pch="")
          plot(density(x$Imp[i,]), xlab="Value",
               main=paste("Imp[", i, ",]"))
          polygon(density(x$Imp[i,]), col="black", border="black")
          ### Only plot an ACF if there's > 1 unique values
          if(!is.constant(x$Imp[i,])) {
               z <- acf(x$Imp[i,], plot=FALSE)
               se <- 1/sqrt(length(x$Imp[i,]))
               plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
                    main=paste("Imp[", i, ",]"), xlab="Lag",
                    ylab="Correlation")
               abline(h=(2*se), col="red", lty=2)
               abline(h=(-2*se), col="red", lty=2)
               }
          else plot(0, 0, main=paste("Imp[", i, ",]"),
               "is a constant.")}
     if(PDF == TRUE) dev.off()
     }

#End
