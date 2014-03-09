###########################################################################
# IAT                                                                     #
#                                                                         #
# The purpose of the IAT function is to estimate the integrated           #
# autocorrelation time of a chain, given its samples. Although the code   #
# is slightly different, it is essentially the same as the IAT function   #
# in the Rtwalk package, which is currently unavailable on CRAN.          #
###########################################################################

IAT <- function(x)
     {
     if(missing(x)) stop("The x argument is required.")
     if(!is.vector(x)) x <- as.vector(x)
     dt <- x
     n <- length(x)
     mu <- mean(dt)
     s2 <- var(dt)
     ### The maximum lag is half the sample size
     maxlag <- max(3, floor(n/2))
     #### The gammas are sums of two consecutive autocovariances
     Ga <- rep(0,2)
     Ga[1] <- s2
     lg <- 1
     Ga[1] <- Ga[1] + sum((dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu)) / n
     m <- 1
     lg <- 2*m
     Ga[2] <- sum((dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu)) / n
     lg <- 2*m + 1
     Ga[2] <- Ga[2] + sum((dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu)) / n
     IAT <- Ga[1]/s2 # Add the autocorrelations
     ### RULE: while Gamma stays positive and decreasing
     while ((Ga[2] > 0.0) & (Ga[2] < Ga[1])) {
          m <- m + 1
          if(2*m + 1 > maxlag) {
               cat("Not enough data, maxlag=", maxlag, "\n")
               break}
          Ga[1] <- Ga[2]
          lg <- 2*m
          Ga[2] <- sum((dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu)) / n
          lg <- 2*m + 1
          Ga[2] <- Ga[2] + sum((dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu)) / n
          IAT <- IAT + Ga[1] / s2
          }
     IAT <- -1 + 2*IAT #Calculates the IAT from the gammas
     return(IAT)
     }

#End
