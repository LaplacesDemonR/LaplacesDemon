###########################################################################
# MCSE                                                                    #
#                                                                         #
# The purpose of the MCSE function is to estimate the Monte Carlo         #
# Standard Error of a vector of posterior samples. Multiple methods are   #
# provided. The purpose of the MCSS function is to calculate the required #
# sample size `n' to achieve acceptable error `a' given a vector `x' of   #
# Monte Carlo samples and the sample variance (rather than the asymptotic #
# variance).                                                              #
###########################################################################

MCSE <- function(x, method="IMPS", batch.size="sqrt", warn=FALSE)
     {
     if(missing(x)) stop("The x argument is required.")
     if(method == "sample.variance") {
          ess <- try(ESS(x), silent=TRUE)
          if(inherits(ess, "try-error")) ess <- length(x)
          se <- sd(x) / sqrt(ess)
          return(se)}
     else if(method == "batch.means") {
          N <- length(x)
          if(N < 1000) if(warn) warning("Samples must be >= 1000.")
          if(N < 10) return(NA)
          if(batch.size == "sqrt") {
               b <- floor(sqrt(N)) # batch size
               a <- floor(N/b) # number of batches
               }
          else if(batch.size == "cuberoot") {
                    b <- floor(N^(1/3)) # batch size
                    a <- floor(N/b) # number of batches
                    }
          else { #Batch size is provided numerically
               stopifnot(is.numeric(batch.size))  
               b <- floor(batch.size) # batch size
               if(b > 1) a <- floor(N/b) # number of batches
               else stop("batch.size is invalid.")
               }
          Ys <- sapply(1:a, function(k) return(mean(x[((k-1)*b+1):(k*b)])))
          muhat <- mean(Ys)
          sigmahatsq <- b*sum((Ys - muhat)^2) / (a-1)
          bmse <- sqrt(sigmahatsq / N)
          return(list(est=muhat, se=bmse))
          }
     else if(method == "IMPS") {
          chainAC <- acf(x, type="covariance" ,plot=FALSE)$acf
          AClen <- length(chainAC)
          gammaAC <- chainAC[1:(AClen-1)] + chainAC[2:AClen]
          m <- 1
          currgamma <- gammaAC[1]
          k <- 1
          while ((k < length(gammaAC)) && (gammaAC[k+1] > 0) &&
               (gammaAC[k] >= gammaAC[k+1]))
               k <- k + 1
          if({k == length(gammaAC)} & {warn == TRUE})
               warning("May need to compute more autocovariances for IMPS.")
          options(warn=-1)
          sigmasq <- -chainAC[1] + 2*sum(gammaAC[1:k])
          se <- sqrt(sigmasq / length(x))
          options(warn=0)
          return(se)
          }
     else stop("The method is unknown.")
     }
MCSS <- function(x, a)
     {
     if(missing(x)) stop("The x argument is required.")
     if(missing(a)) stop("The a argument is required.")
     ess <- ESS(x)
     ratio <- length(x) / ess
     fx <- function(sdx, n, a) ((sdx / sqrt(n)) - a)^2
     optimized <- optimize(f=fx, interval=c(0, .Machine$integer.max),
          maximum=FALSE, a=a, sdx=sd(x))
     return(round(optimized$minimum * ratio))
     }

#End
