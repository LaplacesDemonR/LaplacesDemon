###########################################################################
# burnin                                                                  #
#                                                                         #
# The purpose of the burnin function is to estimate the duration of       #
# burn-in in iterations for one or more MCMC chains.                      #
###########################################################################

burnin <- function(x, method="BMK")
     {
     if(missing(x)) stop("The x argument is required.")
     if(is.vector(x)) x <- matrix(x, length(x), 1)
     n <- nrow(x)
     burn <- rep(0,ncol(x))
     if(method == "BMK") {
          if(n %% 10 == 0) x2 <- x
          if(n %% 10 != 0) x2 <- x[1:(10*trunc(n/10)),]
          HD <- BMK.Diagnostic(x2, 10)
          Ind <- 1 * (HD > 0.5)
          burn <- n
          batch.list <- seq(from=1, to=nrow(x2), by=floor(nrow(x2)/10))
          for (i in 1:9) {
               if(sum(Ind[,i:9]) == 0) {
                    burn <- batch.list[i] - 1
                    break
                    }
               }
          }
     else {
          for (i in 1:ncol(x)) {
               iter <- 1
               stationary <- 0
               jump <- round(n/10)
               while(stationary == 0) {
                    if(method == "KS") {
                         p <- KS.Diagnostic(x[iter:n])
                         if(p <= 0.05) stationary <- 1}
                    else { #method == Geweke
                         z <- try(Geweke.Diagnostic(x[iter:n]),
                              silent=TRUE)
                         if(inherits(z, "try-error")) z <- 3
                         if(abs(z < 2)) stationary <- 1}
                    if(stationary == 0) iter <- iter + jump
                    if(iter >= n) stationary <- 1}
               if(iter > 1) iter <- iter - 1
               if(iter > n) iter <- n
               burn[i] <- iter
               }
          }
     return(burn)
     }

#End


