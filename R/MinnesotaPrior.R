###########################################################################
# MinnesotaPrior                                                          #
#                                                                         #
# The purpose of the MinnesotaPrior function is to return prior           #
# covariance matrices for autoregressive parameters in vector             #
# autoregression (VAR) models.                                            #
###########################################################################

MinnesotaPrior <- function(J, lags=c(1,2), lambda=1, theta=0.5, sigma)
     {
     theta <- max(min(theta, 1), 0)
     Iden <- diag(J)
     L <- length(lags)
     V <- array(0, dim=c(J,J,length(lags)))
     for (l in 1:L) {
          ### Diagonal elements
          V[,,l] <- V[,,l] + Iden * (lambda/lags[l])^2
          ### Off-diagonal elements
          V[,,l] <- V[,,l] + (1 - Iden) *
               ((lambda*theta*matrix(sigma, J, J, byrow=TRUE)) /
               (lags[l]*matrix(sigma, J, J)))^2}
     return(V)
     }

#End
