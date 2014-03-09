###########################################################################
# Stick                                                                   #
#                                                                         #
# The purpose of the Stick function is provide the utility of truncated   #
# stick-breaking regarding the vector theta.                              #
###########################################################################

Stick <- function(theta)
     {
     M <- length(theta) + 1
     theta <- c(theta, 1)
     p <- rep(theta[1], length(theta))
     for (m in 1:(M-1)) {
          p[m+1] <- theta[m+1] * (1-theta[m]) * p[m] / theta[m]}
     return(p)
     }

#End
