###########################################################################
# LossMatrix                                                              #
#                                                                         #
# The purpose of the LossMatrix function is facilitate Bayesian decision  #
# theory for discrete actions among discrete states.                      #
###########################################################################

LossMatrix <- function(L, p.theta)
     {
     ### Initial Checks
     if(missing(L)) stop("L is a required argument.")
     if(missing(p.theta)) stop("p.theta is a required argument.")
     if(!is.matrix(L) & !is.array(L)) stop("L must be a matrix or array.")
     if(!is.array(p.theta))
          stop("p.theta must be a vector or matrix.")
     d.L <- dim(L)
     d.p.theta <- dim(p.theta)
     if(any(d.L[1:2] != d.p.theta[1:2]))
          stop("The rows or columns of L and p.theta differ.")
     if(length(d.L) == 3 & length(d.p.theta) == 3)
          if(d.L[3] != d.p.theta[3])
               stop("The number of samples in L and p.theta differ.")
     if(length(d.L) > 3 | length(d.p.theta) > 3)
          stop("L and p.theta may have no more than 3 dimensions.")
     if(length(d.p.theta) == 2) {
          if(any(colSums(p.theta) != 1))
               stop("Each column in p.theta must sum to one.")
          }
     else {
          for (i in 1:d.p.theta[3]) {
               if(any(colSums(p.theta[,,i]) != 1))
                    stop("Each column in p.theta must sum to one.")}}
     ### Expected Loss
     if(length(d.L) == 2 & length(d.p.theta) == 2) {
          E.Loss <- colSums(L * p.theta)
          }
     else if(length(d.L) == 3 & length(d.p.theta) == 2) {
          E.Loss <- rep(0, d.L[3])
          for (i in 1:d.L[3]) {
               E.Loss <- E.Loss + colSums(L[,,i] * p.theta)}
          E.Loss <- E.Loss / d.L[3]
          }
     else if(length(d.L) == 2 & length(d.p.theta) == 3) {
          E.Loss <- rep(0, d.p.theta[3])
          for (i in 1:d.p.theta[3]) {
               E.Loss <- E.Loss + colSums(L * p.theta[,,i])}
          E.Loss <- E.Loss / d.p.theta[3]
          }
     else {
          E.Loss <- rep(0, d.L[3])
          for (i in 1:d.L[3]) {
               E.Loss <- E.Loss + colSums(L[,,i] * p.theta[,,i])}}
     action <- which.min(E.Loss)
     cat("\nAction", action, "minimizes expected loss.\n")
     out <- list(BayesAction=action, E.Loss=E.Loss)
     return(out)
     }

#End
