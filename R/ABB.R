###########################################################################
# Approximate Bayesian Bootstrap (ABB)                                    #
#                                                                         #
# The purpose of the ABB function is to perform Multiple Imputation (MI)  #
# with the Approximate Bayesian Bootstrap (ABB).                          #
###########################################################################

ABB <- function(X, K=1)
     {
     ### Initial Checks
     if(missing(X)) stop("X is a required argument.")
     if(!is.matrix(X)) X <- as.matrix(X)
     J <- ncol(X)
     N <- nrow(X)
     ### Missingness Indicator
     M <- X*0
     M[which(is.na(X))] <- 1
     if(sum(M) == 0) stop("There are no missing values to impute.")
     M.sums <- colSums(M)
     ### Approximate Bayesian Bootstrap
     MI <- list()
     for (k in 1:K) {
          imp <- NULL
          for (j in 1:J) {
               if(M.sums[j] > 0) {
                    ### Sample X.star.obs | X.obs
                    X.obs <- X[which(M[,j] == 0),j]
                    X.star.obs <- sample(X.obs, length(X.obs),
                         replace=TRUE)
                    ### Sample X.star.mis | X.star.obs
                    X.star.mis <- sample(X.star.obs, M.sums[j],
                         replace=TRUE)
                    if(length(imp) > 0) imp <- c(imp, X.star.mis)
                    else imp <- X.star.mis}
               }
          MI[[k]] <- imp
          }
     return(MI)
     }

#End
