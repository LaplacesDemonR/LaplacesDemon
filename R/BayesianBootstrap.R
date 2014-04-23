###########################################################################
# BayesianBootstrap                                                       #
#                                                                         #
# The purpose of the BayesianBootstrap is to allow the user to produce    #
# either bootstrapped weights or statistics.                              #
###########################################################################

BayesianBootstrap <- function(X, n=1000, Method="weights", Status=NULL)
     {
     ### Initial Checks
     if(missing(X)) stop("X is a required argument.")
     if(!is.matrix(X)) X <- as.matrix(X)
     if(any(!is.finite(X))) stop("Non-finite values found in X.")
     S <- round(abs(n))
     if(S < 1) S <- 1
     if(!(is.numeric(Status) & (length(Status) == 1))) Status <- S + 1
     else {
          Status <- round(abs(Status))
          if(Status < 1 | Status > S) Status <- S + 1}
     N <- nrow(X)
     J <- ncol(X)
     if(identical(Method, "weights")) {
          BB <- replicate(S, diff(c(0, sort(runif(N-1)), 1)))
          return(BB)}
     ### Bayesian Bootstrap: Statistics
     BB <- vector("list", S)
     for (s in 1:S) {
          if(s %% Status == 0) cat("\nBootstrapped Samples:", s)
          u <- c(0, sort(runif(N - 1)), 1)
          g <- diff(u)
          BB[[s]] <- Method(X, g)}
     if(Status < S) cat("\n\nThe Bayesian Bootstrap has finished.\n\n")
     ### Output
     BB <- lapply(BB, identity)
     if(is.vector(BB[[1]])) 
          if(length(BB[[1]]) == 1) BB <- as.matrix(BB)
          else {
               B <- matrix(unlist(BB), S, length(BB[[1]]), byrow=TRUE)
               colnames(B) <- names(BB[[1]])
               BB <- B
               }
     else {
          if(is.null(dim(BB[[1]])))
               stop("Method must return a vector, matrix or array")
          B <- array(NA, dim=c(S, dim(BB[[1]])))
          for (s in 1:S) {B[s,,] <- BB[[s]]}
          BB <- B
          }
     return(BB)
     }

#End
