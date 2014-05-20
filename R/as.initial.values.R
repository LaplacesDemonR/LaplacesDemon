###########################################################################
# as.initial.values                                                       #
#                                                                         #
# The purpose of the as.initial.values function is to retrieve the last   #
# posterior samples from an object of class demonoid, demonoid.hpc,       #
# iterquad, laplace, or pmc to serve as initial values for future         #
# updating.                                                               #
###########################################################################

as.initial.values <- function(x)
     {
     if(!identical(class(x), "demonoid") &
        !identical(class(x), "demonoid.hpc") &
        !identical(class(x), "iterquad") &
        !identical(class(x), "laplace") &
        !identical(class(x), "pmc") &
        !identical(class(x), "vb"))
          stop("The class of x is unknown.")
     if(identical(class(x), "demonoid")) {
          initial.values <- as.vector(x$Posterior1[x$Thinned.Samples,])
          }
     else if(identical(class(x), "demonoid.hpc")) {
          Chains <- length(x)
          LIV <- x[[1]][["Parameters"]]
          initial.values <- matrix(0, Chains, LIV)
          for (i in 1:Chains) {
               initial.values[i,] <- as.vector(x[[i]][["Posterior1"]][x[[i]][["Thinned.Samples"]],])}}
     else if(identical(class(x), "iterquad"))
          initial.values <- as.vector(x$Summary1[,"Mean"])
     else if(identical(class(x), "laplace"))
          initial.values <- as.vector(x$Summary1[,"Mode"])
     else if(identical(class(x), "vb"))
          initial.values <- as.vector(x$Summary1[,"Mean"])
     else if(x$M == 1)
          initial.values <- colMeans(x$Posterior2)
     else if(x$M > 1)
          initial.values <- t(x$Mu[dim(x$Mu)[1],,])
     return(initial.values)
     }

#End
