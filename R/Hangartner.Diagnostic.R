###########################################################################
# Hangartner.Diagnostic                                                   #
###########################################################################

Hangartner.Diagnostic <- function(x, J=2) {
     x <- as.vector(x)
     if(!all(x == round(x))) stop("x is not discrete.")
     N <- length(x)
     j <- rep(1:J, each=N/J)
     if(N %% J != 0) stop("N must be divisible by J.")
     tab <- table(x, j)
     out <- chisq.test(x, j)
     class(out) <- "hangartner"
     return(out)
     }

#End
