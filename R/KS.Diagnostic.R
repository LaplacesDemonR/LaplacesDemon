###########################################################################
# KS.Diagnostic                                                           #
#                                                                         #
# The purpose of the KS.Diagnostic is to assess the stationarity of a     #
# MCMC chain, given its posterior samples.                                #
###########################################################################

KS.Diagnostic <- function(x)
     {
     if(missing(x)) stop("The x argument is required.")
     if(!is.vector(x)) x <- as.vector(x)
     n <- length(x)
     half <- round(n/2)
     out <- ks.test(x[1:half], x[-c(1:half)])$p.value
     return(out)
     }

#End


