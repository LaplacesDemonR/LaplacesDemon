###########################################################################
# WAIC                                                                    #
#                                                                         #
# The purpose of the WAIC function is to calculate the Widely Applicable  #
# Information Criterion.                                                  #
###########################################################################

WAIC <- function(x)
     {
     lppd <- sum (log(rowMeans(exp(x))))
     pWAIC1 <- 2*sum(log(rowMeans(exp(x))) - rowMeans(x))
     pWAIC2 <- sum(.rowVars(x))
     WAIC <- -2*lppd + 2*pWAIC2
     return(list(WAIC=WAIC, lppd=lppd, pWAIC=pWAIC2, pWAIC1=pWAIC1))
     }

#End
