###########################################################################
# as.ppc                                                                  #
#                                                                         #
# The purpose of the as.ppc function is to convert an object of class     #
# demonoid.val to an object of class demonoid.ppc, after which the object #
# is ready for posterior predictive checks.                               #
###########################################################################

as.ppc <- function(x, set=3)
     {
     ### Initial Checks
     if(missing(x)) stop("x is required.")
     if(!identical(class(x), "demonoid.val"))
          stop("x is not of class demonoid.val.")
     set <- round(abs(set))
     if(set < 1) set <- 1
     else if(set > 3) set <- 3
     ### ppc
     if(set == 1) ppc <- list(y=x[[1]]$y, yhat=x[[1]]$yhat)
     else if(set == 2) ppc <- list(y=x[[2]]$y, yhat=x[[2]]$yhat)
     else ppc <- list(y=c(x[[1]]$y, x[[2]]$y), yhat=rbind(x[[1]]$yhat,
          x[[2]]$yhat))
     class(ppc) <- "demonoid.ppc"
     return(ppc)
     }

#End
