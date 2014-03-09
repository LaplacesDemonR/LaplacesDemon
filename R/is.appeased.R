###########################################################################
# is.appeased                                                             #
#                                                                         #
# The purpose of the is.appeased function is to perform a logical test of #
# whether or not Laplace's Demon is appeased with an object of class      #
# demonoid.                                                               #
###########################################################################

is.appeased <- function(x)
     {
     appeased <- FALSE
     if(!identical(class(x), "demonoid"))
          stop("x must be of class demonoid.")
     captive <- capture.output(Consort(x))
     z <- grep("has been appeased", captive)
     if(length(z) > 0) appeased <- TRUE
     return(appeased)
     }

#End
