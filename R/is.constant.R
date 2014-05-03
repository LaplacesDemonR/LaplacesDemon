###########################################################################
# is.constant                                                             #
#                                                                         #
# The purpose of the is.constant function is to provide a logical test of #
# whether or not a vector is a constant.                                  #
###########################################################################

is.constant <-  function(x)
     {
     if(missing(x)) stop("The x argument is required.")
     if(!is.vector(x)) x <- as.vector(x)
     uni <- length(unique(x))
     return(uni <= 1)
     }

#End
