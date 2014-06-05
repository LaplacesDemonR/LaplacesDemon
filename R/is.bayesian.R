###########################################################################
# is.bayesian                                                             #
#                                                                         #
# The purpose of the is.bayesian function is to determine whether or not  #
# a model is Bayesian by comparing the log-posterior (LP) and the LL.     #
###########################################################################

is.bayesian <- function(Model, Initial.Values, Data)
     {
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Initial.Values))
          stop("The Initial.Values argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     bayesian <- FALSE
     Mo <- Model(Initial.Values, Data)
     LL <- Mo[["Dev"]] / -2
     if(Mo[["LP"]] != LL) bayesian <- TRUE
     return(bayesian)
     }

#End
