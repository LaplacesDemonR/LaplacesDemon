###########################################################################
# is.constrained                                                          #
#                                                                         #
# The purpose of the is.constrained function is to provide a logical test #
# of whether or not initial values change as they are passed through the  #
# Model specification function, given data.                               #
###########################################################################

is.constrained <- function(Model, Initial.Values, Data)
     {
     if(missing(Model))
          stop("The Model argument is required.")
     if(missing(Initial.Values))
          stop("The Initial.Values argument is required.")
     if(missing(Data))
          stop("The Data argument is required.")
     Mo <- Model(Initial.Values, Data)
     constr <- Initial.Values != Mo[["parm"]]
     return(constr)
     }

#End
