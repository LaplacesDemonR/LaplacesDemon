###########################################################################
# Model.Spec.Time                                                         #
#                                                                         #
# The purpose of the Model.Spec.Time function is to return three things:  #
# the amount of time in minutes that it took to evaluate a model          #
# specification a number of times, the evaluations per minute, and the    #
# componentwise iterations per minute.                                    #
###########################################################################

Model.Spec.Time <- function(Model, Initial.Values, Data, n=1000)
     {
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Initial.Values))
          stop("The Initial.Values argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     t <- as.vector(system.time(for (i in 1:n) {Model(Initial.Values, Data)})[3])
     out <- list(Time=round(t/60,3),
          Evals.per.Minute=round(n/(t/60),3),
          Componentwise.Iters.per.Minute=round(n/(t/60)/length(Initial.Values),3))
     return(out)
     }

#End
