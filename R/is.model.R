###########################################################################
# is.model                                                                #
#                                                                         #
# The purpose of the is.model function is to estimate if a model          #
# specification function meets some minimum criteria.                     #
###########################################################################

is.model <- function(Model, Initial.Values, Data)
     {
     if(missing(Model)) stop("The Model argument is required.")
     ismodel <- TRUE
     if(!is.function(Model)) {
          cat("\nModel must be a function.\n")
          ismodel <- FALSE}
     if(missing(Initial.Values))
          stop("Initial.Values argument is required.")
     if(!is.vector(Initial.Values))
          stop("Initial.Values must be a vector.")
     if(missing(Data)) stop("The Data argument is required.")
     if(!is.data(Data)) stop("The Data argument is not Data.")
     if(!identical(length(Initial.Values), length(Data[["parm.names"]])))
          stop("Lengths of Initial.Values and parm.names differ.")
     Mo <- try(Model(Initial.Values, Data), silent=TRUE)
     if(inherits(Mo, "try-error")) stop("Error in executing the Model.")
     if(!is.list(Mo)) {
          cat("\nModel must return a list.\n")
          ismodel <- FALSE}
     else if(length(Mo) != 5) {
          cat("\nModel must return 5 list components.\n")
          ismodel <- FALSE}
     else if(!identical(Mo[[1]], Mo[["LP"]])) {
          cat("\nThe first output component must be named LP.\n")
          ismodel <- FALSE}
     else if(length(Mo[["LP"]]) != 1) {
          cat("\nThe length of LP must be 1.\n")
          ismodel <- FALSE}
     else if(!identical(Mo[[2]], Mo[["Dev"]])) {
          cat("\nThe second output component must be named Dev.\n")
          ismodel <- FALSE}
     else if(length(Mo[["Dev"]]) != 1) {
          cat("\nThe length of Dev must be 1.\n")
          ismodel <- FALSE}
     else if(!identical(Mo[[3]], Mo[["Monitor"]])) {
          cat("\nThe third output component must be named Monitor.\n")
          ismodel <- FALSE}
     else if(!identical(length(Mo[["Monitor"]]),
          length(Data[["mon.names"]]))) {
          cat("\nThe lengths of Monitor values and mon.names differ.\n")
          ismodel <- FALSE}
     else if(!identical(Mo[[4]], Mo[["yhat"]])) {
          cat("\nThe fourth output component must be named yhat.\n")
          ismodel <- FALSE}
     else if(!identical(Mo[[5]], Mo[["parm"]])) {
          cat("\nThe fifth output component must be named parm.\n")
          ismodel <- FALSE}
     else if(!identical(length(Mo[["parm"]]),
          length(Data[["parm.names"]]))) {
          cat("\nThe lengths of parm and parm.names differ.\n")
          ismodel <- FALSE}
     return(ismodel)
     }

#End
