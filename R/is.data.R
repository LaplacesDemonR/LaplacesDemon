###########################################################################
# is.data                                                                 #
#                                                                         #
# The purpose of the is.data function is to estimate if a list of data is #
# data as far as IterativeQuadrature, LaplaceApproximation,               #
# LaplacesDemon, PMC, and VariationalBayes are concerned.                 #
###########################################################################

is.data <- function(Data)
     {
     if(missing(Data))
          stop("The Data argument is required.")
     isdata <- TRUE
     if(!is.list(Data)) {
          cat("\nData must be a list.\n")
          isdata <- FALSE}
     if(is.null(Data[["mon.names"]])) {
          cat("\nmon.names is NULL.\n")
          isdata <- FALSE}
     if(is.null(Data[["parm.names"]])) {
          cat("\nparm.names is NULL.\n")
          isdata <- FALSE}
     return(isdata)
     }

#End
