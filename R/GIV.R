###########################################################################
# GIV                                                                     #
#                                                                         #
# GIV stands for ``generate initial values'', and the purpose of the GIV  #
# function is to generate initial values for the IterativeQuadrature,     #
# LaplaceApproximation, LaplacesDemon, PMC, and VariationalBayes          #
# functions.                                                              #
###########################################################################

GIV <- function(Model, Data, n=1000, PGF=FALSE)
     {
     if(missing(Model)) stop("The Model argument is required.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data)) stop("The Data argument is required.")
     if(!is.list(Data)) stop("Data must be a list.")
     if(is.null(Data[["parm.names"]])) stop("parm.names missing in Data.")
     LIV <- length(Data[["parm.names"]])
     iv <- rep(NA, LIV)
     if(PGF == TRUE) {
          if(is.null(Data[["PGF"]])) stop("PGF missing in Data.")
          for (i in 1:n) {
               IV <- as.vector(Data$PGF(Data))
               M <- try(Model(IV, Data), silent=TRUE)
               if(!inherits(M, "try-error") & is.finite(M[["LP"]]) &
                    is.finite(M[["Dev"]]) & 
                    identical(as.vector(M[["parm"]]), as.vector(IV))) {
                    iv <- IV; break}
               }
          }
     else if(PGF == FALSE) {
          high <- 100; low <- -100
          a <- try(Model(rep(low, LIV), Data)[["parm"]], silent=TRUE)
          b <- try(Model(rep(high, LIV), Data)[["parm"]], silent=TRUE)
          if(inherits(a, "try-error") | inherits(b, "try-error")) {
               for (i in 1:n) {
                    IV <- rnorm(LIV, runif(1,-100,100), runif(1,0.1,1000))
                    M <- try(Model(IV, Data), silent=TRUE)
                    if(!inherits(M, "try-error") & is.finite(M[["LP"]]) &
                         is.finite(M[["Dev"]]) & 
                         identical(as.vector(M[["parm"]]), as.vector(IV))) {
                         iv <- IV; break}
                    }
               }
          else {
               ab.range <- b - a
               ab.mu <- a + ab.range / 2
               ab.mu <- ifelse({a == 0} & {b == high}, 10, ab.mu)
               Scale <- 1 / ab.range
               Scale[which(ab.range == 0)] <- 0
               for (i in 1:n) {
                    IV <- rnorm(LIV, ab.mu, ab.range * Scale)
                    M <- try(Model(IV, Data), silent=TRUE)
                    if(!inherits(M, "try-error") & is.finite(M[["LP"]]) &
                         is.finite(M[["Dev"]]) & 
                         identical(as.vector(M[["parm"]]), as.vector(IV))) {
                         iv <- IV; break}
                    Scale <- Scale + ab.range / n / 2}
               }
          }
     if((i == n) | any(is.na(iv)))
          cat("\nWARNING: Acceptable initial values were not generated.\n")
     return(iv)
     }
#End
