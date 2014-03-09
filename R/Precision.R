###########################################################################
# Precision                                                               #
#                                                                         #
# The purpose of these functions is to facilitate conversions between the #
# precision, standard deviation, and variance of scalars, vectors, and    #
# matrices.                                                               #
###########################################################################

Cov2Prec <- function(Cov)
     {
     if(any(!is.finite(Cov))) stop("Cov must be finite.")
     if(is.matrix(Cov)) {
          if(!is.positive.definite(Cov))
               stop("Cov is not positive-definite.")
          Prec <- as.inverse(Cov)}
     else if(is.vector(Cov)) {
          k <- as.integer(sqrt(length(Cov)))
          Cov <- matrix(Cov, k, k)
          if(!is.positive.definite(Cov))
               stop("Cov is not positive-definite.")
          Prec <- as.inverse(Cov)}
     return(Prec)
     }
Prec2Cov <- function(Prec)
     {
     if(any(!is.finite(Prec))) stop("Prec must be finite.")
     if(is.matrix(Prec)) {
          if(!is.positive.definite(Prec))
               stop("Prec is not positive-definite.")
          Cov <- as.inverse(Prec)}
     else if(is.vector(Prec)) {
          k <- as.integer(sqrt(length(Prec)))
          Prec <- matrix(Prec, k, k)
          if(!is.positive.definite(Prec))
               stop("Prec is not positive-definite.")
          Cov <- as.inverse(Prec)}
     return(Cov)
     }
prec2sd <- function(prec=1)
     {
     prec <- as.vector(prec)
     if(prec <=0) stop("prec must be positive.")
     return(sqrt(1/prec))
     }
prec2var <- function(prec=1)
     {
     prec <- as.vector(prec)
     if(prec <=0) stop("prec must be positive.")
     return(1/prec)
     }
sd2prec <- function(sd=1)
     {
     sd <- as.vector(sd)
     if(sd <=0) stop("sd must be positive.")
     return(1/sd^2)
     }
sd2var <- function(sd=1)
     {
     sd <- as.vector(sd)
     if(sd <=0) stop("sd must be positive.")
     return(sd^2)
     }
var2prec <- function(var=1)
     {
     var <- as.vector(var)
     if(var <=0) stop("var must be positive.")
     return(1/var)
     }
var2sd <- function(var=1)
     {
     var <- as.vector(var)
     if(var <=0) stop("var must be positive.")
     return(sqrt(var))
     }

#End
