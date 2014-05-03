###########################################################################
# BayesTheorem                                                            #
###########################################################################

BayesTheorem <- function(PrA, PrBA)
     {
     if(missing(PrA)) stop("The PrA argument is required.")
     if(missing(PrBA)) stop("The PrBA argument is required.")
     if(any(PrA < 0) | any(PrA > 1))
          stop("PrA is not in the interval [0,1].")
     if(any(PrBA < 0) | any(PrBA > 1))
          stop("PrBA is not in the interval [0,1].")
     PrAB <- (PrBA * PrA) / sum(PrBA * PrA)
     class(PrAB) <- "bayestheorem"
     return(PrAB)
     }

#End
