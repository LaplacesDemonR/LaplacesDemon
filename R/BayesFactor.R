###########################################################################
# BayesFactor                                                             #
#                                                                         #
# The purpose of the BayesFactor function is to estimate a Bayes factor   #
# from two objects, either of class demonoid, laplace, or pmc.            #
###########################################################################

BayesFactor <- function(x)
     {
     ### Initial Checks
     if(missing(x)) stop("x is required.")
     Model.num <- length(x)
     for (i in 1:Model.num) {
          if(!identical(class(x[[i]]), "demonoid") &
               !identical(class(x[[i]]), "laplace") &
               !identical(class(x[[i]]), "pmc") &
               !identical(class(x[[i]]), "vb"))
               stop("x is not of class demonoid, laplace, pmc, or vb.")
          if(identical(class(x[[i]]), "laplace") &
               identical(x[[i]]$Converged, FALSE)) { 
               stop("LaplaceApproximation() did not converge in ",
                    "M[",i,"].\n", sep="")}
          if(identical(class(x[[i]]), "vb") &
               identical(x[[i]]$Converged, FALSE)) { 
               stop("VariationalBayes() did not converge in ",
                    "M[",i,"].\n", sep="")}
          if(is.na(x[[i]]$LML))
               stop(cat("LML is missing in M[",i,"].", sep=""))
          }
     ### Bayes factor
     B <- matrix(NA, Model.num, Model.num)
     for (i in 1:Model.num) {for (j in 1:Model.num) {
          B[i,j] <- exp(x[[i]]$LML - x[[j]]$LML)}}
     strength <- rep(NA,6)
     strength[1] <- "-Inf  <  B <= 0.1   Strong against"
     strength[2] <- "0.1   <  B <= (1/3) Substantial against"
     strength[3] <- "(1/3) <  B < 1      Barely worth mentioning against"
     strength[4] <- "1     <= B < 3      Barely worth mentioning for"
     strength[5] <- "3     <= B < 10     Substantial for"
     strength[6] <- "10    <= B < Inf    Strong for"
     ### Posterior Probability
     ML <- rep(NA, Model.num)
     for (i in 1:Model.num) {ML[i] <- exp(x[[i]]$LML)}
     Posterior.Probability <- ML / sum(ML)
     ### Output
     BF.out <- list(B=B, Hypothesis="row > column", 
          Strength.of.Evidence=strength,
          Posterior.Probability=Posterior.Probability)
     class(BF.out) <- "bayesfactor"
     return(BF.out)
     }

#End
