###########################################################################
# print.demonoid                                                          #
#                                                                         #
# The purpose of the print.demonoid function is to print the contents of  #
# an object of class demonoid to the screen.                              #
###########################################################################

print.demonoid <- function(x, ...)
     {
     if(missing(x)) stop("The x argument is required.")
     cat("Call:\n")
     print(x$Call)
     cat("\nAcceptance Rate: ", round(x$Acceptance.Rate,5),
          "\n", sep="")
     cat("Algorithm: ", x$Algorithm, "\n", sep="")
     cat("Covariance Matrix: (NOT SHOWN HERE; diagonal shown instead)\n")
     if(is.matrix(x$Covar)) {
          print(diag(x$Covar))
          }
     else if(!is.list(x$Covar) & is.vector(x$Covar)) {
          print(x$Covar)
          }
     else for (i in 1:length(x$Covar)) {
          cat("Block:", i, "\n")
          print(diag(x$Covar[[i]]))
          cat("\n")}
     cat("\nCovariance (Diagonal) History: (NOT SHOWN HERE)\n")
     cat("Deviance Information Criterion (DIC):\n")
     DIC <- matrix(c(round(x$DIC1[1],3), round(x$DIC1[2],3),
          round(x$DIC1[3],3), round(x$DIC2[1],3), round(x$DIC2[2],3),
          round(x$DIC2[3],3)), 3, 2,
          dimnames=list(c("Dbar","pD","DIC"),c("All","Stationary")))
     print(DIC)
     cat("Initial Values:\n")
     print(x$Initial.Values)
     cat("\nIterations: ", x$Iterations, "\n", sep="")
     cat("Log(Marginal Likelihood): ", x$LML, "\n", sep="")
     cat("Minutes of run-time: ", round(x$Minutes,2), "\n",
          sep="")
     cat("Model: (NOT SHOWN HERE)\n")
     cat("Monitor: (NOT SHOWN HERE)\n")
     cat("Parameters (Number of): ", x$Parameters, "\n",
          sep="")
     cat("Posterior1: (NOT SHOWN HERE)\n")
     cat("Posterior2: (NOT SHOWN HERE)\n")
     cat("Recommended Burn-In of Thinned Samples: ",
          x$Rec.BurnIn.Thinned, "\n", sep="")
     cat("Recommended Burn-In of Un-thinned Samples: ",
          x$Rec.BurnIn.UnThinned, "\n", sep="")
     cat("Recommended Thinning: ", x$Rec.Thinning, "\n", sep="")
     cat("Specs: (NOT SHOWN HERE)\n")
     cat("Status is displayed every ", x$Status,
          " iterations\n", sep="")
     cat("Summary1: (SHOWN BELOW)\n")
     cat("Summary2: (SHOWN BELOW)\n")
     cat("Thinned Samples: ", x$Thinned.Samples, "\n",
          sep="")
     cat("Thinning: ", x$Thinning, "\n", sep="")
     cat("\n\nSummary of All Samples\n")
     print(x$Summary1)
     cat("\n\nSummary of Stationary Samples\n")
     print(x$Summary2)
     invisible(x)
     }

#End
