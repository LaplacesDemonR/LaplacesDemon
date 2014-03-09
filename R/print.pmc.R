###########################################################################
# print.pmc                                                               #
#                                                                         #
# The purpose of the print.pmc function is to print the contents of an    #
# object of class pmc to the screen.                                      #
###########################################################################

print.pmc <- function(x, ...)
     {
     if(missing(x)) stop("The x argument is required.")
     cat("Call:\n")
     print(x$Call)
     cat("\nalpha:\n", sep="")
     print(x$alpha)
     cat("Covariance Matrix: (NOT SHOWN HERE)\n")
     cat("Deviance: (NOT SHOWN HERE)\n")
     cat("Deviance Information Criterion (DIC):\n")
     DIC <- matrix(c(round(x$DIC[1],3), round(x$DIC[2],3),
          round(x$DIC[3],3)), 3, 1,
          dimnames=list(c("Dbar","pD","DIC"),c("All")))
     print(DIC)
     cat("ESSN:\n")
     print(x$ESSN)
     cat("Initial Values:\n")
     print(x$Initial.Values)
     cat("\nIterations: ", x$Iterations, "\n", sep="")
     cat("Log(Marginal Likelihood): ", x$LML, "\n", sep="")
     cat("M (Mixture Components): ", x$M, "\n", sep="")
     cat("Minutes of run-time: ", round(x$Minutes,2), "\n",
          sep="")
     cat("Model: (NOT SHOWN HERE)\n")
     cat("Monitor: (NOT SHOWN HERE)\n")
     cat("Mu: (NOT SHOWN HERE)\n")
     cat("Number of Samples: ", x$N, "\n", sep="")
     cat("nu: ", x$nu, "\n", sep="")
     cat("Parameters (Number of): ", x$Parameters, "\n",
          sep="")
     cat("Perpexity, Normalized:\n")
     print(x$Perplexity)
     cat("Posterior1: (NOT SHOWN HERE)\n")
     cat("Posterior2: (NOT SHOWN HERE)\n")
     cat("Summary: (SHOWN BELOW)\n")
     cat("Thinned Samples: ", x$Thinned.Samples, "\n",
          sep="")
     cat("Thinning: ", x$Thinning, "\n", sep="")
     cat("Weights: (NOT SHOWN HERE)\n")
     cat("\n\nSummary:\n")
     print(x$Summary)
     invisible(x)
     }

#End
