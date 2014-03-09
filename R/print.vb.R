###########################################################################
# print.vb                                                                #
#                                                                         #
# The purpose of the print.vb function is to print the contents of an     #
# object of class vb to the screen.                                       #
###########################################################################

print.vb <- function(x, ...)
     {
     if(missing(x)) stop("The x argument is required.")
     cat("\nCall:\n")
     print(x$Call)
     cat("\nConverged: ", x$Converged, "\n", sep="")
     cat("Covariance Matrix: (NOT SHOWN HERE; diagonal shown instead)\n")
     print(diag(x$Covar))
     cat("\nDeviance (Final): ", x$Deviance[length(x$Deviance)], "\n")
     cat("History: (NOT SHOWN HERE)\n")
     cat("Initial Values:\n")
     print(x$Initial.Values)
     cat("\nIterations: ", x$Iterations, "\n", sep="")
     cat("Log(Marginal Likelihood): ", x$LML, "\n", sep="")
     cat("Log-Posterior (Final): ", x$LP.Final, "\n", sep="")
     cat("Log-Posterior (Initial): ", x$LP.Initial, "\n", sep="")
     cat("Minutes of run-time: ", x$Minutes, "\n", sep="")
     cat("Monitor: (NOT SHOWN HERE)\n")
     cat("Posterior: (NOT SHOWN HERE)\n")
     cat("Step Size (Final): ")
     print(x$Step.Size.Final)
     cat("Step Size (Initial): ", x$Step.Size.Initial, "\n", sep="")
     cat("Summary1: (SHOWN BELOW)\n")
     cat("Summary2: (SHOWN BELOW)\n")
     cat("Tolerance (Final): ", x$Tolerance.Final, "\n", sep="")
     cat("Tolerance (Stop): ", x$Tolerance.Stop, "\n", sep="")
     
     cat("\nSummary1:\n")
     print(x$Summary1)
     if({x$Converged == TRUE} && !any(is.na(x$Posterior))) {
          cat("\nSummary2:\n")
          print(x$Summary2)}
     
     invisible(x)
     }

#End
