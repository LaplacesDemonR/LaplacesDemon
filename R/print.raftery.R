###########################################################################
# print.raftery                                                           #
#                                                                         #
# The purpose of the print.raftery function is to print the contents of   #
# an object of class raftery to the screen.                               #
###########################################################################

print.raftery <- function(x, digits=3, ...) 
     {
     cat("\nQuantile (q) =", x$params["q"])
     cat("\nAccuracy (r) = +/-", x$params["r"])
     cat("\nProbability (s) =", x$params["s"], "\n")
     if(x$resmatrix[1] == "Error") 
     cat("\nYou need a sample size of at least", x$resmatrix[2], 
          "with these values of q, r and s.\n")
     else {
          out <- x$resmatrix
          for (i in ncol(out)) out[, i] <- format(out[, i], digits=digits)
          out <- rbind(matrix(c("Burn-in ", "Total", "Lower bound ", 
               "Dependence", "(M)", "(N)", "(Nmin)", "factor (I)"), 
               byrow=TRUE, nrow=2), out)
          if(!is.null(rownames(x$resmatrix))) 
          out <- cbind(c("", "", rownames(x$resmatrix)), out)
          dimnames(out) <- list(rep("", nrow(out)), rep("", ncol(out)))
          print.default(out, quote=FALSE, ...)
          cat("\n")}
     invisible(x)
     }

#End
