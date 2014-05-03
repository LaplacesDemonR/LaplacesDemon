###########################################################################
# print.miss                                                              #
#                                                                         #
# The purpose of the print.miss function is to print the contents of an   #
# object of class miss to the screen.                                     #
###########################################################################

print.miss <- function(x, ...)
     {
     if(missing(x)) stop("The x argument is required.")
     cat("\nAlgorithm:", x$Algorithm)
     cat("\nImp:")
     cat("\n  Missing Values:", nrow(x$Imp))
     cat("\n  Iterations:", ncol(x$Imp))
     cat("\nparm: (NOT SHOWN HERE)")
     cat("\nPostMode: (NOT SHOWN HERE)")
     cat("\nType: (NOT SHOWN HERE)\n")
     invisible(x)
     }

#End
