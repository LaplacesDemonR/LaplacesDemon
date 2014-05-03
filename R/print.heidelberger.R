###########################################################################
# print.heidelberger                                                      #
#                                                                         #
# The purpose of the print.heidelberger function is to print the contents #
# of an object of class raftery to the screen.                            #
###########################################################################

print.heidelberger <- function(x, digits=3, ...) 
     {
     HW.title <- matrix(c("Stationarity", "test", "start", "iteration",
          "p-value", "", "Halfwidth", "test", "Mean", "", "Halfwidth", ""),
          nrow=2)
     y <- matrix("", nrow=nrow(x), ncol=6)
     for (j in 1:ncol(y)) y[,j] <- format(x[,j], digits=digits)
     y[,c(1,4)] <- ifelse(x[,c(1,4)], "passed", "failed")
     y <- rbind(HW.title, y)
     vnames <- if(is.null(rownames(x))) paste("[,", 1:nrow(x), "]", sep="")
     else rownames(x)
     dimnames(y) <- list(c("", "", vnames), rep("", 6))
     print.default(y[, 1:3], quote=FALSE, ...)
     print.default(y[, 4:6], quote=FALSE, ...)
     invisible(x)
     }

#End
