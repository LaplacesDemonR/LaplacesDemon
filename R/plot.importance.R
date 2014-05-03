###########################################################################
# plot.importance                                                         #
#                                                                         #
# The purpose of the plot.importance function is to plot variable         #
# importance according either to BPIC or the L-criterion in an object of  #
# class importance.                                                       #
###########################################################################

plot.importance <- function(x, Style="BPIC", ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(!identical(Style, "BPIC") & !identical(Style, "Concordance") &
          !identical(Style, "Discrep") & !identical(Style, "L-criterion"))
          stop("Style is unrecognized.")
     if(!identical(class(x), "importance"))
          stop("x must be of class importance.")
     if(identical(Style, "BPIC"))
          dotchart(x[,1], main="Variable Importance", xlab="BPIC", pch=20)
     else if(identical(Style, "Concordance"))
          dotchart(x[,2], main="Variable Importance", xlab="Concordance",
               pch=20)
     else if(identical(Style, "Discrep"))
          dotchart(x[,3], main="Variable Importance",
               xlab="Discrepancy Statistic", pch=20)
     else dotchart(x[,4], main="Variable Importance", xlab="L-criterion",
               pch=20)
     return(invisible())
     }

#End
