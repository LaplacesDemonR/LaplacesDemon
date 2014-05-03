###########################################################################
# Thin                                                                    #
#                                                                         #
# The purpose of the Thin function is to facilitate the thinning of a     #
# matrix of posterior samples.                                            #
###########################################################################

Thin <- function(x, By=1)
     {
     ### Initial Checks
     if(!is.matrix(x)) x <- as.matrix(x)
     rownum <- nrow(x)
     By <- abs(round(By))
     if(By > rownum) stop("By exceeds number of rows in x.")
     ### Thin
     keeprows <- which(rep(1:By, len=rownum) == By)
     z <- x[keeprows,]
     return(z)
     }

#End
