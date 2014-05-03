###########################################################################
# AcceptanceRate                                                          #
#                                                                         #
# The purpose of the AcceptanceRate function is to calculate the          #
# acceptance rate of each chain from its samples.                         #
###########################################################################

AcceptanceRate <- function(x)
     {
     if(missing(x)) stop("x is a required argument.")
     if(!is.matrix(x)) x <- as.matrix(x)
     out <- colMeans(x[-nrow(x),] != x[-1,])
     names(out) <- colnames(x)
     return(out)
     }

x <- matrix(rnorm(10*10),10,10)
colnames(x) <- paste("V", 1:10, sep="")
x[2,] <- x[1,]
AcceptanceRate(x)

#End
