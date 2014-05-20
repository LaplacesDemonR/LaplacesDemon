###########################################################################
# summary.miss                                                            #
#                                                                         #
# The purpose of the summary.miss function is to summarize an object of   #
# class miss.                                                             #
###########################################################################

summary.miss <- function(object=NULL, ...)
     {
     if(is.null(object)) stop("The object argument is NULL.")
     x <- object$Imp
     Summ <- matrix(NA, nrow(x), 7, dimnames=list(1:nrow(x),
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ[,1] <- rowMeans(x)
     Summ[,2] <- sqrt(.rowVars(x))
     Summ[,3] <- 0
     Summ[,4] <- 0
     Summ[,5] <- apply(x, 1, quantile, c(0.025), na.rm=TRUE)
     Summ[,6] <- apply(x, 1, quantile, c(0.500), na.rm=TRUE)
     Summ[,7] <- apply(x, 1, quantile, c(0.925), na.rm=TRUE)
     acf.temp <- matrix(1, trunc(10*log10(ncol(x))), nrow(x))
     for (i in 1:nrow(x)) {
          ### MCSE
          temp <- try(MCSE(x[i,]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ[i,3] <- temp
          else Summ[i,3] <- MCSE(x[i,], method="sample.variance")
          ### ESS
          Summ[i,4] <- ESS(x[i,])}
     print(Summ)
     return(invisible(Summ))
     }

#End
