###########################################################################
# BigData                                                                 #
#                                                                         #
# The purpose of the BigData function is to enable the use of a data set  #
# that is larger than the computer memory (RAM).                          #
###########################################################################

BigData <- function(file, nrow, ncol, size=1, Method="add", CPUs=1,
     Type="PSOCK", FUN, ...)
     {
     FUN <- match.fun(FUN)
     N <- trunc(nrow / size)
     ### Non-Parallel Processing
     if(CPUs == 1) {
          con <- file(file, open="r")
          on.exit(close(con))
          for (i in 1:N) {
               ### Read in a Batch
               X <- matrix(scan(file=con, sep=",", #skip=skip.rows[i],
                    nlines=size, quiet=TRUE), size, ncol, byrow=TRUE)
               ### Perform Function
               if(Method == "rbind") {
                    if(i == 1) out <- FUN(X, ...)
                    else out <- rbind(out, FUN(X, ...))}
               else if(Method == "add") {
                    if(i == 1) out <- FUN(X, ...)
                    else out <- out + FUN(X, ...)}}
          }
     else { ### Parallel Processing
          skip.rows <- c(0, size * 1:(N-1))
          batch <- function(x) {
               #seek(con, 0)
               con <- file(file, open="r")
               on.exit(close(con))
               X <- matrix(scan(file=con, sep=",", skip=skip.rows[x],
                    nlines=size, quiet=TRUE), size, ncol, byrow=TRUE)
               if(Method == "add") out <- sum(FUN(X, ...))
               else out <- FUN(X, ...)
               return(out)
               }
          #library(parallel, quietly=TRUE)
          detectedCores <- max(detectCores(),
               as.integer(Sys.getenv("NSLOTS")), na.rm=TRUE)
          if(CPUs > detectedCores) CPUs <- detectedCores
          cl <- makeCluster(CPUs, Type)
          clusterSetRNGStream(cl)
          out <- parLapply(cl, 1:N, function(x) batch(x))
          stopCluster(cl)
          if(Method == "rbind") {
               out <- unlist(out)
               out <- matrix(out, length(out), 1)
               }
          else if(Method == "add") {
               out <- sum(unlist(out))}
          }
     return(out)
     }

#End
