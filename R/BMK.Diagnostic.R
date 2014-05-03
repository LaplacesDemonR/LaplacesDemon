###########################################################################
# BMK.Diagnostic                                                          #
#                                                                         #
# The purpose of the BMK.Diagnostic function is to estimate Hellinger     #
# distances between consecutive batches of posterior samples, so that     #
# when Hellinger distances are below a given threshold the portion of the #
# chain may suggest MCMC convergence. Although the code is slightly       #
# different, it is similar to the bmkconverge function in the BMK         #
# package.                                                                #
###########################################################################

BMK.Diagnostic <- function(X, batches=10)
     {
     HD.Batch <- function(batch1, batch2)
          {
          n1 <- nrow(as.matrix(batch1))
          batches.combined <- c(batch1, batch2)
          batches.min <- min(batches.combined)
          batches.max <- max(batches.combined)
          P1 <- try(density(batch1, from=batches.min, to=batches.max,
               n=n1), silent=TRUE)
          Q1 <- try(density(batch2, from=batches.min, to=batches.max,
               n=n1), silent=TRUE)
          if(inherits(P1, "try-error")) P1 <- density(rnorm(n1,0,1))
          if(inherits(Q1, "try-error")) Q1 <- density(rnorm(n1,1,1))
          step1 <- P1$x[2] - P1$x[1]
          diver1 <- (sqrt(P1$y) - sqrt(Q1$y))^2 * step1
          out <- sqrt(sum(diver1) / 2)
          return(out)
          }
     HD.Diag <- function(x, batch.size, batch.list)
          {
          x <- as.vector(x)
          c1 <- 0
          for (i in 1:(length(batch.list)-1)) {
               batch.label1 <- batch.list[i]:(batch.list[i+1]-1)
               batch.label2 <- batch.list[i+1]:((i+1)*batch.size)
               HD <- try(HD.Batch(x[batch.label1], x[batch.label2]),
                    silent=TRUE)
               if(inherits(HD, "try-error")) HD <- 0
               c1 <- c(c1, HD)}
          c1 <- c1[-1]
          return(c1)
          }
     ### Initial Checks
     if(!is.matrix(X) & !identical(class(X), "demonoid"))
          stop("X must be a matrix or an object of class demonoid.")
     if(identical(class(X), "demonoid")) X <- X$Posterior1
     n.iter <- nrow(X)
     n.par <- ncol(X)
     batch.size <- floor(n.iter / batches)
     if(n.iter %% batch.size != 0)
          stop("Batches of even size are required.")
     batch.list <- seq(from=1, to=n.iter, by=batch.size)
     size <- floor(n.iter / batch.size) - 1
     out <- matrix(0, n.par, size)
     ### Hellinger Distance
     for (i in 1:n.par) {out[i,] <- HD.Diag(X[,i], batch.size, batch.list)}
     ### Constrain to the interval [0,1]
     d <- dim(out)
     out <- as.vector(out)
     out.num <- which(out < 0)
     out[out.num] <- 0
     out.num <- which(out > 1)
     out[out.num] <- 1
     out <- array(out, dim=d)
     ### Output
     if(is.null(colnames(X)))
          rownames(out) <- paste("V", 1:ncol(X), sep="")
     else rownames(out) <- colnames(X)
     colnames(out) <- (1:size)*batch.size
     class(out) <- "bmk"
     return(out)
     }

#End
