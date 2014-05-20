###########################################################################
# ESS                                                                     #
#                                                                         #
# The purpose of the ESS function is to estimate the effective sample     #
# size (ESS) of a target distribution after taking autocorrelation into   #
# account. Although the code is slightly different, it is essentially the #
# same as the effectiveSize function in the coda package.                 #
###########################################################################

ESS <- function(x)
     {
     x <- as.matrix(x)
     v0 <- order <- rep(0, ncol(x))
     names(v0) <- names(order) <- colnames(x)
     N <- nrow(x)
     z <- 1:N
     for (i in 1:ncol(x)) {
          lm.out <- lm(x[, i] ~ z)
          if(!identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
               ar.out <- try(ar(x[,i], aic=TRUE), silent=TRUE)
               if(!inherits(ar.out, "try-error")) {
                    v0[i] <- ar.out$var.pred / {1 - sum(ar.out$ar)}^2
                    order[i] <- ar.out$order}}}
     spec <- list(spec=v0, order=order)
     spec <- spec$spec
     temp <- N * .colVars(x) / spec
     out <- spec
     out[which(spec != 0)] <- temp[which(spec != 0)]
     out[which(out < .Machine$double.eps)] <- .Machine$double.eps
     out[which(out > N)] <- N
     return(out)
     }

#End
