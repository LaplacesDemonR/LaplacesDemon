###########################################################################
# p.interval                                                              #
#                                                                         #
# The purpose of the p.interval function is to estimate either a          #
# quantile-based probability interval, uni-modal highest posterior        #
# density (HPD) interval, or multi-modal HPD intervals of posterior       #
# samples. The code for the uni-modal HPD interval is similar to code in  #
# the coda package, but this is designed to work with matrices or objects #
# of class demonoid, iterquad, laplace, or pmc.                           #
###########################################################################

p.interval <- function(obj, HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE,
     PDF=FALSE, ...)
     {
     ### Initial Checks
     if(missing(obj)) stop("The obj argument is required.")
     if(length(prob) > 1)
          stop("The prob argument must be a scalar.")
     if((prob <= 0) | (prob >= 1))
          stop("The prob argument must be in the interval (0,1).")
     if(identical(class(obj), "demonoid")) {
          thin <- obj$Rec.BurnIn.Thinned
          n <- nrow(obj$Posterior1)
          if(thin < nrow(obj$Posterior1))
               x <- cbind(obj$Posterior1[thin:n,], obj$Monitor[thin:n,])
          if(thin >= nrow(obj$Posterior1))
               x <- cbind(obj$Posterior1, obj$Monitor)
          colnames(x) <- c(colnames(obj$Posterior1),
               colnames(obj$Monitor))
          obj <- x
          }
     else if(identical(class(obj), "iterquad")) {
          if(any(is.na(obj$Posterior)))
               stop("Posterior samples do not exist in obj.")
          obj <- obj$Posterior
          }
     else if(identical(class(obj), "laplace")) {
          if(any(is.na(obj$Posterior)))
               stop("Posterior samples do not exist in obj.")
          obj <- obj$Posterior
          }
     else if(identical(class(obj), "pmc")) {
          x <- cbind(obj$Posterior2, obj$Monitor)
          colnames(x) <- c(colnames(obj$Posterior2), colnames(obj$Monitor))
          obj <- x
          }
     else {
          x <- as.matrix(obj)
          oname <- deparse(substitute(obj))
          colnames(x) <- colnames(obj)#oname
          obj <- x
          }
     if(any(!is.finite(obj)))
          stop("The obj argument must contain finite values.")
     if(any(apply(obj, 2, is.constant)))
          stop("The obj argument has insufficient length.")
     ### Probability Intervals
     if(HPD == TRUE) {
          vals <- apply(obj, 2, sort)
          if(!is.matrix(vals)) stop("obj must have nsamp > 1.")
          nsamp <- nrow(vals)
          npar <- ncol(vals)
          gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
          init <- 1:(nsamp - gap)
          inds <- apply(vals[init + gap, , drop=FALSE] -
               vals[init, , drop=FALSE], 2, which.min)
          ans <- cbind(vals[cbind(inds, 1:npar)],
               vals[cbind(inds + gap, 1:npar)])
          dimnames(ans) <- list(colnames(obj), c("Lower", "Upper"))
          attr(ans, "Probability.Interval") <- gap / nsamp
          }
     if(HPD == FALSE) {
          ans <- apply(obj, 2, quantile,
               probs=c((1-prob)/2,1-((1-prob)/2)))
          ans <- as.matrix(t(ans))
          colnames(ans) <- c("Lower","Upper")
          attr(ans, "Probability.Interval") <- prob}
     if(MM == TRUE) {
          ansmm <- ans
          mm <- apply(obj, 2, is.multimodal)
          if(any(mm)) {
               cat("\n\nPotentially multimodal column vectors:\n",
                    which(mm),"\n")
               vals <- apply(obj, 2, sort)
               if(!is.matrix(vals)) stop("obj must have nsamp > 1.")
               for (m in which(mm)) {
                    kde <- density(vals[,m])
                    dens <- approx(kde$x, kde$y, vals[,m])$y
                    dens.ind <- dens >= as.vector(quantile(dens,
                         probs=1-prob)) * 1
                    ints <- ""
                    count <- 1
                    for (i in 1:nrow(vals)) {
                         if((i == 1) & (dens.ind[i] == 1)) {
                              ints <- paste("(",round(vals[i,m],3),",",sep="")
                              if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
                              ansmm[m,count] <- vals[i,m]
                              count <- count + 1
                              }
                         if(i > 1) {
                              if((dens.ind[i] == 0) & (dens.ind[i-1] == 1)) {
                                   ints <- paste(ints,round(vals[i-1,m],3),")",sep="")
                                   if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
                                   ansmm[m,count] <- vals[i-1,m]
                                   count <- count + 1
                                   }
                              if((dens.ind[i] == 1) & (dens.ind[i-1] == 0))  {
                                   ints <- paste(ints," (",round(vals[i,m],3),",",sep="")
                                   if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
                                   ansmm[m,count] <- vals[i,m]
                                   count <- count + 1
                                   }
                              }
                         }
                    if((dens.ind[i] == 1) & (dens.ind[i-1] == 1)) {
                         ints <- paste(ints,round(vals[i,m],3),")",sep="")
                         if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
                              ansmm[m,count] <- vals[i,m]
                              count <- count + 1
                         }
                    cat("\nColumn", m, "multimodal intervals:", ints, "\n")}
               }
          }
     if(plot == TRUE) {
          if(PDF == TRUE) {
               pdf("P.Interval.Plots.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(MM == FALSE) {
               for (i in 1:nrow(ans)) {
                    kde <- kde.low <- kde.high <- density(obj[,i])
                    kde.low$x <- kde$x[kde$x < ans[i,1]]
                    kde.low$y <- kde$y[which(kde$x < ans[i,1])]
                    kde.high$x <- kde$x[kde$x > ans[i,2]]
                    kde.high$y <- kde$y[which(kde$x > ans[i,2])]
                    plot(kde, xlab="Value", main=colnames(obj)[i])
                    polygon(kde, col="black", border="black")
                    polygon(c(min(kde.low$x), kde.low$x, max(kde.low$x)),
                         c(min(kde.low$y), kde.low$y, min(kde.low$y)),
                         col="gray", border="gray")
                    polygon(c(min(kde.high$x), kde.high$x, max(kde.high$x)),
                         c(min(kde.high$y), kde.high$y, min(kde.high$y)),
                         col="gray", border="gray")
                    abline(v=0, col="red", lty=2)
                    }
               }
          if(MM == TRUE) {
               for (i in 1:nrow(ansmm)) {
                    ### Mode 1
                    kde <- density(obj[,i])
                    x1 <- min(which(kde$x >= ansmm[i,1]))
                    x2 <- max(which(kde$x <= ansmm[i,2]))
                    plot(kde, xlab="Value", main=colnames(obj)[i])
                    polygon(kde, col="gray", border="gray")
                    with(kde, polygon(x=c(x[c(x1,x1:x2,x2)]),
                         y=c(0,y[x1:x2],0), col="black"))
                    ### Mode 2
                    if((ncol(ansmm) > 2) && !is.na(ansmm[i,3])) {
                         x1 <- min(which(kde$x >= ansmm[i,3]))
                         x2 <- max(which(kde$x <= ansmm[i,4]))
                         with(kde, polygon(x=c(x[c(x1,x1:x2,x2)]),
                              y=c(0,y[x1:x2],0), col="black"))
                         }
                    if((ncol(ansmm) > 4) && !is.na(ansmm[i,5])) {
                         x1 <- min(which(kde$x >= ansmm[i,5]))
                         x2 <- max(which(kde$x <= ansmm[i,6]))
                         with(kde, polygon(x=c(x[c(x1,x1:x2,x2)]),
                              y=c(0,y[x1:x2],0), col="black"))
                         }
                    abline(v=0, col="red", lty=2)
                    }
               }
          if(PDF == TRUE) dev.off()
          }
     return(ans)
     }

#End
