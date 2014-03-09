###########################################################################
# Kullback-Leibler Divergence (KLD)                                       #
#                                                                         #
# The purpose of the KLD function is to calculate the Kullback-Leibler    #
# divergences between two probability distributions, p(x) and p(y).       #
###########################################################################

KLD <- function(px, py, base=exp(1))
     {
     ### Initial Checks
     if(!is.vector(px)) px <- as.vector(px)
     if(!is.vector(py)) py <- as.vector(py)
     n1 <- length(px)
     n2 <- length(py)
     if(!identical(n1, n2)) stop("px and py must have the same length.")
     if(any(!is.finite(px)) || any(!is.finite(py))) 
          stop("px and py must have finite values.")
     if(any(px <= 0)) px <- exp(px)
     if(any(py <= 0)) py <- exp(py)
     px[which(px < .Machine$double.xmin)] <- .Machine$double.xmin
     py[which(py < .Machine$double.xmin)] <- .Machine$double.xmin
     ### Normalize
     px <- px / sum(px)
     py <- py / sum(py)
     ### Kullback-Leibler Calculations
     KLD.px.py <- px * (log(px, base=base)-log(py, base=base))
     KLD.py.px <- py * (log(py, base=base)-log(px, base=base))
     sum.KLD.px.py <- sum(KLD.px.py)
     sum.KLD.py.px <- sum(KLD.py.px)
     mean.KLD <- (KLD.px.py + KLD.py.px) / 2
     mean.sum.KLD <- (sum.KLD.px.py + sum.KLD.py.px) / 2
     ### Output
     out <- list(KLD.px.py=KLD.px.py, #KLD[i](p(x[i]) || p(y[i]))
          KLD.py.px=KLD.py.px, #KLD[i](p(y[i]) || p(x[i]))
          mean.KLD=mean.KLD,
          sum.KLD.px.py=sum.KLD.px.py, #KLD(p(x) || p(y))
          sum.KLD.py.px=sum.KLD.py.px, #KLD(p(y) || p(x))
          mean.sum.KLD=mean.sum.KLD,
          intrinsic.discrepancy=min(sum.KLD.px.py, sum.KLD.py.px))
     return(out)
     }

#End
