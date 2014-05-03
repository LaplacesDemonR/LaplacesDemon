###########################################################################
# Log-Log                                                                 #
#                                                                         #
# The logit and probit links are symmetric, because the probabilities     #
# approach zero or one at the same rate. The log-log and complementary    #
# log-log links are asymmetric. Complementary log-log links approach zero #
# slowly and one quickly. Log-log links approach zero quickly and one     #
# slowly. Either the log-log or complementary log-log link will tend to   #
# fit better than logistic and probit, and are frequently used when the   #
# probability of an event is small or large. A mixture of the two links,  #
# the log-log and complementary log-log is often used, where each link is #
# weighted. The reason that logit is so prevalent is because logistic     #
# parameters can be interpreted as odds ratios.                           #
###########################################################################

loglog <- function(p)
     {
     if({any(p < 0)} || {any(p > 1)}) stop("p must be in [0,1].")
     x <- log(-log(p))
     return(x)
     }
invloglog <- function(x)
     {
     p <- exp(-exp(x))
     return(p)
     }
cloglog <- function(p)
     {
     if({any(p < 0)} || {any(p > 1)}) stop("p must be in [0,1].")
     x <- log(-log(1 - p))
     return(x)
     }
invcloglog <- function(x)
     {
     p <- 1 - exp(-exp(x))
     return(p)
     }

#End
