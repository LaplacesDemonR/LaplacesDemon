###########################################################################
# logit                                                                   #
#                                                                         #
# The logit function is the inverse of the sigmoid or logistic function,  #
# and transforms a continuous value (usually probability p) in the        #
# interval [0,1] to the real line (where it usually is the logarithm of   #
# the odds). The invlogit function (called either the inverse logit or    #
# the logistic function) transforms a real number (usually the logarithm  #
# of the odds) to a value (usually probability p) in the interval [0,1].  #
# If p is a probability, then p/(1-p) is the corresponding odds, while    #
# logit of p is the logarithm of the odds. The difference between the     #
# logits of two probabilities is the logarithm of the odds ratio. The     #
# derivative of probability p in a logistic function is:                  #
# (d / dx) = p * (1 - p).                                                 #
###########################################################################

invlogit <- function(x)
     {
     InvLogit <- 1 / {1 + exp(-x)}
     return(InvLogit)
     }
logit <- function(p)
     {
     if({any(p < 0)} || {any(p > 1)}) stop("p must be in [0,1].")
     Logit <- log(p / {1 - p})
     return(Logit)
     }

#End
