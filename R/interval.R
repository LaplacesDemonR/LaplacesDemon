###########################################################################
# interval                                                                #
#                                                                         #
# The purpose of the interval function is to constrain the element(s) of  #
# a scalar, vector, matrix, or array to the interval [a,b].               #
###########################################################################

interval <- function(x, a=-Inf, b=Inf, reflect=TRUE)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(a > b) stop("a > b.")
     if(reflect & is.finite(a) & is.finite(b) & any(!is.finite(x))) {
               if(is.array(x)) {
                    d <- dim(x)
                    x <- as.vector(x)}
               x.inf.pos <- !is.finite(x);
               x[x.inf.pos] <- interval(x[x.inf.pos], a, b, reflect=FALSE)
               if(is.array(x)) x <- array(x, dim=d)
          }
     ### Scalar
     if(is.vector(x) & {length(x) == 1}) {
          if(reflect == FALSE) x <- max(a, min(b, x))
          else if(x < a | x > b) {
               out <- TRUE
               while(out) {
                    if(x < a) x <- a + a - x
                    if(x > b) x <- b + b - x
                    if(x >= a & x <= b) out <- FALSE
                    }}}
     ### Vector
     else if(is.vector(x) & {length(x) > 1}) {
          if(reflect == FALSE) {
               x.num <- which(x < a)
               x[x.num] <- a
               x.num <- which(x > b)
               x[x.num] <- b}
          else if(any(x < a) | any(x > b)) {
               out <- TRUE
               while(out) {
                    x.num <- which(x < a)
                    x[x.num] <- a + a - x[x.num]
                    x.num <- which(x > b)
                    x[x.num] <- b + b - x[x.num]
                    if(all(x >= a) & all(x <= b)) out <- FALSE
                    }}}
     ### Matrix or Array
     else if(is.array(x)) {
          d <- dim(x)
          x <- as.vector(x)
          if(reflect == FALSE) {
               x.num <- which(x < a)
               x[x.num] <- a
               x.num <- which(x > b)
               x[x.num] <- b}
          else if(any(x < a) | any(x > b)) {
               out <- TRUE
               while(out) {
                    x.num <- which(x < a)
                    x[x.num] <- a + a - x[x.num]
                    x.num <- which(x > b)
                    x[x.num] <- b + b - x[x.num]
                    if(all(x >= a) & all(x <= b)) out <- FALSE
               }}
          x <- array(x, dim=d)}
     return(x)
     }

#End
