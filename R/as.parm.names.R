###########################################################################
# as.parm.names                                                           #
#                                                                         #
# The purpose of the as.parm.names function is to create the vector of    #
# names for the parameters from a list of parameters, which may be any    #
# combination of scalars, vectors, matrices, and upper-triangular         #
# matrices.                                                               #
###########################################################################

as.parm.names <- function(x, uppertri=NULL)
     {
     ### Initial Checks
     if(missing(x)) stop("x is required.")
     if(!is.list(x)) stop("x must be a list.")
     parm.length <- length(x)
     if(is.null(uppertri)) uppertri <- rep(0, parm.length)
     if(!identical(length(uppertri), parm.length))
          stop("Length of uppertri and list attributes differs.")
     ### Length of parm.names
     totlen <- 0
     for (i in 1:parm.length) {
          if(uppertri[i] == 0) {
               xlen <- length(as.vector(which(!is.na(x[[i]]))))
               totlen <- totlen + xlen
               }
          else if(uppertri[i] == 1) {
               if(is.vector(x[[i]])) stop("uppertri=1 found for a vector.")
               xlen <- length(which(!is.na(x[[i]][upper.tri(x[[i]],
                    diag=TRUE)])))
               totlen <- totlen + xlen}}
     ### Assign parm.names
     parm.names <- rep(NA, totlen)
     cnt <- 1
     for (i in 1:parm.length) {
          xname <- names(x)[i]
          xlen <- length(as.vector(x[[i]]))
          ### Scalar
          if(xlen == 1) {
               if(is.na(x[[i]])) stop("scalar has NA.")
               parm.names[cnt] <- paste(xname)
               cnt <- cnt + 1
               }
          ### Vector
          else if(is.vector(x[[i]]) & {xlen > 1}) {
               for (j in which(!is.na(x[[i]]))) {
                    parm.names[cnt] <- paste(xname, "[", j, "]", sep="")
                    cnt <- cnt + 1}
               }
          ### Matrix
          else if(is.matrix(x[[i]]) & (uppertri[i] == 0)) {
               for (k in 1:ncol(x[[i]])) {for (j in 1:nrow(x[[i]])) {
                    if(!is.na(x[[i]][j,k])) {
                         parm.names[cnt] <- paste(xname, "[", j, ",",
                              k, "]", sep = "")
                         cnt <- cnt + 1}}}
               }
          ### Matrix, Upper Triangular
          else if(is.matrix(x[[i]]) & (uppertri[i] == 1)) {
               nr <- nrow(x[[i]])
               nc <- ncol(x[[i]])
               U <- upper.tri(x[[i]], diag=TRUE)
               U[which(is.na(matrix(x[[i]], nr, nc)))] <- FALSE
               for (k in 1:nc) {for (j in 1:nr) {
                    if(U[j, k] == TRUE) {
                         parm.names[cnt] <- paste(xname, "[", j, ",",
                              k, "]", sep = "")
                         cnt <- cnt + 1}}}
               }
          ### Array
          else if(is.array(x[[i]]) & (uppertri[i] == 0)) {
               arrayx <- array(1:prod(dim(x[[i]])), dim(x[[i]]))
               for (j in 1:prod(dim(x[[i]]))) {
                    position <- which(arrayx == j, arr.ind=TRUE)
                    if(!is.na(x[[i]][position])) {
                         parm.names[cnt] <- paste(xname, "[",
                              paste(as.vector(position), sep="",
                              collapse=","), "]", sep="")
                         cnt <- cnt + 1}}
               }
          ### Array, Upper Triangular
          else if(is.array(x[[i]]) & (uppertri[i] == 1))
               stop("upper.tri does not function with arrays.")
          }
     return(parm.names)
     }

#End
