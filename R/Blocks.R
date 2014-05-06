###########################################################################
# Blocks                                                                  #
#                                                                         #
# The purpose of the Blocks function is to return a list in which each    #
# component is a block and contains a vector of positions that indicate   #
# parameter membership.                                                   #
###########################################################################

Blocks <- function(Initial.Values, N, PostCor=NULL)
     {
     ### Initial Checks
     if(missing(Initial.Values))
          stop("The Initial.Values argument is required.")
     Initial.Values <- as.vector(Initial.Values)
     LIV <- length(Initial.Values)
     if(LIV <= 1) stop("More Initial.Values are needed to create blocks.")
     if(missing(N)) N <- trunc(sqrt(LIV))
     if(is.null(PostCor)) {
          if(length(N != 1)) N <- N[1]
          N <- max(min(round(N), LIV), 2)
          }
     if(!is.null(PostCor)) {
          if(!is.matrix(PostCor)) stop("PostCor must be a matrix.")
          if(!identical(nrow(PostCor), ncol(PostCor)))
               stop("PostCor must be a square matrix.")
          if(any(diag(PostCor) != 1))
               stop("PostCor requires the diagonal to contain 1s.")
          if(length(N) == 1) N <- max(min(round(N), LIV), 2)
          else {
               N <- N[1:2]
               N[1] <- max(min(round(N[1]), LIV), 2)
               N[2] <- max(min(round(N[2]), LIV), 2)
               if(N[1] >= N[2]) stop("N is incorrect.")
               }
          }
     ### Silhouette (from silhouette.default.R in package cluster)
     silhouette <- function(x, dist, dmatrix, ...)
          {
          cll <- match.call()
          if(is.list(x) && !is.null(cl <- x$clustering)) x <- cl
          n <- length(x)
          if(!all(x == round(x))) stop("'x' must only have integer codes")
          k <- length(clid <- sort(unique(x)))
          if(k <= 1 || k >= n) return(NA)
          if(missing(dist)) {
               if(missing(dmatrix))
                    stop("Need either a dissimilarity 'dist' or diss.matrix 'dmatrix'")
               if(is.null(dm <- dim(dmatrix)) || length(dm) != 2 || !all(n == dm))
                    stop("'dmatrix' is not a dissimilarity matrix compatible to 'x'")
               }
          else {
               dist <- as.dist(dist)
               if(n != attr(dist, "Size"))
                    stop("clustering 'x' and dissimilarity 'dist' are incompatible")
               dmatrix <- as.matrix(dist)}
          wds <- matrix(NA, n,3, dimnames =
               list(names(x), c("cluster","neighbor","sil_width")))
          for (j in 1:k) {
               Nj <- sum(iC <- x == clid[j])
               wds[iC, "cluster"] <- clid[j]
               diC <- rbind(apply(dmatrix[!iC, iC, drop=FALSE], 2,
                    function(r) tapply(r, x[!iC], mean)))
               minC <- apply(diC, 2, which.min)
               wds[iC,"neighbor"] <- clid[-j][minC]
               s.i <- if(Nj > 1) {
                    a.i <- colSums(dmatrix[iC, iC])/(Nj - 1)
                    b.i <- diC[cbind(minC, seq(along=minC))]
                    ifelse(a.i != b.i, (b.i - a.i) / pmax(b.i, a.i), 0)
                    }
               else 0
               wds[iC,"sil_width"] <- s.i}
          return(wds)
          }
     ### Create Blocks
     if(missing(PostCor)) {
          ### Sequential Blocks
          B <- list()
          pos <- 0
          for (i in 1:N) {
               if(i != N) {
                    B[[i]] <- pos + c(1:trunc(LIV / N))
                    pos <- pos + trunc(LIV / N)
                    }
               else B[[i]] <- c((pos + 1):LIV)}
          }
     else {
          ### Hierarchical Clustering
          di <- dist(1-abs(PostCor))
          hc <- hclust(di, "ave")
          if(length(N) == 1) clusters <- as.vector(cutree(hc, N))
          else {
               av.width <- rep(0, (N[2]-N[1])+1)
               names(av.width) <- seq(from=N[1], to=N[2])
               count <- 1
               for (i in N[1]:N[2]) {
                    av.width[count] <- mean(silhouette(cutree(hc, i), di)[,3])
                    count <- count + 1}
               cat("\nMean Silhouette Width per Hierarchical Cluster Solution:\n")
               print(av.width)
               N <- N[1] - 1 + which.max(av.width)[1]
               cat("\n", N, "Blocks will be created.\n")
               clusters <- as.vector(cutree(hc, N))
               }
          B <- list()
          for (i in 1:N) B[[i]] <- which(clusters == i)
          }
     class(B) <- "blocks"
     return(B)
     }

#End
