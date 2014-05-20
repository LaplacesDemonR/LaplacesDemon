###########################################################################
# SensitivityAnalysis                                                     #
#                                                                         #
# The purpose of the SensitivityAnalysis function is to perform a         #
# sensitivity analysis on the posterior distributions and posterior       #
# inferences from two models.                                             #
###########################################################################

SensitivityAnalysis <- function(Fit1, Fit2, Pred1, Pred2)
     {
     ### Initial Checks
     if(missing(Fit1)) stop("Fit1 is a required argument.")
     if(missing(Fit2)) stop("Fit2 is a required argument.")
     if(missing(Pred1)) stop("Pred1 is a required argument.")
     if(missing(Pred2)) stop("Pred2 is a required argument.")
     "%!in%" <- function(x,table) match(x, table, nomatch=0) == 0
     if(class(Fit1) %!in% c("demonoid","iterquad","laplace","pmc","vb"))
          stop("Fit1 is not an object of class demonoid, iterquad, laplace, pmc, or vb.")
     if(class(Fit2) %!in% c("demonoid","iterquad","laplace","pmc","vb"))
          stop("Fit2 is not an object of class demonoid, iterquad, laplace, pmc, or vb.")
     if(identical(class(Fit1), "demonoid")) post1 <- Fit1$Posterior2
     if(identical(class(Fit2), "demonoid")) post2 <- Fit2$Posterior2
     if(identical(class(Fit1), "iterquad")) post1 <- Fit1$Posterior
     if(identical(class(Fit2), "iterquad")) post2 <- Fit2$Posterior
     if(identical(class(Fit1), "laplace")) post1 <- Fit1$Posterior
     if(identical(class(Fit2), "laplace")) post2 <- Fit2$Posterior
     if(identical(class(Fit1), "pmc")) post1 <- Fit1$Posterior2
     if(identical(class(Fit2), "pmc")) post2 <- Fit2$Posterior2
     if(identical(class(Fit1), "vb")) post1 <- Fit1$Posterior
     if(identical(class(Fit2), "vb")) post2 <- Fit2$Posterior
     if(class(Pred1) %!in% c("demonoid.ppc","iterquad.ppc","laplace.ppc","pmc.ppc","vb.ppc"))
          stop("Fit1 is not an object of class demonoid.ppc, iterquad.ppc, laplace.ppc, pmc.ppc, or vb.ppc.")
     if(class(Pred2) %!in% c("demonoid.ppc","iterquad.ppc","laplace.ppc","pmc.ppc","vb.ppc"))
          stop("Fit2 is not an object of class demonoid.ppc, iterquad.ppc, laplace.ppc, pmc.ppc, or vb.ppc.")
     if(nrow(post1) != nrow(post2))
          stop("The number of posterior samples differ between Fit1 and Fit2.")
     keep <- which(colnames(post1) %in% colnames(post2))
     post1 <- post1[,keep]
     keep <- which(colnames(post2) %in% colnames(post1))
     post2 <- post2[,keep]
     if(!identical(colnames(post1), colnames(post2)))
          stop("Posterior names differ between Fit1 and Fit2.")
     yhat1 <- Pred1$yhat
     yhat2 <- Pred2$yhat
     if(!identical(dim(yhat1), dim(yhat2)))
          stop("Dimensions of yhat differ between Pred1 and Pred2.")
     ### Sensitivity Analysis
     Posterior <- matrix(NA, ncol(post1), 2)
     rownames(Posterior) <- colnames(post1)
     colnames(Posterior) <- c("p(Fit1 > Fit2)", "var(Fit1) / var(Fit2)")
     Posterior[,1] <- colMeans(post1 > post2)
     Posterior[,2] <- .colVars(post1) / .colVars(post2)
     Post.Pred.Dist <- matrix(NA, nrow(yhat1), 2)
     rownames(Post.Pred.Dist) <- rownames(yhat1)
     colnames(Post.Pred.Dist) <- c("p(Pred1 > Pred2)",
          "var(Pred1) / var(Pred2)")
     Post.Pred.Dist[,1] <- rowMeans(yhat1 > yhat2)
     Post.Pred.Dist[,2] <- .rowVars(yhat1) / .rowVars(yhat2)
     ### Output
     out <- list(Posterior=Posterior, Post.Pred.Dist=Post.Pred.Dist)
     class(out) <- "sensitivity"
     return(out)
     }

#End
