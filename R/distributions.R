###########################################################################
# Asymmetric Laplace Distribution                                         #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dalaplace <- function(x, location=0, scale=1, kappa=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(x), length(location), length(scale), length(kappa))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     logconst <- 0.5 * log(2) - log(scale) + log(kappa) - log1p(kappa^2)
     temp <- which(x < location); kappa[temp] <- 1/kappa[temp]
     exponent <- -(sqrt(2) / scale) * abs(x - location) * kappa
     dens <- logconst + exponent
     if(log == FALSE) dens <- exp(logconst + exponent)
     return(dens)
     }
palaplace <- function(q, location=0, scale=1, kappa=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if((kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(q), length(location), length(scale), length(kappa))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- k2 <- rep(kappa, len=NN)
     temp <- which(q < location); k2[temp] <- 1/kappa[temp]
     exponent <- -(sqrt(2) / scale) * abs(q - location) * k2
     temp <- exp(exponent) / (1 + kappa^2)
     p <- 1 - temp
     index1 <- (q < location)
     p[index1] <- (kappa[index1])^2 * temp[index1]
     return(p)
     }
qalaplace <- function(p, location=0, scale=1, kappa=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(p), length(location), length(scale), length(kappa))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     q <- p
     temp <- kappa^2 / (1 + kappa^2)
     index1 <- (p <= temp)
     exponent <- p[index1] / temp[index1]
     q[index1] <- location[index1] + (scale[index1] * kappa[index1]) *
          log(exponent) / sqrt(2)
     q[!index1] <- location[!index1] - (scale[!index1] / kappa[!index1]) *
          (log1p((kappa[!index1])^2) + log1p(-p[!index1])) / sqrt(2)
     q[p == 0] = -Inf
     q[p == 1] = Inf
     return(q)
     }
ralaplace <- function(n, location=0, scale=1, kappa=1)
     {
     location <- rep(location, len=n)
     scale <- rep(scale, len=n)
     kappa <- rep(kappa, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     x <- location + scale *
          log(runif(n)^kappa / runif(n)^(1/kappa)) / sqrt(2)
     return(x)
     }

###########################################################################
# Asymmetric Log-Laplace Distribution                                     #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dallaplace <- function(x, location=0, scale=1, kappa=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(x), length(location), length(scale), length(kappa))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     Alpha <- sqrt(2) * kappa / scale
     Beta  <- sqrt(2) / (scale * kappa)
     Delta <- exp(location)
     exponent <- -(Alpha + 1)
     temp <- which(x >= Delta); exponent[temp] <- Beta[temp] - 1
     exponent <- exponent * (log(x) - location)
     dens <- -location + log(Alpha) + log(Beta) -
          log(Alpha + Beta) + exponent
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
pallaplace <- function(q, location=0, scale=1, kappa=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(q), length(location), length(scale), length(kappa))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     Alpha <- sqrt(2) * kappa / scale
     Beta  <- sqrt(2) / (scale * kappa)
     Delta <- exp(location)
     temp <- Alpha + Beta
     p <- (Alpha / temp) * (q / Delta)^(Beta)
     p[q <= 0] <- 0
     index1 <- (q >= Delta)
     p[index1] <- (1 - (Beta/temp) * (Delta/q)^(Alpha))[index1]
     return(p)
     }
qallaplace <- function(p, location=0, scale=1, kappa=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(p), length(location), length(scale), length(kappa))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     Alpha <- sqrt(2) * kappa / scale
     Beta  <- sqrt(2) / (scale * kappa)
     Delta <- exp(location)
     temp <- Alpha + Beta
     q <- Delta * (p * temp / Alpha)^(1/Beta)
     index1 <- (p > Alpha / temp)
     q[index1] <- (Delta * ((1-p) * temp / Beta)^(-1/Alpha))[index1]
     q[p == 0] <- 0
     q[p == 1] <- Inf
     return(q)
     }
rallaplace <- function(n, location=0, scale=1, kappa=1)
     {
     location <- rep(location, len=n); scale <- rep(scale, len=n)
     kappa <- rep(kappa, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     x <- exp(location) *
          (runif(n)^kappa / runif(n)^(1/kappa))^(scale / sqrt(2))
     return(x)
     }

###########################################################################
# Asymmetric Multivariate Laplace Distribution                            #
###########################################################################

daml <- function (x, mu, Sigma, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     x.Omega.mu <- rowSums(x %*% Omega * mu)
     x.Omega.x <- rowSums(x %*% Omega * x)
     mu.Omega.mu <- rowSums(mu %*% Omega * mu)
     dens <- as.vector(log(2) + x.Omega.mu -
          (log(2*pi)*(k/2) + logdet(Sigma)*0.5) +
          (log(x.Omega.x) - (log(2 + mu.Omega.mu)))*((2-k)/4) +
          log(besselK(sqrt((2 + mu.Omega.mu) * x.Omega.x), (2-k)/2)))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
raml <- function(n, mu, Sigma)
     {
     mu <- rbind(mu)
     if(missing(Sigma)) Sigma <- diag(ncol(mu))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     e <- matrix(rexp(n, 1), n, k)
     z <- rmvn(n, rep(0, k), Sigma)
     x <- mu*e + sqrt(e)*z
     return(x)
     }

###########################################################################
# Bernoulli Distribution                                                  #
#                                                                         #
# These functions are similar to those in the Rlab package.               #
###########################################################################

dbern <- function(x, prob, log=FALSE)
     {return(dbinom(x, 1, prob, log))}
pbern <- function(q, prob, lower.tail=TRUE, log.p=FALSE)
     {return(pbinom(q, 1, prob, lower.tail, log.p))}
qbern <- function(p, prob, lower.tail=TRUE, log.p=FALSE)
     {return(qbinom(p, 1, prob, lower.tail, log.p))}
rbern <- function(n, prob)
     {return(rbinom(n, 1, prob))}

###########################################################################
# Categorical Distribution                                                #
###########################################################################

dcat <- function (x, p, log=FALSE)
     {
     if(!is.matrix(p)) {
          if(is.vector(x)) p <- matrix(p, length(x), length(p), byrow=TRUE)
          else p <- matrix(p, nrow(x), length(p), byrow=TRUE)}
     if(is.vector(x)) {
          if(length(x) == 1) {
               temp <- rbind(rep(0, ncol(p)))
               temp[1,x] <- 1
               x <- temp
               }
          else x <- as.indicator.matrix(x)}
     if(!identical(nrow(x), nrow(p)))
          stop("The number of rows of x and p differ.")
     if(!identical(ncol(x), ncol(p))) {
          x.temp <- matrix(0, nrow(p), ncol(p))
          x.temp[, as.numeric(colnames(x))] <- x
          x <- x.temp}
     dens <- x * log(p)
     if(log == FALSE) dens <- x * p
     dens <- as.vector(rowSums(dens))
     return(dens)
     }
qcat <- function(pr, p, lower.tail=TRUE, log.pr=FALSE)
     {
     if(!is.vector(pr)) pr <- as.vector(pr)
     if(!is.vector(p)) p <- as.vector(p)
     if(log.pr == FALSE) {
          if(any(pr < 0) | any(pr > 1))
               stop("pr must be in [0,1].")}
     else if(any(!is.finite(pr)) | any(pr > 0))
          stop("pr, as a log, must be in (-Inf,0].")
     if(sum(p) != 1) stop("sum(p) must be 1.")
     if(lower.tail == FALSE) pr <- 1 - pr
     breaks <- c(0, cumsum(p))
     if(log.pr == TRUE) breaks <- log(breaks)
     breaks <- matrix(breaks, length(pr), length(breaks), byrow=TRUE)
     x <- rowSums(pr > breaks)
     return(x)
     }
rcat <- function(n, p)
     {
     if(is.vector(p)) {
          x <- as.vector(which(rmultinom(n, size=1, prob=p) ==
               1, arr.ind=TRUE)[, "row"])
          }
     else {
          d <- dim(p)
          n <- d[1]
          k <- d[2]
          lev <- dimnames(p)[[2]]
          if(!length(lev)) lev <- 1:k
          z <- colSums(p)
          U <- apply(p, 1, cumsum)
          U[,k] <- 1
          un <- rep(runif(n), rep(k,n))
          x <- lev[1 + colSums(un > U)]}
     return(x)
     }

###########################################################################
# Continuous Relaxation of a Markov Random Field Distribution             #
###########################################################################

dcrmrf <- function(x, alpha, Omega, log=FALSE)
     {
     alpha <- as.vector(alpha)
     if(missing(Omega)) Omega <- diag(length(alpha))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     Omega <- as.symmetric.matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     dens <- as.vector(-0.5*t(x) %*% as.inverse(Omega) %*% x) +
          log(prod(1 + exp(x + alpha)))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rcrmrf <- function(n=1, alpha, Omega)
     {
     alpha <- as.vector(alpha)
     J <- length(alpha)
     if(missing(Omega)) Omega <- diag(J)
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     dens <- rep(0,J)
     x <- rep(0,n)
     for (i in 1:n) {
          for (j in 1:J) {
               z <- as.vector(rmvn(1,
                    as.vector(Omega %*% diag(J)[j,]), Omega))
               dens[j] <- dcrmrf(z, alpha, Omega, log=FALSE)}
          x[i] <- sample(1:J, size=1, replace=TRUE, prob=dens)}
     return(x)
     }

###########################################################################
# Dirichlet Distribution                                                  #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

ddirichlet <- function(x, alpha, log=FALSE)
     {
     if(missing(x)) stop("x is a required argument.")
     if(missing(alpha)) stop("alpha is a required argument.")
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(alpha))
          alpha <- matrix(alpha, nrow(x), length(alpha), byrow=TRUE)
     if(any(rowSums(x) != 1)) x / rowSums(x)
     if(any(x < 0)) stop("x must be non-negative.")
     if(any(alpha <= 0)) stop("alpha must be positive.")
     dens <- as.vector(lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) +
          rowSums((alpha-1)*log(x)))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rdirichlet <- function (n, alpha)
     {
     alpha <- rbind(alpha)
     alpha.dim <- dim(alpha)
     if(n > alpha.dim[1])
          alpha <- matrix(alpha, n, alpha.dim[2], byrow=TRUE)
     x <- matrix(rgamma(alpha.dim[2]*n, alpha), ncol=alpha.dim[2])
     sm <- x %*% rep(1, alpha.dim[2])
     return(x/as.vector(sm))
     }

###########################################################################
# Generalized Pareto Distribution                                         #
###########################################################################

dgpd <- function(x, mu, sigma, xi, log=FALSE)
     {
     x <- as.vector(x)
     mu <- as.vector(mu)
     sigma <- as.vector(sigma)
     xi <- as.vector(xi)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     NN <- max(length(x), length(mu), length(sigma), length(xi))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); xi <- rep(xi, len=NN)
     xi.ge.0 <- which(xi >= 0)
     xi.lt.0 <- which(xi < 0)
     if(any(mu > x) |
     any(mu[xi.lt.0] < x[xi.lt.0] + sigma[xi.lt.0] / xi[xi.lt.0]))
     stop("x is outside of support.")
     z <- (x - mu) / sigma
     xi0 <- which(xi == 0)
     xi1 <- which(xi != 0)
     dens <- rep(NA, NN)
     dens[xi0] <- -z[xi0] - log(sigma[xi0])
     dens[xi1] <- log(1/sigma[xi1]) + log(1 + xi[xi1] *
     z[xi1]) * (-1/xi[xi1] - 1)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rgpd <- function(n, mu, sigma, xi)
     {
     mu <- as.vector(mu)
     sigma <- as.vector(sigma)
     xi <- as.vector(xi)
     if(any(sigma <= 0)) stop("The sigma parameter must be non-negative.")
     mu <- rep(mu, len=n); sigma <- rep(sigma, len=n)
     xi <- rep(xi, len=n)
     x <- rep(NA,n)
     u <- runif(n, min=0.001)
     xi0 <- which(xi == 0)
     xi1 <- which(xi != 0)
     x[xi0] <- mu[xi0] - sigma[xi0] * log(u[xi0])
     x[xi1] <- mu[xi1] + sigma[xi1] *
     (u[xi1]^(-xi[xi1]) - 1) / xi[xi1]
     return(x)
     }

###########################################################################
# Generalized Poisson                                                     #
###########################################################################

dgpois <- function(x, lambda=0, omega=0, log=FALSE)
     {
     x <- as.vector(x); lambda <- as.vector(lambda)
     omega <- as.vector(omega)
     lambda.star <- (1 - omega)*lambda + omega*x
     dens <- log(1 - omega) + log(lambda) + (x - 1)*log(lambda.star) -
          lgamma(x + 1) - lambda.star
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }

###########################################################################
# Half-Cauchy Distribution                                                #
###########################################################################

dhalfcauchy <- function(x, scale=25, log=FALSE)
     {
     x <- as.vector(x); scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(scale))
     x <- rep(x, len=NN); scale <- rep(scale, len=NN)
     dens <- log(2*scale) - log(pi*{x*x + scale*scale})
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
phalfcauchy <- function(q, scale=25)
     {
     q <- as.vector(q); scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(scale))
     q <- rep(q, len=NN); scale <- rep(scale, len=NN)
     z <- {2/pi}*atan(q/scale)
     return(z)
     }
qhalfcauchy <- function(p, scale=25)
     {
     p <- as.vector(p); scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(scale))
     p <- rep(p, len=NN); scale <- rep(scale, len=NN)
     q <- scale*tan({pi*p}/2)
     return(q)
     }
rhalfcauchy <- function(n, scale=25)
     {
     scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     p <- runif(n, 0, 1)
     x <- scale*tan({pi*p}/2)
     return(x)
     }

###########################################################################
# Half-Normal Distribution                                                #
#                                                                         #
# This half-normal distribution has mean=0 and is similar to the halfnorm #
# functions in package fdrtool.                                           #
###########################################################################

dhalfnorm <- function(x, scale=sqrt(pi/2), log=FALSE)
     {
     dens <- log(2) + dnorm(x, mean=0, sd=sqrt(pi/2) / scale, log=TRUE)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
phalfnorm <- function(q, scale=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(scale))
     q <- rep(q, len=NN); scale <- rep(scale, len=NN)
     p <- 2*pnorm(q, mean=0, sd=sqrt(pi/2) / scale) - 1
     if(lower.tail == FALSE) p <- 1-p
     if(log.p == TRUE) p <- log.p(p)
     return(p)
     }
qhalfnorm <- function(p, scale=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(scale))
     p <- rep(p, len=NN); scale <- rep(scale, len=NN)
     if(log.p == TRUE) p <- exp(p)
     if(lower.tail == FALSE) p <- 1-p
     q <- qnorm((p+1)/2, mean=0, sd=sqrt(pi/2) / scale)
     return(q)
     }
rhalfnorm <- function(n, scale=sqrt(pi/2))
     {
     scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     x <- abs(rnorm(n, mean=0, sd=sqrt(pi/2) / scale))
     return(x)
     }

###########################################################################
# Half-t Distribution                                                     #
###########################################################################

dhalft <- function(x, scale=25, nu=1, log=FALSE)
     {
     x <- as.vector(x); scale <- as.vector(scale); nu <- as.vector(nu)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(scale), length(nu))
     x <- rep(x, len=NN); scale <- rep(scale, len=NN)
     nu <- rep(nu, len=NN)
     dens <- log(2) -log(scale) + lgamma((nu + 1)/2) - lgamma(nu/2) -
       .5*log(pi*nu) - (nu + 1)/2 * log(1 + (1/nu) * (x/scale) * (x/scale))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
phalft <- function(q, scale=25, nu=1)
     {
     q <- as.vector(q); scale <- as.vector(scale); nu <- as.vector(nu)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(scale), length(nu))
     q <- rep(q, len=NN); scale <- rep(scale, len=NN)
     p <- ptrunc(q, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(p)
     }
qhalft <- function(p, scale=25, nu=1)
     {
     p <- as.vector(p); scale <- as.vector(scale); nu <- as.vector(nu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(scale), length(nu))
     p <- rep(p, len=NN); scale <- rep(scale, len=NN)
     q <- rtrunc(p, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(q)
     }
rhalft <- function(n, scale=25, nu=1)
     {
     scale <- rep(scale, len=n); nu <- rep(nu, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     x <- rtrunc(n, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(x)
     }

###########################################################################
# Horseshoe Distribution                                                  #
###########################################################################

dhs <- function(x, lambda, tau, log=FALSE)
     {
     dens <- dnorm(x, 0, lambda*tau, log=log)
     return(dens)
     }
rhs <- function(n, lambda, tau)
     {
     x <- rnorm(n, 0, lambda*tau)
     return(x)
     }

###########################################################################
# Huang-Wand Distribution                                                 #
###########################################################################

dhuangwand <- function(x, nu=2, a, A, log=FALSE)
     {
     if(!is.matrix(x)) x <- matrix(x)
     if(!is.positive.definite(x))
          stop("Matrix x is not positive-definite.")
     a <- as.vector(a)
     A <- as.vector(A)
     k <- nrow(x)
     if(!identical(length(a), length(A), k, ncol(x)))
          stop("Dimensions of x, a, and A do not agree.")
     dens <- sum(dinvgamma(a, 0.5, 1/A^2, log=TRUE)) +
          dinvwishart(x, nu + k - 1, 2*nu*diag(1/a), log=TRUE)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rhuangwand <- function(nu=2, a, A)
     {
     if(missing(A)) k <- length(a)
     else {
          k <- length(A)
          a <- rinvgamma(k, 0.5, 1/A^2)}
     x <- rinvwishart(nu + k - 1, 2*nu*diag(1/a))
     return(x)
     }

###########################################################################
# Huang-Wand Distribution (Cholesky Parameterization)                     #
###########################################################################

dhuangwandc <- function(x, nu=2, a, A, log=FALSE)
     {
     if(!is.matrix(x)) x <- matrix(x)
     a <- as.vector(a)
     A <- as.vector(A)
     k <- nrow(x)
     if(!identical(length(a), length(A), k, ncol(x)))
          stop("Dimensions of x, a, and A do not agree.")
     dens <- sum(dinvgamma(a, 0.5, 1/A^2, log=TRUE)) +
          dinvwishartc(x, nu + k - 1, 2*nu*diag(1/a), log=TRUE)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rhuangwandc <- function(nu=2, a, A)
     {
     if(missing(A)) k <- length(a)
     else {
          k <- length(A)
          a <- rinvgamma(k, 0.5, 1/A^2)}
     x <- rinvwishartc(nu + k - 1, 2*nu*diag(1/a))
     return(x)
     }

###########################################################################
# Inverse Beta Distribution                                               #
###########################################################################

dinvbeta <- function(x, a, b, log=FALSE)
     {
     const <- lgamma(a + b) - lgamma(a) - lgamma(b)
     dens <- const + (a - 1) * log(x) - (a + b) * log(1 + x)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rinvbeta <- function(n, a, b)
     {
     x <- rbeta(n, a, b)
     x <- x / (1-x)
     return(x)
     }

###########################################################################
# Inverse Chi-Squared Distribution                                        #
#                                                                         #
# These functions are similar to those in the GeoR package.               #
###########################################################################

dinvchisq <- function(x, df, scale=1/df, log=FALSE)
     {
     x <- as.vector(x); df <- as.vector(df); scale <- as.vector(scale)
     if(any(x <= 0)) stop("x must be positive.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(df), length(scale))
     x <- rep(x, len=NN); df <- rep(df, len=NN);
     scale <- rep(scale, len=NN)
     nu <- df / 2
     dens <- nu*log(nu) - lgamma(nu) + nu*log(scale) -
          (nu+1)*log(x) - (nu*scale/x)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }

rinvchisq <- function(n, df, scale=1/df)
     {
     df <- rep(df, len=n); scale <- rep(scale, len=n)
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     z <- rchisq(n, df=df)
     z[which(z == 0)] <- 1e-100
     x <- (df*scale) / z
     return(x)
     }

###########################################################################
# Inverse Gamma Distribution                                              #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

dinvgamma <- function(x, shape=1, scale=1, log=FALSE)
     {
     x <- as.vector(x); shape <- as.vector(shape)
     scale <- as.vector(scale)
     if(any(shape <= 0) | any(scale <=0))
          stop("The shape and scale parameters must be positive.")
     NN <- max(length(x), length(shape), length(scale))
     x <- rep(x, len=NN); shape <- rep(shape, len=NN)
     scale <- rep(scale, len=NN)
     alpha <- shape; beta <- scale
     dens <- alpha * log(beta) - lgamma(alpha) -
          {alpha + 1} * log(x) - {beta/x}
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rinvgamma <- function(n, shape=1, scale=1)
     {
     x <- rgamma(n=n, shape=shape, rate=scale)
     x[which(x < 1e-300)] <- 1e-300
     return(1 / x)
     }

###########################################################################
# Inverse Gaussian Distribution                                           #
###########################################################################

dinvgaussian <- function(x, mu, lambda, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); lambda <- as.vector(lambda)
     if(any(x <= 0)) stop("x must be positive.")
     if(any(mu <= 0)) stop("The mu parameter must be positive.")
     if(any(lambda <= 0)) stop("The lambda parameter must be positive.")
     NN <- max(length(x), length(mu), length(lambda))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     lambda <- rep(lambda, len=NN)
     dens <- log(lambda^0.5/(2 * pi * x^3)^0.5) - ((lambda * (x - mu)^2)/(2 * mu^2 * x))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rinvgaussian <- function(n, mu, lambda)
     {
     mu <- rep(mu, len=n); lambda <- rep(lambda, len=n)
     if(any(mu <= 0)) stop("The mu parameter must be positive.")
     if(any(lambda <= 0)) stop("The lambda parameter must be positive.")
     nu <- rnorm(n)
     y <- nu^2
     x <- mu + ((mu^2*y)/(2*lambda)) - (mu/(2*lambda)) *
          sqrt(4*mu*lambda*y + mu^2*y^2)
     z <- runif(n)
     temp <- which(z > {mu / (mu + x)})
     x[temp] <- mu[temp] * mu[temp] / x[temp]
     return(x)
     }

###########################################################################
# Inverse Matrix Gamma Distribution                                       #
###########################################################################

dinvmatrixgamma <- function(X, alpha, beta, Psi, log=FALSE)
     {
     if(!is.matrix(X)) stop("X is not a matrix.")
     if(alpha <= 2) stop("The alpha parameter must be greater than 2.")
     if(beta <= 0) stop("The beta parameter must be positive.")
     if(missing(Psi)) Psi <- diag(ncol(X))
     if(!is.positive.definite(Psi))
          stop("Matrix Psi is not positive-definite.")
     k <- nrow(Psi)
     gamsum <- 0
     for (i in 1:k) gamsum <- gamsum + lgamma(alpha - 0.5*(i-1))
     gamsum <- gamsum + log(pi)*(k*(k-1)/4)
     Omega <- as.inverse(Psi)
     dens <- logdet(Psi)*alpha - (log(beta)*(k*alpha) + gamsum) +
          logdet(X)*(-alpha-(k+1)/2) +
          tr(-(1/beta)*(Psi %*% as.inverse(X)))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }

rinvmatrixgamma <- function(alpha, beta, Psi){
  rinvwishart(2*alpha, 2/beta*Psi)
}

###########################################################################
# Inverse Wishart Distribution                                            #
###########################################################################

dinvwishart <- function(Sigma, nu, S, log=FALSE)
     {
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(!identical(dim(S), dim(Sigma)))
          stop("The dimensions of Sigma and S differ.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     k <- nrow(Sigma)
     gamsum <- 0
     for (i in 1:k) {gamsum <- gamsum + lgamma((nu + 1 - i)/2)}
     dens <- -(nu*k/2)*log(2) - ((k*(k - 1))/4)*log(pi) - gamsum +
          (nu/2)*log(det(S)) - ((nu + k + 1)/2)*logdet(Sigma) -
          0.5*tr(S %*% as.inverse(Sigma))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rinvwishart <- function(nu, S)
     {return(chol2inv(rwishartc(nu, chol2inv(chol(S)))))}

###########################################################################
# Inverse Wishart Distribution (Cholesky Parameterization)                #
###########################################################################

dinvwishartc <- function(U, nu, S, log=FALSE)
     {
     if(missing(U)) stop("Upper triangular U is required.")
     Sigma <- t(U) %*% U
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(!identical(dim(S), dim(Sigma)))
          stop("The dimensions of Sigma and S differ.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     k <- nrow(Sigma)
     gamsum <- 0
     for (i in 1:k) {gamsum <- gamsum + lgamma((nu + 1 - i)/2)}
     dens <- -(nu*k/2)*log(2) - ((k*(k - 1))/4)*log(pi) - gamsum +
          (nu/2)*log(det(S)) - ((nu + k + 1)/2)*logdet(Sigma) -
          0.5*tr(S %*% as.inverse(Sigma))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rinvwishartc <- function(nu, S)
     {return(chol(as.inverse(rwishart(nu, as.inverse(S)))))}

###########################################################################
# Laplace Distribution                                                    #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dlaplace <- function(x, location=0, scale=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(location), length(scale))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     dens <- (-abs(x - location) / scale) - log(2 * scale)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
plaplace <- function(q, location=0, scale=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     z <- {q - location} / scale
     NN <- max(length(q), length(location), length(scale))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     p <- q
     temp <- which(q < location); p[temp] <- 0.5 * exp(z[temp])
     temp <- which(q >= location); p[temp] <- 1 - 0.5 * exp(-z[temp])
     return(p)
     }
qlaplace <- function(p, location=0, scale=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(location), length(scale))
     p <- p2 <- rep(p, len=NN); location <- rep(location, len=NN)
     temp <- which(p > 0.5); p2[temp] <- 1 - p[temp]
     q <- location - sign(p - 0.5) * scale * log(2 * p2)
     return(q)
     }
rlaplace <- function(n, location=0, scale=1)
     {
     location <- rep(location, len=n); scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     r <- r2 <- runif(n)
     temp <- which(r > 0.5); r2[temp] <- 1 - r[temp]
     x <- location - sign(r - 0.5) * scale * log(2 * r2)
     return(x)
     }

###########################################################################
# Laplace Distribution (Precision Parameterization)                       #
###########################################################################

dlaplacep <- function(x, mu=0, tau=1, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(x), length(mu), length(tau))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     dens <- log(tau/2) + (-tau*abs(x-mu))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
plaplacep <- function(q, mu=0, tau=1)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(q), length(mu), length(tau))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     p <- plaplace(q, mu, scale=1/tau)
     return(p)
     }
qlaplacep <- function(p, mu=0, tau=1)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(p), length(mu), length(tau))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     q <- qlaplace(p, mu, scale=1/tau)
     return(q)
     }
rlaplacep <- function(n, mu=0, tau=1)
     {
     mu <- rep(mu, len=n); tau <- rep(tau, len=n)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     x <- rlaplace(n, mu, scale=1/tau)
     return(x)
     }

###########################################################################
# Laplace Distribution Mixture                                            #
###########################################################################

dlaplacem <- function(x, p, location, scale, log=FALSE)
     {
     if(missing(x)) stop("x is a required argument.")
     x <- as.vector(x)
     n <- length(x)
     if(missing(p)) stop("p is a required argument.")
     p <- as.vector(p)
     if(any(p <= 0) | any(p > 1)) stop("p must be in (0,1].")
     if(sum(p) != 1) stop("p must sum to 1 for all components.")
     m <- length(p)
     p <- matrix(p, n, m, byrow=TRUE)
     if(missing(location)) stop("location is a required argument.")
     location <- as.vector(location)
     if(!identical(m, length(location)))
          stop("p and location differ in length.")
     location <- matrix(location, n, m, byrow=TRUE)
     if(missing(scale)) stop("scale is a required argument.")
     scale <- as.vector(scale)
     if(!identical(m, length(scale)))
          stop("p and scale differ in length.")
     scale <- matrix(scale, n, m, byrow=TRUE)
     dens <- matrix(dlaplace(x, location, scale, log=TRUE), n, m)
     dens <- dens + log(p)
     if(log == TRUE) dens <- apply(dens, 1, logadd)
     else dens <- rowSums(exp(dens))
     return(dens)
     }
plaplacem <- function(q, p, location, scale)
     {
     n <- length(q)
     m <- length(p)
     q <- matrix(q, n, m)
     p <- matrix(p, n, m, byrow=TRUE)
     location <- matrix(location, n, m, byrow=TRUE)
     scale <- matrix(scale, n, m, byrow=TRUE)
     cdf <- matrix(plaplace(q, location, scale), n, m)
     cdf <- rowSums(cdf * p)
     return(cdf)
     }
rlaplacem <- function(n, p, location, scale)
     {
     if(missing(p)) stop("p is a required argument.")
     p <- as.vector(p)
     if(any(p <= 0) | any(p > 1)) stop("p must be in (0,1].")
     if(sum(p) != 1) stop("p must sum to 1 for all components.")
     m <- length(p)
     p <- matrix(p, n, m, byrow=TRUE)
     if(missing(location)) stop("location is a required argument.")
     location <- as.vector(location)
     if(!identical(m, length(location)))
          stop("p and location differ in length.")
     if(missing(scale)) stop("scale is a required argument.")
     scale <- as.vector(scale)
     if(!identical(m, length(scale)))
          stop("p and scale differ in length.")
     if(any(scale <= 0)) stop("scale must be positive.")
     z <- rcat(n, p)
     x <- rlaplace(n, location=location[z], scale=scale[z])
     return(x)
     }

###########################################################################
# LASSO Distribution                                                      #
###########################################################################

dlasso <- function(x, sigma, tau, lambda, a=1, b=1, log=FALSE)
     {
     if(any(c(sigma, tau, lambda, a, b) <= 0))
          stop("Scale parameters must be positive.")
     dens <- sum(dmvn(x, 0, sigma^2*diag(tau^2), log=TRUE)) +
          log(1/sigma^2) + sum(dexp(tau^2, lambda^2/2, log=TRUE)) +
          dgamma(lambda^2, a, b, log=TRUE)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rlasso <- function(n, sigma, tau, lambda, a=1, b=1)
     {
     a <- as.vector(a)[1]
     b <- as.vector(b)[1]
     if(missing(sigma)) sigma <- runif(1, 1e-100, 1000)
     if(missing(lambda)) lambda <- sqrt(rgamma(1, a, b))
     if(missing(tau)) stop("The tau parameter is a required argument.")
     sigma <- as.vector(sigma)[1]
     lambda <- as.vector(lambda)[1]
     x <- rnorm(n, 0, sigma^2*diag(tau^2))
     return(x)
     }

###########################################################################
# Log-Laplace Distribution                                                #
###########################################################################

dllaplace <- function(x, location=0, scale=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(location), length(scale))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     Alpha <- sqrt(2) * scale
     Beta  <- sqrt(2) / scale
     Delta <- exp(location)
     exponent <- -(Alpha + 1)
     temp <- which(x < Delta); exponent[temp] <- Beta[temp] - 1
     exponent <- exponent * (log(x) - location)
     dens <- -location + log(Alpha) + log(Beta) -
          log(Alpha + Beta) + exponent
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
pllaplace <- function(q, location=0, scale=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(location), length(scale))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     Alpha <- sqrt(2) * scale
     Beta  <- sqrt(2) / scale
     Delta <- exp(location)
     temp <- Alpha + Beta
     p <- (Alpha / temp) * (q / Delta)^(Beta)
     p[q <= 0] <- 0
     index1 <- (q >= Delta)
     p[index1] <- (1 - (Beta/temp) * (Delta/q)^(Alpha))[index1]
     return(p)
     }
qllaplace <- function(p, location=0, scale=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(location), length(scale))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     Alpha <- sqrt(2) * scale
     Beta  <- sqrt(2) / scale
     Delta <- exp(location)
     temp <- Alpha + Beta
     q <- Delta * (p * temp / Alpha)^(1/Beta)
     index1 <- (p > Alpha / temp)
     q[index1] <- (Delta * ((1-p) * temp / Beta)^(-1/Alpha))[index1]
     q[p == 0] <- 0
     q[p == 1] <- Inf
     return(q)
     }
rllaplace <- function(n, location=0, scale=1)
     {
     location <- rep(location, len=n); scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     x <- exp(location) * (runif(n) / runif(n))^(scale / sqrt(2))
     return(x)
     }

###########################################################################
# Log-Normal Distribution (Precision Parameterization)                    #
###########################################################################

dlnormp <- function(x, mu, tau = NULL, var = NULL, log = FALSE)
{
	if (!is.null(tau) & !is.null(var))
		stop("use either tau or var, but not both")

	x <- as.vector(x); mu <- as.vector(mu);
	if (!is.null(tau))
	{
		tau <- as.vector(tau)
		if(any(tau <= 0))
			stop("The tau parameter must be positive.")

		NN <- max(length(x), length(mu), length(tau))
		x <- rep(x, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
		dens <- 0.5*(log(tau) - log(2*pi)) - log(x) - tau/2*(log(x) - mu)^2
	}

	if (!is.null(var))
	{
		var <- as.vector(var)
		if(any(var <= 0))
			stop("The var parameter must be positive.")

		NN <- max(length(x), length(mu), length(var))
		x <- rep(x, len=NN); mu <- rep(mu, len=NN); var <- rep(var, len=NN)
		dens <- -0.5*log(2*pi*var) - log(x) - (log(x) - mu)^2/(2*var)
	}

	if(log == FALSE)
		dens <- exp(dens)

	return(dens)
}
plnormp <- function(q, mu, tau, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(q), length(mu), length(tau))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     p <- pnorm(q, mu, sqrt(1/tau), lower.tail, log.p)
     return(p)
     }
qlnormp <- function(p, mu, tau, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(p), length(mu), length(tau))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     q <- qnorm(p, mu, sqrt(1/tau), lower.tail, log.p)
     return(q)
}

rlnormp <- function(n, mu, tau = NULL, var = NULL)
{
	if (!is.null(tau) & !is.null(var))
		stop("use either tau or var, but not both")

	mu <- rep(mu, len=n);
	if (!is.null(tau))
	{
		if(any(tau <= 0))
			stop("The tau parameter must be positive.")

		tau <- rep(tau, len=n)
		x <- rlnorm(n, log(mu^2 / sqrt(1/tau^2 + mu^2)), sqrt(log(1 + 1/(mu^2*tau^2))))
	}
	if (!is.null(var))
	{
		var = rep(var, len=n)
		if(any(var <= 0))
			stop("The var parameter must be positive.")

		var <- rep(var, len=n)
		x <- rlnorm(n, log(mu^2 / sqrt(var^2 + mu^2)), sqrt(log(1 + (var^2/mu^2))))
	}
	return(x)
}

###########################################################################
# Matrix Gamma Distribution                                               #
###########################################################################

dmatrixgamma <- function(X, alpha, beta, Sigma, log=FALSE)
     {
     if(!is.matrix(X)) stop("X is not a matrix.")
     if(alpha <= 2) stop("The alpha parameter must be greater than 2.")
     if(beta <= 0) stop("The beta parameter must be positive.")
     if(missing(Sigma)) Sigma <- diag(ncol(X))
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- nrow(Sigma)
     gamsum <- 0
     for (i in 1:k) gamsum <- gamsum + lgamma(alpha - 0.5*(i-1))
     gamsum <- gamsum + log(pi)*(k*(k-1)/4)
     Omega <- as.inverse(Sigma)
     dens <- alpha*logdet(Omega) + logdet(X)*(alpha - 0.5*(k+1)) -
          (log(beta)*(k*alpha) + gamsum) + (-1/beta)*tr(Omega %*% X)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }

rmatrixgamma <- function(alpha, beta, Sigma){
  rwishart(2*alpha, beta/2*Sigma)
}
###########################################################################
# Matrix Normal Distribution                                              #
###########################################################################

dmatrixnorm <- function(X, M, U, V, log=FALSE)
     {
     if(!is.matrix(X)) X <- rbind(X)
     if(!is.matrix(M)) M <- matrix(M, nrow(X), ncol(X), byrow=TRUE)
     if(missing(U)) U <- diag(nrow(X))
     if(!is.matrix(U)) U <- matrix(U)
     if(!is.positive.definite(U))
          stop("Matrix U is not positive-definite.")
     if(missing(V)) V <- diag(ncol(X))
     if(!is.matrix(V)) V <- matrix(V)
     if(!is.positive.definite(V))
          stop("Matrix V is not positive-definite.")
     n <- nrow(X)
     k <- ncol(X)
     ss <- X - M
     dens <- -0.5*tr(as.inverse(V) %*% t(ss) %*%
          as.inverse(U) %*% ss) - 
       (log(2*pi)*(n*k/2) - logdet(V) * (n/2) - logdet(U) * (k/2))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmatrixnorm <- function(M, U, V)
     {
     if(missing(M)) stop("Matrix M is missing.")
     if(!is.matrix(M)) M <- matrix(M)
     if(missing(U)) stop("Matrix U is missing.")
     if(!is.matrix(U)) U <- matrix(U)
     if(!is.positive.definite(U))
          stop("Matrix U is not positive-definite.")
     if(missing(V)) stop("Matrix V is missing.")
     if(!is.matrix(V)) V <- matrix(V)
     if(!is.positive.definite(V))
          stop("Matrix V is not positive-definite.")
     if(nrow(M) != nrow(U))
          stop("Dimensions of M and U are incorrect.")
     if(ncol(M) != ncol(V))
          stop("Dimensions of M and V are incorrect.")
     n <- nrow(U)
     k <- ncol(V)
     Z <- matrix(rnorm(n * k), n, k)
     X <- M + t(chol(U)) %*% Z %*% chol(V)
     return(X)
     }

###########################################################################
# Multivariate Cauchy Distribution                                        #
###########################################################################

dmvc <- function(x, mu, S, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(S)) S <- diag(ncol(x))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     k <- nrow(S)
     ss <- x - mu
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((1 + k)/2) - (lgamma(0.5) +
          (k/2)*log(pi) + 0.5*logdet(S) + ((1+k)/2)*log(1+z)))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvc <- function(n=1, mu=rep(0,k), S)
     {
     mu <- rbind(mu)
     if(missing(S)) S <- diag(ncol(mu))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     k <- ncol(S)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     x <- rchisq(n,1)
     x[which(x == 0)] <- 1e-100
     z <- rmvn(n, rep(0,k), S)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Multivariate Cauchy Distribution (Cholesky Parameterization)            #
###########################################################################

dmvcc <- function(x, mu, U, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- nrow(U)
     S <- t(U) %*% U
     ss <- x - mu
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma(k/2) - (lgamma(0.5) + log(1^(k/2)) +
          (k/2)*log(pi) + 0.5*logdet(S) + ((1+k)/2)*log(1+z)))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvcc <- function(n=1, mu=rep(0,k), U)
     {
     mu <- rbind(mu)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     S <- t(U) %*% U
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     x <- rchisq(n,1)
     x[which(x == 0)] <- 1e-100
     z <- rmvnc(n, rep(0,k), U)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Multivariate Cauchy Distribution (Precision Parameterization)           #
###########################################################################

dmvcp <- function(x, mu, Omega, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     k <- nrow(Omega)
     logdetOmega <- logdet(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((1+k)/2) - (lgamma(0.5) + log(1^(k/2)) +
          (k/2)*log(pi)) + 0.5*logdetOmega + (-(1+k)/2)*log(1 + z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvcp <- function(n=1, mu, Omega)
     {
     mu <- rbind(mu)
     if(missing(Omega)) Omega <- diag(ncol(mu))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     Sigma <- as.inverse(Omega)
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     x <- rchisq(n,1)
     x[which(x == 0)] <- 1e-100
     z <- rmvn(n, rep(0,k), Sigma)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Multivariate Cauchy Distribution (Precision-Cholesky Parameterization)  #
###########################################################################

dmvcpc <- function(x, mu, U, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- nrow(U)
     Omega <- t(U) %*% U
     logdetOmega <- logdet(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((1+k)/2) - (lgamma(0.5) + log(1^(k/2)) +
          (k/2)*log(pi)) + 0.5*logdetOmega + (-(1+k)/2)*log(1 + z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvcpc <- function(n=1, mu, U)
     {
     mu <- rbind(mu)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     x <- rchisq(n,1)
     x[which(x == 0)] <- 1e-100
     z <- rmvnc(n, rep(0,k), U)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Multivariate Laplace Distribution                                       #
###########################################################################

dmvl <- function(x, mu, Sigma, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma))
         stop("Matrix Sigma is not positive-definite.")
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     z[which(z == 0)] <- 1e-300
     dens <- as.vector(log(2) - log(2 * pi) * (k/2) - logdet(Sigma) * 0.5 + 
                      (log(pi) - log(2) - log(2 * z) * 0.5) * 0.5 - 
                      sqrt(2 * z) - log(z/2) * 0.5 * (k/2 - 1))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvl <- function(n, mu, Sigma)
     {
     mu <- rbind(mu)
     if(missing(Sigma)) Sigma <- diag(ncol(mu))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     e <- matrix(rexp(n, 1), n, k)
     z <- rmvn(n, rep(0, k), Sigma)
     x <- mu + sqrt(e)*z
     return(x)
     }

###########################################################################
# Multivariate Laplace Distribution (Cholesky Parameterization)           #
###########################################################################

dmvlc <- function(x, mu, U, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     Sigma <- t(U) %*% U
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     z[which(z == 0)] <- 1e-300
     dens <- as.vector(log(2) - log(2*pi)*(k/2) + logdet(Sigma)*0.5 +
          (log(pi) - log(2) + log(2*z)*0.5)*0.5 - log(2*z)*0.5 -
          log(z/2)*0.5*(k/2 - 1))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvlc <- function(n, mu, U)
     {
     mu <- rbind(mu)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     e <- matrix(rexp(n, 1), n, k)
     z <- rmvnc(n, rep(0, k), U)
     x <- mu + sqrt(e)*z
     return(x)
     }

###########################################################################
# Multivariate Normal Distribution                                        #
###########################################################################

dmvn <- function(x, mu, Sigma, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(-0.5 * (k * log(2 * pi) + logdet(Sigma)) - (0.5 * z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvn <- function(n=1, mu=rep(0,k), Sigma)
     {
     mu <- rbind(mu)
     if(missing(Sigma)) Sigma <- diag(ncol(mu))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     z <- matrix(rnorm(n*k),n,k) %*% chol(Sigma)
     x <- mu + z
     return(x)
     }

###########################################################################
# Multivariate Normal Distribution (Cholesky Parameterization)            #
###########################################################################

dmvnc <- function(x, mu, U, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     Sigma <- t(U) %*% U
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(-0.5 * (k * log(2 * pi) + logdet(Sigma)) - (0.5 * z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvnc <- function(n=1, mu=rep(0,k), U)
     {
     mu <- rbind(mu)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     z <- matrix(rnorm(n*k),n,k) %*% U
     x <- mu + z
     return(x)
     }

###########################################################################
# Multivariate Normal Distribution (Precision Parameterization)           #
###########################################################################

dmvnp <- function(x, mu, Omega, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     k <- nrow(Omega)
     logdetOmega <- logdet(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector((-k/2)*log(2*pi) + 0.5*logdetOmega - 0.5*z)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvnp <- function(n=1, mu=rep(0, k), Omega)
     {
     mu <- rbind(mu)
     if(missing(Omega)) Omega <- diag(ncol(mu))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     k <- ncol(Omega)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     z <- matrix(rnorm(n*k),n,k) %*% solve(t(chol(Omega)))
     x <- mu + z
     return(x)
     }

###########################################################################
# Multivariate Normal Distribution (Precision-Cholesky Parameterization)  #
###########################################################################

dmvnpc <- function(x, mu, U, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     Omega <- t(U) %*% U
     logdetOmega <- logdet(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector((-k/2)*log(2*pi) + 0.5*logdetOmega - 0.5*z)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvnpc <- function(n=1, mu=rep(0,k), U)
     {
     mu <- rbind(mu)
     if(missing(U)) stop("Upper triangular U is required.")
     Sigma <- as.inverse(t(U) %*% U)
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     z <- matrix(rnorm(n*k),n,k) %*% chol(Sigma)
     x <- mu + z
     return(x)
     }

###########################################################################
# Multivariate Polya Distribution                                         #
###########################################################################

dmvpolya <- function(x, alpha, log=FALSE)
     {
     x <- as.vector(x)
     alpha <- as.vector(alpha)
     if(!identical(length(x), length(alpha)))
          stop("x and alpha differ in length.")
     dens <- (log(factorial(sum(x))) - sum(log(factorial(x)))) +
          (log(factorial(sum(alpha)-1)) -
          log(factorial(sum(x) + sum(alpha)-1))) +
          (sum(log(factorial(x + alpha - 1)) - log(factorial(alpha - 1))))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvpolya <- function(n=1, alpha)
     {
     p <- rdirichlet(n, alpha)
     x <- rcat(n,p)
     return(x)
     }

###########################################################################
# Multivariate Power Exponential Distribution                             #
###########################################################################

dmvpe <- function(x=c(0,0), mu=c(0,0), Sigma=diag(2), kappa=1, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     temp <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(((log(k)+lgamma(k/2)) - ((k/2)*log(pi) +
          0.5*logdet(Sigma) + lgamma(1 + k/(2*kappa)) +
          (1 + k/(2*kappa))*log(2))) + kappa*(-0.5*temp))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvpe <- function(n, mu=c(0,0), Sigma=diag(2), kappa=1)
     {
     mu <- rbind(mu)
     k <- ncol(mu)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     if(k != nrow(Sigma)) {
          stop("mu and Sigma have non-conforming size.")}
     ev <- eigen(Sigma, symmetric=TRUE)
     if(!all(ev$values >= -sqrt(.Machine$double.eps) *
          abs(ev$values[1]))) {
          stop("Sigma must be positive-definite.")}
     SigmaSqrt <- ev$vectors %*% diag(sqrt(ev$values),
          length(ev$values)) %*% t(ev$vectors)
     radius <- (rgamma(n, shape=k/(2*kappa), scale=1/2))^(1/(2*kappa))
     runifsphere <- function(n, k)
          {
          p <- as.integer(k)
          if(!is.integer(k))
               stop("k must be an integer in [2,Inf)")
          if(k < 2) stop("k must be an integer in [2,Inf).")
          Mnormal <- matrix(rnorm(n*k,0,1), nrow=n)
          rownorms <- sqrt(rowSums(Mnormal^2))
          unifsphere <- sweep(Mnormal,1,rownorms, "/")
          return(unifsphere)
          }
     un <- runifsphere(n=n, k=k)
     x <- mu + radius * un %*% SigmaSqrt
     return(x)
     }

###########################################################################
# Multivariate Power Exponential Distribution (Cholesky Parameterization) #
###########################################################################

dmvpec <- function(x=c(0,0), mu=c(0,0), U, kappa=1, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     Sigma <- t(U) %*% U
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     temp <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(((log(k)+lgamma(k/2)) - ((k/2)*log(pi) +
          0.5*logdet(Sigma) + lgamma(1 + k/(2*kappa)) +
          (1 + k/(2*kappa))*log(2))) + kappa*(-0.5*temp))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvpec <- function(n, mu=c(0,0), U, kappa=1)
     {
     mu <- rbind(mu)
     k <- ncol(mu)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     if(k != nrow(U)) {
          stop("mu and U have non-conforming size.")}
     Sigma <- t(U) %*% U
     ev <- eigen(Sigma, symmetric=TRUE)
     if(!all(ev$values >= -sqrt(.Machine$double.eps) *
          abs(ev$values[1]))) {
          stop("Sigma must be positive-definite.")}
     SigmaSqrt <- ev$vectors %*% diag(sqrt(ev$values),
          length(ev$values)) %*% t(ev$vectors)
     radius <- (rgamma(n, shape=k/(2*kappa), scale=1/2))^(1/(2*kappa))
     runifsphere <- function(n, k)
          {
          p <- as.integer(k)
          if(!is.integer(k))
               stop("k must be an integer in [2,Inf)")
          if(k < 2) stop("k must be an integer in [2,Inf).")
          Mnormal <- matrix(rnorm(n*k,0,1), nrow=n)
          rownorms <- sqrt(rowSums(Mnormal^2))
          unifsphere <- sweep(Mnormal,1,rownorms, "/")
          return(unifsphere)
          }
     un <- runifsphere(n=n, k=k)
     x <- mu + radius * un %*% SigmaSqrt
     return(x)
     }

###########################################################################
# Multivariate t Distribution                                             #
###########################################################################

dmvt <- function(x, mu, S, df=Inf, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(S)) S <- diag(ncol(x))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(df > 10000)) return(dmvn(x, mu, S, log))
     k <- nrow(S)
     ss <- x - mu
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((df+k)/2) - lgamma(df/2) - (k/2)*log(df) -
          (k/2)*log(pi) - 0.5*logdet(S) - ((df+k)/2)*log(1 + (1/df) * z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvt <- function(n=1, mu=rep(0,k), S, df=Inf)
     {
     mu <- rbind(mu)
     if(missing(S)) S <- diag(ncol(mu))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     k <- ncol(S)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     if(df==Inf) x <- 1 else x <- rchisq(n,df) / df
     x[which(x == 0)] <- 1e-100
     z <- rmvn(n, rep(0,k), S)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Multivariate t Distribution (Cholesky Parameterization)                 #
###########################################################################

dmvtc <- function(x, mu, U, df=Inf, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(df > 10000)) return(dmvnc(x, mu, U, log))
     k <- nrow(U)
     ss <- x - mu
     S <- t(U) %*% U
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((df+k)/2) - lgamma(df/2) + (k/2)*df +
          (k/2)*log(pi) + 0.5*logdet(S) + ((df+k)/2)*log(1 + (1/df) * z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvtc <- function(n=1, mu=rep(0,k), U, df=Inf)
     {
     mu <- rbind(mu)
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     k <- ncol(U)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     if(df==Inf) x <- 1 else x <- rchisq(n,df) / df
     x[which(x == 0)] <- 1e-100
     z <- rmvnc(n, rep(0,k), U)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Multivariate t Distribution (Precision Parameterization)                #
###########################################################################

dmvtp <- function(x, mu, Omega, nu=Inf, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     if(any(nu > 10000)) return(dmvnp(x, mu, Omega, log))
     k <- ncol(Omega)
     logdetOmega <- logdet(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((nu+k)/2) - (lgamma(nu/2) + (k/2)*log(nu) +
          (k/2)*log(pi)) + 0.5*logdetOmega +
          (-(nu+k)/2)*log(1 + (1/nu) * z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvtp <- function(n=1, mu, Omega, nu=Inf)
     {
     mu <- rbind(mu)
     if(missing(Omega)) Omega <- diag(ncol(mu))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     Sigma <- as.inverse(Omega)
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     if(nu == Inf) x <- 1 else x <- rchisq(n,nu) / nu
     x[which(x == 0)] <- 1e-100
     z <- rmvn(n, rep(0,k), Sigma)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Multivariate t Distribution (Precision-Cholesky Parameterization)       #
###########################################################################

dmvtpc <- function(x, mu, U, nu=Inf, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- matrix(mu, nrow(x), ncol(x), byrow=TRUE)
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     if(any(nu > 10000)) return(dmvnpc(x, mu, U, log))
     k <- ncol(U)
     Omega <- t(U) %*% U
     logdetOmega <- logdet(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((nu+k)/2) - (lgamma(nu/2) + (k/2)*log(nu) +
          (k/2)*log(pi)) + 0.5*logdetOmega +
          (-(nu+k)/2)*log(1 + (1/nu) * z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvtpc <- function(n=1, mu, U, nu=Inf)
     {
     mu <- rbind(mu)
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     k <- ncol(U)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     if(nu == Inf) x <- 1 else x <- rchisq(n,nu) / nu
     x[which(x == 0)] <- 1e-100
     z <- rmvnpc(n, rep(0,k), U)
     x <- mu + z/sqrt(x)
     return(x)
     }

###########################################################################
# Normal Distribution Mixture                                             #
###########################################################################

dnormm <- function(x, p, mu, sigma, log=FALSE)
     {
     if(missing(x)) stop("x is a required argument.")
     x <- as.vector(x)
     n <- length(x)
     if(missing(p)) stop("p is a required argument.")
     p <- as.vector(p)
     if(any(p <= 0) | any(p > 1)) stop("p must be in (0,1].")
     if(sum(p) != 1) stop("p must sum to 1 for all components.")
     m <- length(p)
     p <- matrix(p, n, m, byrow=TRUE)
     if(missing(mu)) stop("mu is a required argument.")
     mu <- as.vector(mu)
     if(!identical(m, length(mu)))
          stop("p and mu differ in length.")
     mu <- matrix(mu, n, m, byrow=TRUE)
     if(missing(sigma)) stop("sigma is a required argument.")
     sigma <- as.vector(sigma)
     if(!identical(m, length(sigma)))
          stop("p and sigma differ in length.")
     sigma <- matrix(sigma, n, m, byrow=TRUE)
     dens <- matrix(dnorm(x, mu, sigma, log=TRUE), n, m)
     dens <- dens + log(p)
     if(log == TRUE) dens <- apply(dens, 1, logadd)
     else dens <- rowSums(exp(dens))
     return(dens)
     }
pnormm <- function(q, p, mu, sigma, lower.tail=TRUE, log.p=FALSE)
     {
     n <- length(q)
     m <- length(p)
     q <- matrix(q, n, m)
     p <- matrix(p, n, m, byrow=TRUE)
     mu <- matrix(mu, n, m, byrow=TRUE)
     sigma <- matrix(sigma, n, m, byrow=TRUE)
     cdf <- matrix(pnorm(q, mu, sigma, lower.tail=lower.tail,
          log.p=log.p), n, m)
     if(log.p == FALSE) cdf <- rowSums(cdf * p)
     else stop("The log.p argument does not work yet.")
     return(cdf)
     }
rnormm <- function(n, p, mu, sigma)
     {
     if(missing(p)) stop("p is a required argument.")
     p <- as.vector(p)
     if(any(p <= 0) | any(p > 1)) stop("p must be in (0,1].")
     if(sum(p) != 1) stop("p must sum to 1 for all components.")
     m <- length(p)
     p <- matrix(p, n, m, byrow=TRUE)
     if(missing(mu)) stop("mu is a required argument.")
     mu <- as.vector(mu)
     if(!identical(m, length(mu)))
          stop("p and mu differ in length.")
     if(missing(sigma)) stop("sigma is a required argument.")
     sigma <- as.vector(sigma)
     if(!identical(m, length(sigma)))
          stop("p and sigma differ in length.")
     if(any(sigma <= 0)) stop("sigma must be positive.")
     z <- rcat(n, p)
     x <- rnorm(n, mean=mu[z], sd=sigma[z])
     return(x)
     }

###########################################################################
# Normal Distribution (Precision Parameterization)                        #
###########################################################################

dnormp <- function(x, mean=0, prec=1, log=FALSE)
     {
     #dens <- sqrt(prec/(2*pi)) * exp(-(prec/2)*(x-mu)^2)
     dens <- dnorm(x, mean, sqrt(1/prec), log)
     return(dens)
     }
pnormp <- function(q, mean=0, prec=1, lower.tail=TRUE, log.p=FALSE)
     {return(pnorm(q, mean=mean, sd=sqrt(1/prec), lower.tail, log.p))}
qnormp <- function(p, mean=0, prec=1, lower.tail=TRUE, log.p=FALSE)
     {return(qnorm(p, mean=mean, sd=sqrt(1/prec), lower.tail, log.p))}
rnormp <- function(n, mean=0, prec=1)
     {return(rnorm(n, mean=mean, sd=sqrt(1/prec)))}

###########################################################################
# Normal Distribution (Variance Parameterization)                         #
###########################################################################

dnormv <- function(x, mean=0, var=1, log=FALSE)
     {
     #dens <- (1/(sqrt(2*pi*var))) * exp(-((x-mu)^2/(2*var)))
     dens <- dnorm(x, mean, sqrt(var), log)
     return(dens)
     }
pnormv <- function(q, mean=0, var=1, lower.tail=TRUE, log.p=FALSE)
     {return(pnorm(q, mean=mean, sd=sqrt(var), lower.tail, log.p))}
qnormv <- function(p, mean=0, var=1, lower.tail=TRUE, log.p=FALSE)
     {return(qnorm(p, mean=mean, sd=sqrt(var), lower.tail, log.p))}
rnormv <- function(n, mean=0, var=1)
     {return(rnorm(n, mean=mean, sd=sqrt(var)))}

###########################################################################
# Normal-Inverse-Wishart Distribution                                     #
###########################################################################

dnorminvwishart <- function(mu, mu0, lambda, Sigma, S, nu, log=FALSE)
     {
     dens <- dinvwishart(Sigma, nu, S, log=TRUE) +
          dmvn(mu, mu0, 1/lambda*Sigma, log=TRUE)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rnorminvwishart <- function(n=1, mu0, lambda, S, nu)
     {
     Sigma <- rinvwishart(nu, S)
     mu <- rmvn(n, mu0, 1/lambda*Sigma)
     return(list(mu=mu, Sigma=Sigma))
     }

###########################################################################
# Normal-Laplace Distribution                                             #
###########################################################################

dnormlaplace <- function(x, mu=0, sigma=1, alpha=1, beta=1, log=FALSE)
     {
     x <- as.vector(x)
     mu <- as.vector(mu)
     sigma <- as.vector(sigma)
     alpha <- as.vector(alpha)
     beta <- as.vector(beta)
     if(any(sigma <= 0))
          stop("The sigma parameter must be positive.")
     if(any(alpha <= 0))
          stop("The alpha parameter must be positive.")
     if(any(beta <= 0))
          stop("The beta parameter must be positive.")
    NN <- max(length(x), length(mu), length(sigma), length(alpha),
          length(beta))
    x <- rep(x, len=NN)
    mu <- rep(mu, len=NN)
    sigma <- rep(sigma, len=NN)
    alpha <- rep(alpha, len=NN)
    beta <- rep(beta, len=NN)
     a <- dnorm((x - mu) / sigma, log=TRUE)
     z <- alpha*sigma - (x - mu) / sigma
     b <- pnorm(z, lower.tail=FALSE, log.p=TRUE) -
          dnorm(z, log=TRUE)
     z <- beta*sigma + (x - mu) / sigma
     cc <- pnorm(z, lower.tail=FALSE, log.p=TRUE) -
          dnorm(z, log=TRUE)
     z <- log(alpha*beta / (alpha + beta))
     d <- pnorm(z, lower.tail=FALSE, log.p=TRUE) -
          dnorm(z, log=TRUE)
     if(log == FALSE) dens <- exp(d + a + b) + exp(d + a + cc)
     else {
          dens <- rep(NA, NN)
          e <- d + a + b
          f <- d + a + cc
          for (i in 1:NN) dens[i] <- logadd(c(e[i], f[i]))
      }
     return(dens)
     }
rnormlaplace <- function(n, mu=0, sigma=1, alpha=1, beta=1)
     {
     z <- rnorm(n, 0, sigma^2)
     w <- rslaplace(n, mu, 1/beta, 1/alpha)
     return(z + w)
     }

###########################################################################
# Normal-Wishart Distribution                                             #
###########################################################################

dnormwishart <- function(mu, mu0, lambda, Omega, S, nu, log=FALSE)
     {
     dens <- dwishart(Omega, nu, S, log=TRUE) +
          dmvnp(mu, mu0, lambda*Omega, log=TRUE)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rnormwishart <- function(n=1, mu0, lambda, S, nu)
     {
     Omega <- rwishart(nu, S)
     mu <- rmvn(n, mu0, as.inverse(lambda*Omega))
     return(list(mu=mu, Omega=Omega))
     }

###########################################################################
# Pareto Distribution                                                     #
###########################################################################

dpareto <- function(x, alpha, log=FALSE)
     {
     x <- as.vector(x); alpha <- as.vector(alpha)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(x), length(alpha))
     x <- dens <- rep(x, len=NN); alpha <- rep(alpha, len=NN)
     dens[which(x < 1)] <- -Inf
     temp <- which(x >= 1)
     dens[temp] <- log(alpha[temp]) - (alpha[temp] + 1)*log(x[temp])
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
ppareto <- function(q, alpha)
     {
     q <- as.vector(q); alpha <- as.vector(alpha)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(q), length(alpha))
     p <- q <- rep(q, len=NN); alpha <- rep(alpha, len=NN)
     p[which(q < 1)] <- 0
     temp <- which(q >= 1)
     p[temp] <- 1 - 1 / q[temp]^alpha[temp]
     return(p)
     }
qpareto <- function(p, alpha)
     {
     p <- as.vector(p); alpha <- as.vector(alpha)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(p), length(alpha))
     p <- rep(p, len=NN); alpha <- rep(alpha, len=NN)
     q <- (1-p)^(-1/alpha)
     return(q)
     }
rpareto <- function(n, alpha)
     {
     alpha <- rep(alpha, len=n)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     x <- runif(n)^(-1/alpha)
     return(x)
     }

###########################################################################
# Power Exponential Distribution                                          #
#                                                                         #
# These functions are similar to those in the normalp package.            #
###########################################################################

dpe <- function(x, mu=0, sigma=1, kappa=2, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); sigma <- as.vector(sigma)
     kappa <- as.vector(kappa)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(x), length(mu), length(sigma), length(kappa))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); kappa <- rep(kappa, len=NN)
     cost <- 2 * kappa^(1/kappa) * gamma(1 + 1/kappa) * sigma
     expon1 <- (abs(x - mu))^kappa
     expon2 <- kappa * sigma^kappa
     dens <- log(1/cost) + (-expon1 / expon2)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
ppe <- function(q, mu=0, sigma=1, kappa=2, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu); sigma <- as.vector(sigma)
     kappa <- as.vector(kappa)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(q), length(mu), length(sigma), length(kappa))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); kappa <- rep(kappa, len=NN)
     z <- (q - mu) / sigma
     zz <- abs(z)^kappa
     p <- pgamma(zz, shape=1/kappa, scale=kappa)
     p <- p / 2
     temp <- which(z < 0); p[temp] <- 0.5 - p[temp]
     temp <- which(z >= 0); p[temp] <- 0.5 + p[temp]
     if(lower.tail == FALSE) p <- 1 - p
     if(log.p == TRUE) p <- log(p)
     return(p)
     }
qpe <- function(p, mu=0, sigma=1, kappa=2, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu); sigma <- as.vector(sigma)
     kappa <- as.vector(kappa)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(p), length(mu), length(sigma), length(kappa))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); kappa <- rep(kappa, len=NN)
     if(log.p == TRUE) p <- log(p)
     if(lower.tail == FALSE) p <- 1 - p
     zp <- 0.5 - p
     temp <- which(p >= 0.5); zp[temp] <- p[temp] - 0.5
     zp <- 2 * zp
     qg <- qgamma(zp, shape=1/kappa, scale=kappa)
     z <- qg^(1/kappa)
     temp <- which(p < 0.5); z[temp] <- -z[temp]
     q <- mu + z * sigma
     return(q)
     }
rpe <- function(n, mu=0, sigma=1, kappa=2)
     {
     mu <- rep(mu, len=n); sigma <- rep(sigma, len=n)
     kappa <- rep(kappa, len=n)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     qg <- rgamma(n, shape=1/kappa, scale=kappa)
     z <- qg^(1/kappa)
     u <- runif(n)
     temp <- which(u < 0.5); z[temp] <- -z[temp]
     x <- mu + z * sigma
     return(x)
     }

###########################################################################
# Scaled Inverse Wishart Distribution                                     #
###########################################################################

dsiw <- function(Q, nu, S, zeta, mu, delta, log=FALSE)
     {
     dens <- dinvwishart(Q, nu, S, log=TRUE) +
          sum(dmvn(log(zeta), mu, diag(delta), log=TRUE))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rsiw <- function(nu, S, mu, delta)
     {
     Q <- rinvwishart(nu, S)
     Zeta <- diag(as.vector(exp(rmvn(1, mu, diag(delta)))))
     x <- Zeta %*% Q %*% Zeta
     return(x)
     }

###########################################################################
# Skew Discrete Laplace Distribution                                      #
###########################################################################

dsdlaplace <- function(x, p, q, log=FALSE)
     {
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(q < 0) || any(q > 1)) stop("q must be in [0,1].")
     NN <- max(length(x), length(p), length(q))
     x <- rep(x, len=NN); p <- rep(p, len=NN); q <- rep(q, len=NN)
     dens <- log(1-p) + log(1-q) - (log(1-p*q) + x*log(p))
     temp <- which(x < 0)
     dens[temp] <- log(1-p[temp]) + log(1-q[temp]) -
          (log(1-p[temp]*q[temp]) + abs(x[temp])*log(q[temp]))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
psdlaplace <- function(x, p, q)
     {
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(q < 0) || any(q > 1)) stop("q must be in [0,1].")
     NN <- max(length(x), length(p), length(q))
     x <- rep(x, len=NN); p <- rep(p, len=NN); q <- rep(q, len=NN)
     pr <- 1-(1-q)*p^(floor(x)+1)/(1-p*q)
     temp <- which(x < 0)
     pr[temp] <- (1-p[temp])*q[temp]^(-floor(x[temp]))/(1-p[temp]*q[temp])
     return(pr)
     }
qsdlaplace <- function(prob, p, q)
     {
     if(any(prob < 0) || any(prob > 1)) stop("prob must be in [0,1].")
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(q < 0) || any(q > 1)) stop("q must be in [0,1].")
     NN <- max(length(prob), length(p), length(q))
     prob <- rep(prob, len=NN); p <- rep(p, len=NN); q <- rep(q, len=NN)
     x <- numeric(NN)
     for (i in 1:NN) {
          k <- 0
          if(prob[i] >= psdlaplace(k, p[i], q[i])) {
               while(prob[i] >= psdlaplace(k, p[i], q[i])) {
                    k <- k + 1}}
          else if(prob[i] < psdlaplace(k, p[i], q[i])) {
               while(prob[i] < psdlaplace(k, p[i], q[i])) {
                    k <- k - 1}
               k <- k + 1}
          x[i] <- k
          }
     return(x)
     }
rsdlaplace <- function(n, p, q)
     {
     if(length(p) > 1) stop("p must have a length of 1.")
     if(length(q) > 1) stop("q must have a length of 1.")
     if((p < 0) || (p > 1)) stop("p must be in [0,1].")
     if((q < 0) || (q > 1)) stop("q must be in [0,1].")
     u <- runif(n)
     return(qsdlaplace(u,p,q))
     }

###########################################################################
# Skew-Laplace Distribution                                               #
###########################################################################

dslaplace <- function(x, mu, alpha, beta, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu)
     alpha <- as.vector(alpha); beta <- as.vector(beta)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     NN <- max(length(x), length(mu), length(alpha), length(beta))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     alpha <- rep(alpha, len=NN); beta <- rep(beta, len=NN)
     ab <- alpha + beta
     dens <- belowMu <- log(1/ab) + ((x - mu)/alpha)
     aboveMu <- log(1/ab) + ((mu - x)/beta)
     temp <- which(x > mu); dens[temp] <- aboveMu[temp]
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
pslaplace <- function(q, mu, alpha, beta)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     alpha <- as.vector(alpha); beta <- as.vector(beta)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     NN <- max(length(q), length(mu), length(alpha), length(beta))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     alpha <- rep(alpha, len=NN); beta <- rep(beta, len=NN)
     ab <- alpha + beta
     p <- belowMu <- (alpha/ab) * exp((q - mu)/alpha)
     aboveMu <- 1 - (beta/ab) * exp((mu - q)/beta)
     temp <- which(q >= mu); p[temp] <- aboveMu[temp]
     return(p)
     }
qslaplace <- function(p, mu, alpha, beta)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     alpha <- as.vector(alpha); beta <- as.vector(beta)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     NN <- max(length(p), length(mu), length(alpha), length(beta))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     alpha <- rep(alpha, len=NN); beta <- rep(beta, len=NN)
     ab <- alpha + beta
     q <- belowMu <- alpha*log(p*ab/alpha) + mu
     aboveMu <- mu - beta*log(ab*(1 - p)/beta)
     temp <- which(p >= alpha/ab); q[temp] <- aboveMu[temp]
     return(q)
     }
rslaplace <- function(n, mu, alpha, beta)
     {
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     mu <- rep(mu, len=n)
     alpha <- rep(alpha, len=n)
     beta <- rep(beta, len=n)
     y <- rexp(n, 1)
     probs <- alpha / (alpha + beta)
     signs <- runif(n)
     temp <- which(signs <= probs)
     signs[temp] <- -1
     signs[-temp] <- 1
     mult <- signs * beta
     temp <- which(signs < 0)
     mult[temp] <- signs[temp] * alpha[temp]
     x <- mult * y + mu
     return(x)
     }

###########################################################################
# Stick-Breaking Prior Distribution                                       #
###########################################################################

dStick <- function(theta, gamma, log=FALSE)
     {
     if(log == FALSE) dens <- prod(dbeta(theta, 1, gamma, log=FALSE))
     else dens <- sum(dbeta(theta, 1, gamma, log=TRUE))
     return(dens)
     }
rStick <- function(M, gamma)
     {
     betas <- rbeta(M+1, 1, gamma)
     remaining <- c(1, cumprod(1 - betas))[1:(M+1)]
     w <- remaining * betas
     w <- w / sum(w)
     return(w)
     #return(Stick(rbeta(M,1,gamma))) #Deprecated
     }

###########################################################################
# Student t Distribution (3-parameter)                                    #
#                                                                         #
# The pst and qst functions are similar to the TF functions in the        #
# gamlss.dist package, but dst and rst have been refined.                 #
###########################################################################

dst <- function(x, mu=0, sigma=1, nu=10, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu)
     sigma <- as.vector(sigma); nu <- as.vector(nu)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     else if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(x), length(mu), length(sigma), length(nu))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); nu <- rep(nu, len=NN)
     const <- lgamma((nu+1)/2) - lgamma(nu/2) - log(sqrt(pi*nu) * sigma)
     dens <- const + log((1 + (1/nu)*((x-mu)/sigma)^2)^(-(nu+1)/2))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     #return({1/sigma} * dt({x-mu}/sigma, df=nu, log)) #Deprecated
     }
pst <- function(q, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     sigma <- as.vector(sigma); nu <- as.vector(nu)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(q), length(mu), length(sigma), length(nu))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); nu <- rep(nu, len=NN)
     p <- pt({q-mu}/sigma, df=nu, lower.tail=lower.tail, log.p=log.p)
     temp <- which(nu > 1e6)
     p[temp] <- pnorm(q[temp], mu[temp], sigma[temp],
          lower.tail=lower.tail, log.p=log.p)
     return(p)
     }
qst <- function(p, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     sigma <- as.vector(sigma); nu <- as.vector(nu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(p), length(mu), length(sigma), length(nu))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); nu <- rep(nu, len=NN)
     q <- mu + sigma * qt(p, df=nu, lower.tail=lower.tail)
     temp <- which(nu > 1e6)
     q[temp] <- qnorm(p[temp], mu[temp], sigma[temp],
          lower.tail=lower.tail, log.p=log.p)
     return(q)
     }
rst <- function(n, mu=0, sigma=1, nu=10)
     {
     mu <- rep(mu, len=n); sigma <- rep(sigma, len=n); nu <- rep(nu, len=n)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     n <- ceiling(n)
     y <- rnorm(n)
     z <- rchisq(n, nu)
     x <- mu + sigma*y*sqrt(nu/z)
     return(x)
     }

###########################################################################
# Student t Distribution (Precision Parameterization)                     #
###########################################################################

dstp <- function(x, mu=0, tau=1, nu=10, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu)
     tau <- as.vector(tau); nu <- as.vector(nu)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(x), length(mu), length(tau), length(nu))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     tau <- rep(tau, len=NN); nu <- rep(nu, len=NN)
     dens <- (lgamma((nu+1)/2) - lgamma(nu/2)) + 0.5*log(tau/(nu*pi)) +
          (-(nu+1)/2)*log(1 + (tau/nu)*(x-mu)^2)
     temp <- which(nu > 1e6)
     dens[temp] <- dnorm(x[temp], mu[temp], sqrt(1/tau[temp]), log=TRUE)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
pstp <- function(q, mu=0, tau=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     tau <- as.vector(tau); nu <- as.vector(nu)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(q), length(mu), length(tau), length(nu))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     tau <- rep(tau, len=NN); nu <- rep(nu, len=NN)
     p <- pst(q, mu, sqrt(1/tau), nu, lower.tail, log.p)
     return(p)
     }
qstp <- function(p, mu=0, tau=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     tau <- as.vector(tau); nu <- as.vector(nu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(p), length(mu), length(tau), length(nu))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     tau <- rep(tau, len=NN); nu <- rep(nu, len=NN)
     q <- qst(p, mu, sqrt(1/tau), nu, lower.tail, log.p)
     return(q)
     }
rstp <- function(n, mu=0, tau=1, nu=10)
     {
     mu <- rep(mu, len=n); tau <- rep(tau, len=n); nu <- rep(nu, len=n)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     x <- rst(n, mu, sqrt(1/tau), nu)
     n <- ceiling(n)
     p <- runif(n)
     x <- qst(p, mu=mu, sqrt(1/tau), nu=nu)
     return(x)
     }

###########################################################################
# Truncated Distribution                                                  #
#                                                                         #
# These functions are similar to those from Nadarajah, S. and Kotz, S.    #
# (2006). ``R Programs for Computing Truncated Distributions''. Journal   #
# of Statistical Software, 16, Code Snippet 2, 1-8. These functions have  #
# been corrected to work with log-densities.                              #
###########################################################################

dtrunc <- function(x, spec, a=-Inf, b=Inf, log=FALSE, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     if(any(x < a) | any(x > b))
          stop("At least one instance of (x < a) or (x > b) found.")
     dens <- rep(0, length(x))
     g <- get(paste("d", spec, sep=""), mode="function")
     G <- get(paste("p", spec, sep=""), mode="function")
     if(log == TRUE) {
          dens <- g(x, log=TRUE, ...) - log(G(b, ...) - G(a, ...))
          }
     else {
          dens <- g(x, ...) / (G(b, ...) - G(a, ...))}
     return(dens)
     }
extrunc <- function(spec, a=-Inf, b=Inf, ...)
     {
     f <- function(x) x * dtrunc(x, spec, a=a, b=b, log=FALSE, ...)
     return(integrate(f, lower=a, upper=b)$value)
     }
ptrunc <- function(x, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     if(any(x < a) | any(x > b))
          stop("At least one instance of (x < a) or (x > b) found.")
     p <- x
     aa <- rep(a, length(x))
     bb <- rep(b, length(x))
     G <- get(paste("p", spec, sep=""), mode="function")
     p <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
     p <- p - G(aa, ...)
     p <- p / {G(bb, ...) - G(aa, ...)}
     return(p)
     }
qtrunc <- function(p, spec, a=-Inf, b=Inf, ...)
     {
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     q <- p
     G <- get(paste("p", spec, sep=""), mode="function")
     Gin <- get(paste("q", spec, sep=""), mode="function")
     q <- Gin(G(a, ...) + p*{G(b, ...) - G(a, ...)}, ...)
     return(q)
     }
rtrunc <- function(n, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     x <- u <- runif(n)
     x <- qtrunc(u, spec, a=a, b=b,...)
     return(x)
     }
vartrunc <- function(spec, a=-Inf, b=Inf, ...)
     {
     ex <- extrunc(spec, a=a, b=b, ...)
     f <- function(x) {
          {x - ex}^2 * dtrunc(x, spec, a=a, b=b, log=FALSE, ...)}
     sigma2 <- integrate(f, lower=a, upper=b)$value
     return(sigma2)
     }

###########################################################################
# Wishart Distribution                                                    #
###########################################################################

dwishart <- function(Omega, nu, S, log=FALSE)
     {
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(!identical(dim(Omega), dim(S)))
          stop("The dimensions of Omega and S differ.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     k <- nrow(Omega)
     gamsum <- 0
     for (i in 1:k) {gamsum <- gamsum + lgamma((nu + 1 - i)/2)}
     dens <- -((nu*k)/2) * log(2) - ((k*(k - 1))/4) * log(pi) - gamsum -
          (nu/2) * log(det(S)) + ((nu - k - 1)/2) * logdet(Omega) -
          (tr(as.inverse(S) %*% Omega)/2)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rwishart <- function(nu, S)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     k <- nrow(S)
     Z <- matrix(0, k, k)
     x <- rchisq(k, nu:{nu - k + 1})
     x[which(x == 0)] <- 1e-100
     diag(Z) <- sqrt(x)
     if(k > 1) {
          kseq <- 1:(k-1)
          Z[rep(k*kseq, kseq) +
               unlist(lapply(kseq, seq))] <- rnorm(k*{k - 1}/2)}
     return(crossprod(Z %*% chol(S)))
     }

###########################################################################
# Wishart Distribution (Cholesky Parameterization)                        #
###########################################################################

dwishartc <- function(U, nu, S, log=FALSE)
     {
     if(missing(U)) stop("Upper triangular U is required.")
     Omega <- t(U) %*% U
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(!identical(dim(Omega), dim(S)))
          stop("The dimensions of Omega and S differ.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     k <- nrow(Omega)
     gamsum <- 0
     for (i in 1:k) {gamsum <- gamsum + lgamma((nu + 1 - i)/2)}
     dens <- -((nu*k)/2) * log(2) - ((k*(k - 1))/4) * log(pi) - gamsum -
          (nu/2) * log(det(S)) + ((nu - k - 1)/2) * logdet(Omega) -
          (tr(as.inverse(S) %*% Omega)/2)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rwishartc <- function(nu, S)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(nu < nrow(S)) {
          stop("The nu parameter is less than the dimension of S.")}
     k <- nrow(S)
     Z <- matrix(0, k, k)
     x <- rchisq(k, nu:{nu - k + 1})
     x[which(x == 0)] <- 1e-100
     diag(Z) <- sqrt(x)
     if(k > 1) {
          kseq <- 1:(k-1)
          Z[rep(k*kseq, kseq) +
               unlist(lapply(kseq, seq))] <- rnorm(k*{k - 1}/2)}
     return(Z %*% chol(S))
     }

###########################################################################
# Yang-Berger Distribution                                                #
###########################################################################

dyangberger <- function(x, log=FALSE)
     {
     if(missing(x)) stop("Matrix x is a required argument.")
     if(!is.matrix(x)) x <- as.matrix(x)
     x <- as.symmetric.matrix(x)
     if(!is.positive.definite(x))
          stop("Matrix x is not positive-definite.")
     d <- sort(eigen(x)$values)
     dens <- log(1) - logdet(x) * prod(diff(d))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
dyangbergerc <- function(x, log=FALSE)
     {
     if(missing(x))
          stop("Upper triangular matrix x is a required argument.")
     if(!is.matrix(x)) x <- as.matrix(x)
     x <- t(x) %*% x
     d <- sort(eigen(x)$values)
     dens <- log(1) - logdet(x) * prod(diff(d))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }

###########################################################################
# Zellner's g-Prior                                                       #
###########################################################################

dhyperg <- function(g, alpha=3, log=FALSE)
     {
     if(g <= 0) stop("The g parameter must be positive.")
     if(alpha <= 0) stop("The alpha parameter must be positive.")
     dens <- log((alpha - 2)/2) -(alpha/2)*log(1 + g)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
dzellner <- function(beta, g, sigma, X, log=FALSE)
     {
     if(g <= 0) stop("The g parameter must be positive.")
     if(sigma <= 0) stop("The sigma parameter must be positive.")
  ## HS (01/2020): The original formulation was:
  ### g*sigma*sigma* as.inverse(t(X) %*% X)
  ### the latter bit is equivalent to: solve(crossprod(X))
  ### which should be equivalent to the more stable and accurate
  ### chol2inv(qr.R(qr(X))) 
  ## see: strucchange::solveCrossprod
     dens <- dmvn(beta, rep(0, length(beta)),
          g*sigma*sigma*chol2inv(qr.R(qr(X))), log=log)
     return(dens)
     }
rzellner <- function(n, g, sigma, X)
     {
     x <- rmvn(n, rep(0, ncol(X)), g*sigma*sigma*chol2inv(qr.R(qr(X))))
     return(x)
     }

#End
