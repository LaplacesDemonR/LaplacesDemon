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
     exponent <- -(sqrt(2) / scale) * abs(x - location) *
          ifelse(x >= location, kappa, 1/kappa)
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
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     exponent <- -(sqrt(2) / scale) * abs(q - location) *
          ifelse(q >= location, kappa, 1/kappa)
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
     exponent <- ifelse(x >= Delta, -(Alpha+1),
          (Beta-1)) * (log(x) - location)
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

dcat <- function(x, p, log=FALSE)
     {
     if(is.vector(x) & !is.matrix(p))
          p <- matrix(p, length(x), length(p), byrow=TRUE)
     if(is.matrix(x) & !is.matrix(p))
          p <- matrix(p, nrow(x), length(p), byrow=TRUE)
     if(is.vector(x) & {length(x) == 1}) {
          temp <- rep(0, ncol(p))
          temp[x] <- 1
          x <- t(temp)}
     else if(is.vector(x) & (length(x) > 1))
          x <- as.indicator.matrix(x)
     if(!identical(nrow(x), nrow(p)))
          stop("The number of rows of x and p differ.")
     if(!identical(ncol(x), ncol(p))) {
          x.temp <- matrix(0, nrow(p), ncol(p))
          x.temp[,as.numeric(colnames(x))] <- x
          x <- x.temp}
     dens <- x*log(p)
     if(log == FALSE) dens <- x*p
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
          x <- as.vector(which(rmultinom(n, size=1, prob=p) == 1,
               arr.ind=TRUE)[, "row"])}
     else {
          n <- nrow(p)
          x <- apply(p, 1, function(x) {
               as.vector(which(rmultinom(1, size=1, prob=x) == 1,
               arr.ind=TRUE)[, "row"])})
          }
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
     dens <- (-(nu+1)/2)*log(1 + (1/nu)*(x/scale)*(x/scale))
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

dhs <- function(x, lambda, tau, sigma, log=FALSE)
     {
     NN <- max(length(x), length(lambda))
     x <- rep(x, len=NN); lambda <- rep(lambda, len=NN)
     p.theta.lambda <- dnorm(x, 0, lambda, log=TRUE)
     p.lambda.tau <- dhalfcauchy(lambda, tau, log=TRUE)
     p.tau <- dhalfcauchy(tau, sigma, log=TRUE)
     dens <- p.theta.lambda - {p.lambda.tau + p.tau} / NN
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rhs <- function(n, lambda, tau, sigma)
     {
     if(missing(tau)) tau <- rhalfcauchy(1, sigma)
     if(missing(lambda)) lambda <- rhalfcauchy(n, tau)
     theta <- rnorm(n, 0, lambda)
     return(theta)
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
     dens <- nu*log(nu) - log(gamma(nu)) + nu*log(scale) -
          (nu+1)*log(x) - (nu*scale/x)
     if(log == FALSE) dens <- exp(dens)
     }

rinvchisq <- function(n, df, scale=1/df)
     {
     df <- rep(df, len=n); scale <- rep(scale, len=n)
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     z <- rchisq(n, df=df)
     z <- ifelse(z == 0, 1e-100, z)
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
     {return(1 / rgamma(n=n, shape=shape, rate=scale))}

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
     dens <- log(lambda / (2*pi*x^3)^0.5) -
          ((lambda*(x - mu)^2) / (2*mu^2*x))
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
     x <- ifelse(z > (mu / (mu+x)), mu^2/x, x)
     return(x)
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
          (nu/2)*log(det(S)) - ((nu + k + 1)/2)*log(det(Sigma)) -
          0.5*tr(S %*% as.inverse(Sigma))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rinvwishart <- function(nu, S)
     {return(as.inverse(rwishart(nu, as.inverse(S))))}

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
          (nu/2)*log(det(S)) - ((nu + k + 1)/2)*log(det(Sigma)) -
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
     p <- ifelse(q < location, 0.5 * exp(z), 1 - 0.5 * exp(-z))
     return(p)
     }
qlaplace <- function(p, location=0, scale=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(location), length(scale))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     q <- location - sign(p - 0.5) * scale * log(2 * ifelse(p < 0.5,
          p, 1 - p))
     return(q)
     }
rlaplace <- function(n, location=0, scale=1)
     {
     location <- rep(location, len=n); scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     r <- runif(n)
     x <- location - sign(r - 0.5) * scale * log(2 * ifelse(r < 0.5,
          r, 1 - r))
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
     exponent <- ifelse(x >= Delta, -(Alpha+1),
          (Beta-1)) * (log(x) - location)
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

dlnormp <- function(x, mu, tau, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(x), length(mu), length(tau))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     dens <- log(sqrt(tau/(2*pi))) + log(1/x) + (-(tau/2)*(log(x-mu)^2))
     if(log == FALSE) dens <- exp(dens)
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
rlnormp <- function(n, mu, tau)
     {
     mu <- rep(mu, len=n); tau <- rep(tau, len=n)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     x <- rnorm(n, mu, sqrt(1/tau))
     return(x)
     }

###########################################################################
# Multivariate Cauchy Distribution                                        #
###########################################################################

dmvc <- function(x, mu, S, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(S)) S <- diag(ncol(x))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     k <- nrow(S)
     ss <- x - mu
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma(k/2) - (lgamma(0.5) + log(1^(k/2)) +
          (k/2)*log(pi) + 0.5*log(det(S)) + ((1+k)/2)*log(1+z)))
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
     x <- ifelse(x == 0, 1e-100, x)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     k <- nrow(U)
     S <- t(U) %*% U
     ss <- x - mu
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma(k/2) - (lgamma(0.5) + log(1^(k/2)) +
          (k/2)*log(pi) + 0.5*log(det(S)) + ((1+k)/2)*log(1+z)))
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
     x <- ifelse(x == 0, 1e-100, x)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     k <- nrow(Omega)
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((1+k)/2) - (lgamma(0.5) + log(1^(k/2)) +
          (k/2)*log(pi)) + 0.5*log(detOmega) + (-(1+k)/2)*log(1 + z))
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
     x <- ifelse(x == 0, 1e-100, x)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     k <- nrow(U)
     Omega <- t(U) %*% U
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((1+k)/2) - (lgamma(0.5) + log(1^(k/2)) +
          (k/2)*log(pi)) + 0.5*log(detOmega) + (-(1+k)/2)*log(1 + z))
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
     x <- ifelse(x == 0, 1e-100, x)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma)) 
         stop("Matrix Sigma is not positive-definite.")
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- log(2 / ((2*pi)^(k/2) * sqrt(det(Sigma)))) +
          log((sqrt(pi / (2*sqrt(2*z))) * exp(-sqrt(2*z))) /
          sqrt(z/2)^(k/2 - 1))
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     Sigma <- t(U) %*% U
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- log(2 / ((2*pi)^(k/2) * sqrt(det(Sigma)))) +
          log((sqrt(pi / (2*sqrt(2*z))) * exp(-sqrt(2*z))) /
          sqrt(z/2)^(k/2 - 1))
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     d <- eigen(Sigma, symmetric=TRUE)$values
     dens <- as.vector(-0.5 * (k * log(2 * pi) + sum(log(d))) - (0.5 * z))
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     Sigma <- t(U) %*% U
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     d <- eigen(Sigma, symmetric=TRUE)$values
     dens <- as.vector(-0.5 * (k * log(2 * pi) + sum(log(d))) - (0.5 * z))
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     k <- nrow(Omega)
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector((-k/2)*(log(2) + log(pi)) + 0.5*log(detOmega) -
          0.5*z)
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
     z <- matrix(rnorm(n*k),n,k) %*% as.inverse(t(chol(Omega)))
     x <- mu + z
     return(x)
     }

###########################################################################
# Multivariate Normal Distribution (Precision-Cholesky Parameterization)  #
###########################################################################

dmvnpc <- function(x, mu, U, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     k <- ncol(U)
     Omega <- t(U) %*% U
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector((-k/2)*(log(2) + log(pi)) + 0.5*log(detOmega) -
          0.5*z)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
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
          0.5*log(det(Sigma)) + lgamma(1 + k/(2*kappa)) +
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     Sigma <- t(U) %*% U
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     temp <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(((log(k)+lgamma(k/2)) - ((k/2)*log(pi) +
          0.5*log(det(Sigma)) + lgamma(1 + k/(2*kappa)) +
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
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
     dens <- as.vector(lgamma((df+k)/2) - lgamma(df/2) + (k/2)*df +
          (k/2)*log(pi) + 0.5*log(det(S)) + ((df+k)/2)*log(1 + (1/df) * z))
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
     x <- ifelse(x == 0, 1e-100, x)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(df > 10000)) return(dmvnc(x, mu, U, log))
     k <- nrow(U)
     ss <- x - mu
     S <- t(U) %*% U
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((df+k)/2) - lgamma(df/2) + (k/2)*df +
          (k/2)*log(pi) + 0.5*log(det(S)) + ((df+k)/2)*log(1 + (1/df) * z))
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
     x <- ifelse(x == 0, 1e-100, x)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     if(any(nu > 10000)) return(dmvnp(x, mu, Omega, log))
     k <- ncol(Omega)
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((nu+k)/2) - (lgamma(nu/2) + (k/2)*log(nu) +
          (k/2)*log(pi)) + 0.5*log(detOmega) +
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
     x <- ifelse(x == 0, 1e-100, x)
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
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(U)) stop("Upper triangular U is required.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     if(any(nu > 10000)) return(dmvnpc(x, mu, U, log))
     k <- ncol(U)
     Omega <- t(U) %*% U
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(lgamma((nu+k)/2) - (lgamma(nu/2) + (k/2)*log(nu) +
          (k/2)*log(pi)) + 0.5*log(detOmega) +
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
     x <- ifelse(x == 0, 1e-100, x)
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
# Pareto Distribution                                                     #
###########################################################################

dpareto <- function(x, alpha, log=FALSE)
     {
     x <- as.vector(x); alpha <- as.vector(alpha)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(x), length(alpha))
     x <- rep(x, len=NN); alpha <- rep(alpha, len=NN)
     dens <- ifelse(x < 1, -Inf, log(alpha) - (alpha + 1)*log(x))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
ppareto <- function(q, alpha)
     {
     q <- as.vector(q); alpha <- as.vector(alpha)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(q), length(alpha))
     q <- rep(q, len=NN); alpha <- rep(alpha, len=NN)
     p <- ifelse(q < 1, 0, 1 - 1/q^alpha)
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
     p <- ifelse(z < 0, 0.5 - p, 0.5 + p)
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
     zp <- ifelse(p < 0.5, 0.5 - p, p - 0.5)
     zp <- 2 * zp
     qg <- qgamma(zp, shape=1/kappa, scale=kappa)
     z <- qg^(1/kappa)
     z <- ifelse(p < 0.5, -z, z)
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
     z <- ifelse(runif(n) < 0.5, -z, z)
     x <- mu + z * sigma
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
     dens <- ifelse(x >= 0, log(1-p) + log(1-q) - (log(1-p*q) +
          x*log(p)), log(1-p) + log(1-q) - (log(1-p*q) + abs(x)*log(q)))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
psdlaplace <- function(x, p, q)
     {
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(q < 0) || any(q > 1)) stop("q must be in [0,1].")
     NN <- max(length(x), length(p), length(q))
     x <- rep(x, len=NN); p <- rep(p, len=NN); q <- rep(q, len=NN)
     pr <- ifelse(x >= 0, 1-(1-q)*p^(floor(x)+1)/(1-p*q),
          (1-p)*q^(-floor(x))/(1-p*q))
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
     belowMu <- log(1/ab) + ((x - mu)/alpha)
     aboveMu <- log(1/ab) + ((mu - x)/beta)
     dens <- ifelse(x <= mu, belowMu, aboveMu)
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
     belowMu <- (alpha/ab) * exp((q - mu)/alpha)
     aboveMu <- 1 - (beta/ab) * exp((mu - q)/beta)
     p <- ifelse(q < mu, belowMu, aboveMu)
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
     belowMu <- alpha*log(p*ab/alpha) + mu
     aboveMu <- mu - beta*log(ab*(1 - p)/beta)
     q <- ifelse(p < alpha/ab, belowMu, aboveMu)
     return(q)
     }
rslaplace <- function(n, mu, alpha, beta)
     {
     #mu <- rep(mu, len=n); alpha <- rep(alpha, len=n)
     #beta <- rep(beta, len=n)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     ab <- alpha + beta
     y <- rexp(n,1)
     probs <- c(alpha,beta) / ab
     signs <- sample(c(-1,1), n, replace=TRUE, prob=probs)
     mult <- ifelse(signs < 0, signs*alpha, signs*beta)
     x <- mult*y + mu
     return(x)
     }

###########################################################################
# Stick-Breaking Prior Distribution                                       #
###########################################################################

dStick <- function(theta, gamma, log=FALSE)
     {
     dens <- sum(dbeta(theta, 1, gamma, log=log))
     return(dens)
     }
rStick <- function(M, gamma)
     {
     return(Stick(rbeta(M,1,gamma)))
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
     if(length(nu) > 1) p <- ifelse(nu > 1000000,
          pnorm(q, mu, sigma, lower.tail=lower.tail, log.p=log.p),
          pt({q-mu}/sigma, df=nu, lower.tail=lower.tail, log.p=log.p))
     else p <- if(nu > 1000000) {pnorm(q, mu, sigma,
          lower.tail=lower.tail, log.p=log.p)}
          else {pt({q-mu}/sigma, df=nu, lower.tail=lower.tail,
               log.p=log.p)}
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
     if(length(nu) > 1) q <- ifelse(nu > 1000000,
          qnorm(p, mu, sigma, lower.tail=lower.tail, log.p=log.p),
          mu + sigma * qt(p, df=nu, lower.tail=lower.tail))
     else q <- if(nu > 1000000) {qnorm(p, mu, sigma,
          lower.tail=lower.tail, log.p=log.p)}
          else {mu + sigma * qt(p, df=nu, lower.tail=lower.tail)}
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
     dens <- ifelse(nu > 1000000,
          dnorm(x, mu, sqrt(1/tau), log=TRUE),
          (lgamma((nu+1)/2) - lgamma(nu/2)) + 0.5*log(tau/(nu*pi)) +
          (-(nu+1)/2)*log(1 + (tau/nu)*(x-mu)^2))
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
          (nu/2) * log(det(S)) + ((nu - k - 1)/2) * log(det(Omega)) -
          (tr(as.inverse(S) %*% Omega)/2)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rwishart <- function(nu, S)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(nu < nrow(S)) {
          stop("The nu parameter is less than the dimension of S.")}
     k <- nrow(S)
     Z <- matrix(0, k, k)
     x <- rchisq(k, nu:{nu - k + 1})
     x <- ifelse(x == 0, 1e-100, x)
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
          (nu/2) * log(det(S)) + ((nu - k - 1)/2) * log(det(Omega)) -
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
     x <- ifelse(x == 0, 1e-100, x)
     diag(Z) <- sqrt(x)
     if(k > 1) {
          kseq <- 1:(k-1)
          Z[rep(k*kseq, kseq) +
               unlist(lapply(kseq, seq))] <- rnorm(k*{k - 1}/2)}
     return(chol(crossprod(Z %*% chol(S))))
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
     dens <- dmvn(beta, rep(0, length(beta)),
          g*sigma*sigma*as.inverse(t(X) %*% X), log=log)
     return(dens)
     }

#End
