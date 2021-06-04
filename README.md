[![Travis-CI Build Status](https://travis-ci.org/LaplacesDemonR/LaplacesDemon.svg?branch=master)](https://travis-ci.org/LaplacesDemonR/LaplacesDemon)

LaplacesDemon
=============

A complete environment for Bayesian inference within R

The goal of `LaplacesDemon`, often referred to as LD, is to provide a complete and self-contained Bayesian environment within R. For example, this package includes dozens of MCMC algorithms, Laplace Approximation, iterative quadrature, Variational Bayes, parallelization, big data, PMC, over 100 examples in the Examples vignette, dozens of additional probability distributions, numerous MCMC diagnostics, Bayes factors, posterior predictive checks, a variety of plots, elicitation, parameter and variable importance, Bayesian forms of test statistics (such as Durbin-Watson, Jarque-Bera, etc.), validation, and numerous additional utility functions, such as functions for multimodality, matrices, or timing your model specification. Other vignettes include an introduction to Bayesian inference, as well as a tutorial.

There are many plans for the growth of this package, and many are long-term plans such as to cotinuously stockpile distributions, examples, samplers, and optimization algorithms. Contributions to this package are welcome.

The main function in this package is the `LaplacesDemon` function, and the best place to start is probably with the LaplacesDemon Tutorial vignette.

# Installation #
---

From `CRAN`
    install.packages("LaplacesDemon")


Using the 'devtools' package:

    install.packages("devtools")
    library(devtools)
    install_github("LaplacesDemonR/LaplacesDemon")


Package History
=============

`LaplacesDemon` was initially developed and uploaded to CRAN by Byron Hall, the owner of Statisticat, LLC. Later on, the maintainer of the package changed to Martina Hall. 

The last version available on CRAN from the original authors and maintainers was version 13.03.04, which was removed from CRAN on 2013-07-16 at the request of the maintainer. 

After removal from CRAN, the development of `LaplacesDemon` continued for some time on GitHub under the name of Statisticat LLC (presumably still run by Byron Hall). The last commit by Statisticat for `LaplacesDemon` on GitHub was performed on 25. Mar 2015. After that Statisticat deleted their account on GitHub and ceased further development of the package. 

As Statisticat could not be reached, neither by e-mail nor by snail-mail (the latter was attempted by Rasmus Bååth), Henrik Singmann took over as maintainer of `LaplacesDemon` in July 2016 with the goal to resubmit the package to CRAN (as version 16.0.x). Henrik Singmann does not actively continue the development of `LaplacesDemon` but only retains it on CRAN in its current state.

Note that in order to resubmit the package to CRAN all links to the now defunct website of Statisticat (formerly: http://www.bayesian-inference.com) were replaced with links to versions of this website on the web archive (https://web.archive.org/web/20141224051720/http://www.bayesian-inference.com/index).

To contribute to the development of `LaplacesDemon` or discuss the development please visit its new repository: https://github.com/LaplacesDemonR/LaplacesDemon

