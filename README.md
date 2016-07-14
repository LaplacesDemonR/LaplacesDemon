LaplacesDemon
=============

A complete environment for Bayesian inference within R

The goal of `LaplacesDemon`, often referred to as LD, is to provide a complete and self-contained Bayesian environment within R. For example, this package includes dozens of MCMC algorithms, Laplace Approximation, iterative quadrature, Variational Bayes, parallelization, big data, PMC, over 100 examples in the Examples vignette, dozens of additional probability distributions, numerous MCMC diagnostics, Bayes factors, posterior predictive checks, a variety of plots, elicitation, parameter and variable importance, Bayesian forms of test statistics (such as Durbin-Watson, Jarque-Bera, etc.), validation, and numerous additional utility functions, such as functions for multimodality, matrices, or timing your model specification. Other vignettes include an introduction to Bayesian inference, as well as a tutorial.

There are many plans for the growth of this package, and many are long-term plans such as to cotinuously stockpile distributions, examples, samplers, and optimization algorithms. Contributions to this package are welcome.

The main function in this package is the `LaplacesDemon` function, and the best place to start is probably with the LaplacesDemon Tutorial vignette.

# Installation #
---

Using the 'devtools' package:

    install.packages("devtools")
    library(devtools)
    install_github("LaplacesDemonR/LaplacesDemon")


Important Note
=============

`LaplacesDemon` was initially developed and uploaded to CRAN by Byron Hall, most likely the owner of Statisticat, LLC. Later on, the maintainer of the package (or the name of the maintainer) changed to Martina Hall. 

The last version available on CRAN from the original authors and maintainers (i.e., Byron or Martina Hall) was version 13.03.04, which was removed from CRAN on 2013-07-16 at the request of the maintainer. 

After removal from CRAN, the development of `LaplacesDemon` continued for some time on GitHub, most likely still by the original author(s) which now only went by their company name Statisticat, LLC. The last commit by Statisticat for `LaplacesDemon` on GitHub was performed on 25. Mar 2015. After that Statisticat deleted their account on GitHub and ceased further development of the package. 

As Statisticat could not be reached, neither by e-mail nor by snail-mail (the latter was attempted by Rasmus Bååth), Henrik Singmann took over as maintainer of `LaplacesDemon` in July 2016 with the goal to resubmit the package to CRAN. 

To contribute to the development of `LaplacesDemon` or discuss the development please visit its new repository: https://github.com/LaplacesDemonR/LaplacesDemon
