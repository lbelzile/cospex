
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cospex: Conditional Spatial Extremes

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License](https://img.shields.io/badge/license-GPL(%3E=%203)-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

<!-- badges: end -->

The `cospex` packages proposes S4 classes for defining and simulating
observations from the conditional spatial extremes (and extensions
thereof) of Wadsworth and Tawn. At maturity, it is expected that the
package functionalities will include

-   [ ] `S4` classes for models
-   [x] `S4` classes for scaling functions
-   [ ] predefined scaling functions
-   [ ] maximum likelihood estimation routines
-   [ ] simulation: unconditional (forward sampling)
-   [ ] simulation: conditional
-   [ ] simulation: importance sampling scheme
-   [ ] diagnostic plot: lag independence
-   [ ] diagnostic plot: Kendall’s *τ*
-   [ ] utilities: covariance models
-   [ ] utilities: distance with geometric anisotropy
-   [ ] `S4` methods for summary objects
-   [ ] extension: addition of nugget term
-   [ ] marginal transformation to standardized margins (with
    semiparametric transformation and point mass for zero inflation)
-   [ ] left-censoring
-   [ ] unit tests
-   [ ] vignettes
-   [ ] linear combination basis functions with gaussian weights in
    place of residual Gaussian process (INLA, multiresolution)

## Installation

<!-- You can install the released version of `cospex` from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("cospex") -->
<!-- ``` -->

You can install the development version of \`cospex from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("lbelzile/cospex")
```

<!-- `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
