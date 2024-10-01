
# DirectEffects <a href='https://mattblackwell.github.io/DirectEffects'><img src='man/figures/logo.png' align="right" height="138" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/mattblackwell/DirectEffects/workflows/R-CMD-check/badge.svg)](https://github.com/mattblackwell/DirectEffects/actions)
<!-- badges: end -->

## Overview

DirectEffects is an R package to estimate controlled direct effects
(CDEs), which are the effect of a treatment fixing a set of downstream
mediators to particular values. As of now, the package supports general
doubly robust estimation of CDEs under both selection-on-observables and
difference-in-differences assumptions. sequential g-estimation and a
two-stage matching approach called telescope matching. For more
information on how CDEs can be useful for applied research and a brief
introduction to sequential g-estimation, see [Acharya, Blackwell, and
Sen
(2016)](http://www.mattblackwell.org/files/papers/direct-effects.pdf).
For more on the telescope matching procedure, see [Blackwell and
Strezhnev
(2022)](https://www.mattblackwell.org/files/papers/telescope_matching.pdf).
For more on the difference-in-differences approach, see [Blackwell,
Glynn, Hilbig, and Phillips
(2024)](https://www.mattblackwell.org/files/papers/did_cde.pdf).

## Installation

You can install DirectEffects via CRAN for the current stable version or
via GitHub for the development version.

``` r
# Installing from CRAN
install.packages("DirectEffects")

# Installing development version from Github:
# install.packages("devtools")
devtools::install_github("mattblackwell/DirectEffects", build_vignettes = TRUE)
```

## Usage

DirectEffects uses a modular approach to specifying the models and
estimators used to estimate the average controlled direct effect (ACDE).
You can specify a propensity score and outcome regression model for each
treatment variable (that is, each causal variable of interest), along
with a set of parametric and machine learning estimators for these
nuisance functions. The package has a workflow for [models that assume
the
selection-on-observables](https://mattblackwell.github.io/DirectEffects/articles/DirectEffects.html)
and for models that utilize a [difference-in-differences
approach](https://mattblackwell.github.io/DirectEffects/articles/did_cde.html).

``` r
library(DirectEffects)

data(jobcorps)
my_aipw <- cde_aipw() |>
  set_treatment(treat, ~ female + age_cat) |>
  treat_model(engine = "logit") |>
  outreg_model(engine = "lm") |>
  set_treatment(work2year2q, ~ emplq4 + pemplq4) |>
  treat_model(engine = "logit") |>
  outreg_model(engine = "lm") |>
  estimate(exhealth30 ~ treat + work2year2q, data = jobcorps)

broom::tidy(my_aipw)
```

    ##                                            term      estimate  std.error
    ## treat_1_0             treat [(1, 0) vs. (0, 0)]  0.0324855302 0.02248379
    ## treat_1_1             treat [(1, 1) vs. (0, 1)]  0.0301271316 0.01394355
    ## work2year2q_0_1 work2year2q [(0, 1) vs. (0, 0)]  0.0038019794 0.02057364
    ## work2year2q_1_1 work2year2q [(1, 1) vs. (1, 0)] -0.0006775739 0.01663843
    ##                   statistic    p.value     conf.low  conf.high   df
    ## treat_1_0        1.44484233 0.14858424 -0.011595861 0.07656692 3818
    ## treat_1_1        2.16064973 0.03076055  0.002792942 0.05746132 6207
    ## work2year2q_0_1  0.18479854 0.85339644 -0.036533854 0.04413781 3991
    ## work2year2q_1_1 -0.04072342 0.96751774 -0.033294846 0.03193970 6034

## Specific estimator implementations

DirectEffects also implements a number of specific estimators outside
the modular framework. These are often specific estimators proposed in
the literature:

- [`sequential_g()`](https://mattblackwell.github.io/DirectEffects/articles/sequential_g.html):
  estimate controlled direct effects using two-stage linear models.
- [`telescope_match()`](https://mattblackwell.github.io/DirectEffects/articles/telescope_matching.html):
  estimated controlled direct effects using a two-stage matching
  procedure with bias correction.

DirectEffects also provides diagnostics for these two approaches,
including sensitivity analyses and balance checks.
