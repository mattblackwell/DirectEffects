# DirectEffects

`DirectEffects` is an R package to estimate controlled direct effects (CDEs). As of now, the package supports sequential g-estimation and a two-stage matching approach called telescope matching. For more information on how CDEs can be useful for applied research and a brief introduction to sequential g-estimation, see our [2016 APSR][de-paper]. For more on the telescope matching procedure, see our [working paper][tm-paper]. To install the development version of `DirectEffects`, run the following code in R:
```R
require(devtools)
install_github("mattblackwell/DirectEffects", build_vignettes = TRUE)
```

[de-paper]: http://www.mattblackwell.org/files/papers/direct-effects.pdf
[tm-paper]: https://www.mattblackwell.org/files/papers/telescope_matching.pdf
