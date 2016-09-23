# DirectEffects

`DirectEffects` is an R package to estimate controlled direct effects (CDEs). As of now, the only model supported is sequential g-estimation, but we plan to expand to other models, including doubly robust estimators, in the future. For more information on how CDEs can be useful for applied research and a brief introduction to sequential g-estimation, see our [2016 APSR][de-paper]. To install the development version of `DirectEffects`, run the following code in R:
```R
require(devtools)
install_github("mattblackwell/DirectEffects", ref = "master")
```

[de-paper]: http://www.mattblackwel.org/files/papers/direct-effects.pdf
