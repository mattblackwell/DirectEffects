## ----loadpkg, echo = F, include = F--------------------------------------
library(dplyr)
library(ggrepel)
library(scales)
library(DirectEffects) # devtools::install_github("mattblackwell/DirectEffects", ref = "develop")

## ---- echo = TRUE, eval = F----------------------------------------------
#  install_github("mattblackwell/DirectEffects", ref = "master")
#  library(DirectEffects)

## ---- echo = F, out.width = "600px"--------------------------------------
knitr::include_graphics("figures/ABS_fig3.png")

## ------------------------------------------------------------------------
data(ploughs)

## ---- echo = FALSE, eval = TRUE------------------------------------------
tbl_df(ploughs)

## ---- echo = F, fig.width= 7,fig.height=4--------------------------------
ggplot(ploughs, aes(x = plow, women_politics/100, size = exp(ln_income))) + 
  labs(
    title = "Bivariate relationship between agricultural practice and modern gender equity",
    caption = "Source: Alesina, Giuliano, and Nunn (2013)",
    y = "% political positions held by women (2000)",
    x = "% ethnic groups that traditionally used the plough"
  ) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(label = percent) +
  scale_y_continuous(label = percent) +
  scale_size_continuous(name = "Country Income (USD)") +
  theme_light()

## ------------------------------------------------------------------------
ate.mod <- lm(women_politics ~ plow + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged, data = ploughs)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
summary(ate.mod)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
ploughs$centered_ln_inc <- ploughs$ln_income - mean(ploughs$ln_income, na.rm = TRUE)
ploughs$centered_ln_incsq <- ploughs$centered_ln_inc^2

## ---- echo = TRUE, eval = TRUE-------------------------------------------
first <- lm(women_politics~ plow + centered_ln_inc + centered_ln_incsq + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged +years_civil_conflict +  years_interstate_conflict  + oil_pc + european_descent + communist_dummy +polity2_2000+serv_va_gdp2000, data = ploughs)

## ------------------------------------------------------------------------
form.main <- women_politics ~ plow + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged | centered_ln_inc + centered_ln_incsq

## ---- echo = TRUE, eval = TRUE-------------------------------------------
direct <- sequential_g(formula = form.main,
                       first_mod = first,
                       data = ploughs,
                       subset = rownames(ploughs) %in% rownames(model.matrix(first)))

## ---- eval = F-----------------------------------------------------------
#  mf <- match.call(expand.dots = FALSE)
#  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
#  mf <- mf[c(1L, m)]

## ---- eval = F-----------------------------------------------------------
#  bt <- terms(formula, data = data, rhs = 2)
#  bnames <- attr(bt, "term.labels")
#  bvars <- match(bnames, names(fcoefs), 0L)
#  mf$formula <- formula

## ---- eval = F-----------------------------------------------------------
#  Y <- model.response(mf, "numeric")
#  M <- model.matrix(bt, mf, contrasts)
#  Y <- Y - M %*% fcoefs[bvars]

## ---- eval = F-----------------------------------------------------------
#  mtX <- terms(formula, data = data, rhs = 1)
#  X <- model.matrix(mtX, mf, contrasts)
#  out <- lm.fit(X, Y, offset = offset, ...)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
summary(direct)

