## ----loadpkg, echo = F, include = F--------------------------------------
library(dplyr)
library(ggrepel)
library(scales)
library(reshape2)
library(DT)
library(DirectEffects) # devtools::install_github("mattblackwell/DirectEffects", ref = "develop")

## ---- echo = TRUE, eval = F----------------------------------------------
#  install_github("mattblackwell/DirectEffects", ref = "master")
#  library(DirectEffects)

## ---- echo = F, out.width = "600px"--------------------------------------
knitr::include_graphics("figures/ABS_fig3.png")

## ------------------------------------------------------------------------
data(ploughs)

## ------------------------------------------------------------------------
tbl_df(ploughs)

## ---- echo = F, fig.width= 7, fig.height=4, warning=F--------------------
ggplot(ploughs, aes(x = plow, women_politics/100, size = exp(ln_income))) + 
  labs(
    title = "Agricultural practice and modern gender equity",
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
fit.first <- lm(women_politics~ plow + centered_ln_inc + centered_ln_incsq + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged +years_civil_conflict +  years_interstate_conflict  + oil_pc + european_descent + communist_dummy  + polity2_2000 + serv_va_gdp2000, 
                data = ploughs)

## ------------------------------------------------------------------------
form.main <- women_politics ~ plow + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged | centered_ln_inc + centered_ln_incsq

## ---- echo = TRUE, eval = TRUE-------------------------------------------
direct <- sequential_g(formula = form.main,
                       first_mod = fit.first,
                       data = ploughs,
                       subset = rownames(ploughs) %in% rownames(model.matrix(fit.first)))

## ---- echo = F, fig.width= 7, fig.height=4, warning=F--------------------
mnames <- c("centered_ln_inc", "centered_ln_incsq")

curve.first <- function(x) {
  coefs <- coefficients(fit.first)[mnames]
  yhat <- coefs[1]*x + coefs[2]*x^2
  return(yhat)
}

ploughs %>% 
  ggplot(., aes_string(x = mnames[1], y = "women_politics")) + 
  geom_point() +
  stat_function(fun = curve.first, size = 1, color = "indianred", size = 2, se = F) +
  labs(title = "De-mediation function (line) overlayed on bivariate relationship (points)") +
  annotate('text', x = 2, y = 0, 
           label = "gamma(m)== 2.6347851*m + 0.8801396*m^{2}",
           parse = TRUE,
           color = "indianred")  +
  theme_light()

## ---- eval = F-----------------------------------------------------------
#    rawY <- model.response(mf, "numeric")
#    M <- model.matrix(object = bt, data = mf, contrasts)
#    gamma <- M %*% fcoefs[bvars] # fitted values from de-mediation function
#    Y <- rawY - gamma

## ---- echo = F, fig.width= 7, fig.height=4-------------------------------
rawYdat <- ploughs$women_politics
mDat <- select_(ploughs, .dots = mnames) %>% as.matrix()
mCoefs <- coefficients(fit.first)[mnames]
gammaDat <- mDat %*% mCoefs 
Ydat <- rawYdat - as.numeric(gammaDat)

# plot
data_frame(Ydat, rawYdat, id = 1:length(rawYdat)) %>%
  melt(., id = "id") %>% 
  filter(!is.na(value)) %>% 
  ggplot(., aes(y = variable, x = value, group = id)) +
  geom_path(arrow = arrow(angle = 15, length = unit(0.2,"cm"), type = "closed",ends =  "first"),
            color = "darkgray", alpha = 0.4) +
  geom_hline(yintercept = 1:2) +
  geom_point() +
  scale_y_discrete(name = "", labels = c("de-mediated values", "original values")) +
  labs(x = "Dependent Variable: Percent of Women in Politics",
       title = "Implications of de-mediation") +
  theme_light()

## ---- eval = F-----------------------------------------------------------
#    mtX <- terms(formula, data = data, rhs = 1) # Y ~ A + X
#    X <- model.matrix(mtX, mf, contrasts) # numeric matrix
#    out <- lm.fit(X, Y, offset = offset, ...)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
summary(direct)

## ---- echo = F, fig.width= 7, fig.height=4-------------------------------
df.coef <- rbind(summary(ate.mod)$coef["plow", ],
               summary(direct)$coef["plow", ]) %>% 
  as.data.frame()  %>%
  mutate(pos = c("ATE of Plows", "ACDE of Plows, via sequential g-estimation"))

ggplot(df.coef, aes(x = Estimate, y = 0)) +
  facet_wrap(~ pos, ncol = 1) +
  geom_point() +
  geom_segment(data = df.coef, aes(x = Estimate + 1.96*`Std. Error`, xend = Estimate - 1.96*`Std. Error`, y = 0, yend = 0)) +
  scale_y_continuous(labels = NULL, limit = c(-1, 1), breaks = NULL, minor_breaks = NULL) +
  scale_x_continuous(limit = c(-10, 5)) +
  geom_vline(xintercept = 0, color = "red", linetype = 2) +
  geom_label(aes(label = round(Estimate, 1)), y = 0.3) +
  theme_light() +
  labs(title = "Estimated Causal Effects on Percent of Women in Political Office",
      x = "",
      y = "",
      caption = "Lines are 95% confidence intervals")

## ------------------------------------------------------------------------
direct$coefficients

## ------------------------------------------------------------------------
direct$vcov

## ------------------------------------------------------------------------
head(direct$model)
head(direct$y)
head(direct$x) # will be non-NULL if the argument x in sequential_g is set to TRUE

## ------------------------------------------------------------------------
str(direct$terms)

