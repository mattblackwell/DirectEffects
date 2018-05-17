## ----loadpkg, echo = F, include = F--------------------------------------
library(dplyr)
library(scales)
library(reshape2)
library(ggplot2)
library(DirectEffects) # devtools::install_github("mattblackwell/DirectEffects", ref = "develop")

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  install_github("mattblackwell/DirectEffects", ref = "master")
#  library(DirectEffects)

## ---- eval = FALSE-------------------------------------------------------
#  data(ploughs)
#  ploughs$centered_ln_inc <- ploughs$ln_income - mean(ploughs$ln_income, na.rm = TRUE)
#  ploughs$centered_ln_incsq <- ploughs$centered_ln_inc^2
#  
#  first <- lm(women_politics~ plow + centered_ln_inc + centered_ln_incsq + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged + years_civil_conflict +  years_interstate_conflict  + oil_pc + european_descent + communist_dummy + polity2_2000 + serv_va_gdp2000,
#              data = ploughs)
#  
#  direct <- sequential_g(formula = women_politics ~ plow + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged | centered_ln_inc + centered_ln_incsq,
#                         first_mod = first,
#                         data = ploughs, subset = rownames(ploughs) %in% rownames(model.matrix(first)))

## ---- echo = FALSE, out.width = "600px"----------------------------------
knitr::include_graphics("figures/ABS_fig3.png")

## ------------------------------------------------------------------------
data(ploughs)

## ------------------------------------------------------------------------
ate.mod <- lm(women_politics ~ plow + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged, data = ploughs)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
summary(ate.mod)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
ploughs$centered_ln_inc <- ploughs$ln_income - mean(ploughs$ln_income, na.rm = TRUE)
ploughs$centered_ln_incsq <- ploughs$centered_ln_inc^2

## ---- echo = TRUE, eval = TRUE-------------------------------------------
fit_first <- lm(women_politics ~ plow + centered_ln_inc + centered_ln_incsq +
                  agricultural_suitability +  tropical_climate +  large_animals +
                  political_hierarchies + economic_complexity + rugged +
                  years_civil_conflict +  years_interstate_conflict  + oil_pc +
                  european_descent + communist_dummy  + polity2_2000 +
                  serv_va_gdp2000,  data = ploughs)

## ------------------------------------------------------------------------
form_main <- women_politics ~ plow + agricultural_suitability +  tropical_climate +  large_animals + political_hierarchies + economic_complexity + rugged | centered_ln_inc + centered_ln_incsq

## ---- echo = TRUE, eval = TRUE-------------------------------------------
direct <- sequential_g(formula = form_main,
                       first_mod = fit_first,
                       data = ploughs,
                       subset = rownames(ploughs) %in% rownames(model.matrix(fit_first)))

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
coef(direct)

## ------------------------------------------------------------------------
confint(direct, "plow")

## ------------------------------------------------------------------------
vcov(direct)

## ------------------------------------------------------------------------
head(direct$model)
head(direct$y)

