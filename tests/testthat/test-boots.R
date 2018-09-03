
context("test-boots")

test_that("automatic subsetting of second stage model", {
  data(ploughs)
  ploughs$centered_ln_inc <- ploughs$ln_income - mean(ploughs$ln_income, na.rm = TRUE)
  ploughs$centered_ln_incsq <- ploughs$centered_ln_inc^2

  fit_first <- lm(
    women_politics ~ plow + centered_ln_inc + centered_ln_incsq +
      agricultural_suitability + tropical_climate + large_animals +
      political_hierarchies + economic_complexity +
      rugged + years_civil_conflict +
      years_interstate_conflict + oil_pc +
      european_descent + communist_dummy + polity2_2000 +
      serv_va_gdp2000,
    data = ploughs
  )

  form_main <- women_politics ~ plow +
    agricultural_suitability + tropical_climate + large_animals +
    political_hierarchies + economic_complexity +
    rugged | centered_ln_inc + centered_ln_incsq

  # different types of subsetting
  ## default
  d0 <- sequential_g(
    formula = form_main, first_mod = fit_first, data = ploughs,
    bootstrap = "none"
  )

  d0tab <- coef(summary(d0))

  ## specify boot
  d1 <- sequential_g(
    formula = form_main, first_mod = fit_first, data = ploughs,
    bootstrap = "standard", boots_n = 100
  )
  d1tab <- coef(summary(d1))

  # test that rownames are the same
  expect_equal(coef(d1), coef(d0))
  expect_gt(sum(abs(d0tab[, "Std. Error"] - d1tab[, "Std. Error"])), 0)
  expect_gt(sum(abs(d0tab[, "t value"] - d1tab[, "t value"])), 0)
  expect_gt(sum(abs(d0tab[, "Pr(>|t|)"] - d1tab[, "Pr(>|t|)"])), 0)
})
