context("test-subset")

test_that("automatic subsetting of second stage model", {
  
  data(ploughs)
  ploughs$centered_ln_inc <- ploughs$ln_income - mean(ploughs$ln_income, na.rm = TRUE)
  ploughs$centered_ln_incsq <- ploughs$centered_ln_inc^2
  
  fit_first <- lm(
    women_politics ~ plow + centered_ln_inc + centered_ln_incsq + 
      agricultural_suitability + tropical_climate +  large_animals +
      political_hierarchies + economic_complexity +
      rugged + years_civil_conflict +
      years_interstate_conflict  + oil_pc +
      european_descent + communist_dummy + polity2_2000 +
      serv_va_gdp2000, data = ploughs)
  
  form_main <- women_politics ~ plow + 
    agricultural_suitability + tropical_climate + large_animals + 
    political_hierarchies + economic_complexity + 
    rugged | centered_ln_inc +centered_ln_incsq
  
  # different types of subsetting
  ## default
  d0 <- sequential_g(formula = form_main, first_mod = fit_first, data = ploughs)
  
  ## specify
  d1 <- sequential_g(formula = form_main, first_mod = fit_first, data = ploughs, 
                     subset = rownames(ploughs) %in% rownames(model.matrix(fit_first)))
  
  ## custom subset by using symbols
  f2 <- update(fit_first, subset = tropical_climate > 0.5)
  d2 <- sequential_g(formula = form_main, first_mod = f2, data = ploughs,
                     subset = tropical_climate > 0.5)

  # test that rownames are the same  
  expect_equal(nrow(d0$model), nrow(model.matrix(fit_first)))
  expect_equal(nrow(d1$model), nrow(model.matrix(fit_first)))
  expect_equal(nrow(d2$model), nrow(model.matrix(f2)))
})


