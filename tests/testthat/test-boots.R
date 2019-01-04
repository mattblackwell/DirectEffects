
context("test-boots")

test_that("automatic subsetting of second stage model", {

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
  rugged | + years_civil_conflict +
  years_interstate_conflict  + oil_pc +
  european_descent + communist_dummy + polity2_2000 +
  serv_va_gdp2000 | centered_ln_inc +centered_ln_incsq


  # different types of subsetting
  ## default
  d0  <- sequential_g(formula = form_main,
                      data = ploughs)

  d0boots <- boots_g(d0, boots = 100)

  # test that rownames are the same
  expect_equal(nrow(d0boots), 100)
  expect_equal(ncol(d0boots), length(d0$coefficients))
})
