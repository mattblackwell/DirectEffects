context("test_logit")

data("jobcorps")
jc <- jobcorps[
  1:100,
  c("treat", "female", "age_cat", "work2year2q", "pemplq4", "emplq4", "exhealth30")
]

test_that(
  "outreg predicted values", {

    test_outreg <- cde_reg_impute()    
    test_outreg <- set_treatment(test_outreg, treat, ~ female + age_cat)
    test_outreg <- outreg_model(test_outreg, engine = "lm")
    test_outreg <- set_treatment(test_outreg, work2year2q, ~ pemplq4 + emplq4)
    test_outreg <- outreg_model(test_outreg, engine = "lm")

    de_outreg_nocross <- estimate(test_outreg, exhealth30 ~ treat + work2year2q, data = jc, crossfit = FALSE)

    or_m_11_mod <- lm(exhealth30 ~ pemplq4 + emplq4 + female + age_cat,
                      data = jc,
                      subset = treat == 1 & work2year2q == 1)
    pred_m_11 <- predict(or_m_11_mod, newdata = jc)
    expect_equivalent(pred_m_11, de_outreg_nocross$model_fits[[2L]]$outreg_pred[, "1_1"])

    or_m_01_mod <- lm(exhealth30 ~ pemplq4 + emplq4 + female + age_cat,
                      data = jc,
                      subset = treat == 0 & work2year2q == 1)
    pred_m_01 <- predict(or_m_01_mod, newdata = jc)
    expect_equivalent(pred_m_01, de_outreg_nocross$model_fits[[2L]]$outreg_pred[, "0_1"])

    or_m_10_mod <- lm(exhealth30 ~ pemplq4 + emplq4 + female + age_cat,
                      data = jc,
                      subset = treat == 1 & work2year2q == 0)
    pred_m_10 <- predict(or_m_10_mod, newdata = jc)
    expect_equivalent(pred_m_10, de_outreg_nocross$model_fits[[2L]]$outreg_pred[, "1_0"])

    or_m_00_mod <- lm(exhealth30 ~ pemplq4 + emplq4 + female + age_cat,
                      data = jc,
                      subset = treat == 0 & work2year2q == 0)
    pred_m_00 <- predict(or_m_00_mod, newdata = jc)
    expect_equivalent(pred_m_00, de_outreg_nocross$model_fits[[2L]]$outreg_pred[, "0_0"])
    

    or_a_11_mod <- lm(pred_m_11 ~ female + age_cat,
                      data = jc,
                      subset = treat == 1)
    pred_a_11 <- predict(or_a_11_mod, newdata = jc)
    expect_equivalent(pred_a_11, de_outreg_nocross$model_fits[[1L]]$outreg_pred[, "1_1"])

    or_a_01_mod <- lm(pred_m_01 ~ female + age_cat,
                      data = jc,
                      subset = treat == 0)
    pred_a_01 <- predict(or_a_01_mod, newdata = jc)
    expect_equivalent(pred_a_01, de_outreg_nocross$model_fits[[1L]]$outreg_pred[, "0_1"])

    or_a_10_mod <- lm(pred_m_10 ~ female + age_cat,
                      data = jc,
                      subset = treat == 1)
    pred_a_10 <- predict(or_a_10_mod, newdata = jc)
    expect_equivalent(pred_a_10, de_outreg_nocross$model_fits[[1L]]$outreg_pred[, "1_0"])

    or_a_00_mod <- lm(pred_m_00 ~ female + age_cat,
                      data = jc,
                      subset = treat == 0)
    pred_a_00 <- predict(or_a_00_mod, newdata = jc)
    expect_equivalent(pred_a_00, de_outreg_nocross$model_fits[[1L]]$outreg_pred[, "0_0"])

})
