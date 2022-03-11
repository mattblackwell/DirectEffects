context("test_logit")

data("jobcorps")
jc <- jobcorps[
  1:100,
  c("treat", "female", "age_cat", "work2year2q", "pemplq4", "emplq4", "exhealth30")
]

test_that(
  "logit predicted values", {

    test_ipw <- cde_ipw() |>
      set_treatment(treat, ~ female + age_cat) |>
      treat_model(engine = "logit") |>
      set_treatment(work2year2q, ~ pemplq4 + emplq4) |>
      treat_model(engine = "logit")

    de_ipw_nocross <- test_ipw |> estimate(exhealth30 ~ treat + work2year2q, jc, crossfit = FALSE)

    ipw_a1_mod <- glm(treat ~ female + age_cat, data = jc, family = binomial())
    ipw_a1_preds <- ipw_a1_mod$fitted
    expect_equivalent(ipw_a1_mod$fitted, de_ipw_nocross$model_fits[[1]]$treat_pred[, "1"])

    ipw_a2_mod_0 <- glm(work2year2q ~ female + age_cat + pemplq4 + emplq4,
                        data = jc, subset = treat == 0, family = binomial())
    ipw_a2_preds_0 <- predict(ipw_a2_mod_0, newdata = jc, type = "response")
    ipw_a2_mod_1 <- glm(work2year2q ~ female + age_cat + pemplq4 + emplq4,
                        data = jc, subset = treat == 1, family = binomial())
    ipw_a2_preds_1 <- predict(ipw_a2_mod_1, newdata = jc, type = "response")
    
    expect_equivalent(ipw_a2_preds_0, de_ipw_nocross$model_fits[[2]]$treat_pred[, "0_1"])
    expect_equivalent(ipw_a2_preds_1, de_ipw_nocross$model_fits[[2]]$treat_pred[, "1_1"])

    w_11 <- 1 / (ipw_a1_preds * ipw_a2_preds_1)
    w_10 <- 1 / (ipw_a1_preds * (1 - ipw_a2_preds_1))
    w_01 <- 1 / ((1 - ipw_a1_preds) * ipw_a2_preds_0)
    w_00 <- 1 / ((1 - ipw_a1_preds) * (1 - ipw_a2_preds_0))
    d_11 <- as.numeric(jc$treat == 1 & jc$work2year2q == 1)
    d_10 <- as.numeric(jc$treat == 1 & jc$work2year2q == 0)
    d_01 <- as.numeric(jc$treat == 0 & jc$work2year2q == 1)
    d_00 <- as.numeric(jc$treat == 0 & jc$work2year2q == 0) 
    w <- d_11 * w_11 + d_10 * w_10 + d_01 * w_01 + d_00 * w_00
    
    mu_11 <- sum(w_11 * d_11 * jc$exhealth30) / nrow(jc)
    mu_01 <- sum(w_01 * d_01 * jc$exhealth30) / nrow(jc)
    mu_10 <- sum(w_10 * d_10 * jc$exhealth30) / nrow(jc)
    mu_00 <- sum(w_00 * d_00 * jc$exhealth30) / nrow(jc)
    lm_out <- lm(exhealth30 ~ treat * work2year2q, data = jc, weights = w)
    
    expect_equal(
      unname(lm_out$coefficients[2]),
      de_ipw_nocross$estimates[1, "estimate"]
    )
    expect_equivalent(
      unname(lm_out$coefficients[2] + lm_out$coefficients[4]),
      de_ipw_nocross$estimates[2, "estimate"],
      tolerance = 1e-04
    )    
    expect_equivalent(
      unname(lm_out$coefficients[3]),
      de_ipw_nocross$estimates[3, "estimate"]
    )
    expect_equivalent(
      unname(lm_out$coefficients[3] + lm_out$coefficients[4]),
      de_ipw_nocross$estimates[4, "estimate"],
      tolerance = 1e-04
    )

})
