pan.lag <- function(x, ind, lag = 1) {
  unlist(tapply(x,ind, function(x) c(rep(NA, times = lag),x[-((length(x) - lag +1):length(x))])))
}


cdesens <- function(y_mod, m_mod, rho =  seq(-0.9,0.9, by = 0.05), data) {
  
  acde.sens <- rep(NA, times = length(rho))
  acde.sens.se <- rep(NA, times = length(rho))
  res.y <- residuals(lm(onset ~ pan.lag(warl, ccode, 1) +  pan.lag(gdpenl, ccode, 1) + pan.lag(lpop, ccode, 2) + pan.lag(polity2l, ccode, 1) + lmtnest + ncontig + Oil  + ethfrac + relfrac, data = data, subset = !is.na(instab2)))
  res.m <- residuals(lm(instab2 ~  pan.lag(warl, ccode, 1) +  pan.lag(gdpenl, ccode, 1) + pan.lag(lpop, ccode, 2) + pan.lag(polity2l, ccode, 1) + lmtnest + ncontig + Oil  + ethfrac + relfrac, data = data))
  rho.tilde <- cor(res.y, res.m)
  m.fixed <- coef(first)["instab2"] - rho*sd(res.y)*sqrt((1-rho.tilde^2)/(1-rho^2))/sd(res.m)
  n <- nrow(model.matrix(direct))
  
  for (i in 1:length(rho)) {
    sens.direct  <- lm(I(onset - m.fixed[i]*instab2) ~ lmtnest + ncontig  + Oil  + ethfrac + relfrac, data = data, subset = rownames(data) %in% rownames(model.matrix(first)))
    acde.sens[i] <- coef(sens.direct)["ethfrac"]
    sens.vcov <- sequential_g(first, sens.direct, "instab2")
    acde.sens.se[i] <- sqrt(sens.vcov["ethfrac", "ethfrac"])
  }
  
  list(
    direct = sens.direct,
    acde = acde.sens,
    se = acde.sens.se,
    vcov = sens.vcov,
  )
}
