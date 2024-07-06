## Modelling excess zeros in count data: A new perspective on modelling approaches
## John Haslett, Andrew Parnell, John Hinde, Rafael Moral
## Background functions to reproduce figures on the paper

library(tidyverse)
library(gamlss)

## log-likelihood when y = 0
llik_zero <- Vectorize(function(mu, q, phi, ZI_type, family) {
  pi0 <- switch(family,
                "poisson" = dpois(0, lambda = mu),
                "NBlin" = dNBII(0, mu = mu, sigma = phi),
                "NBquad" = dNBI(0, mu = mu, sigma = phi))
  lz <- switch(ZI_type,
               "A" = log(q),
               "B" = log(pi0^q),
               "C" = log(q + (1 - q) * pi0),
               "D" = log(q * pi0/(1 + (q - 1) * pi0)))
  return(lz)
}, c("mu","q","phi"))

## log-likelihood when y > 0
llik_positive <- Vectorize(function(mu, q, phi, y, ZI_type, family) {
  pi0 <- switch(family,
                "poisson" = dpois(0, lambda = mu),
                "NBlin" = dNBII(0, mu = mu, sigma = phi),
                "NBquad" = dNBI(0, mu = mu, sigma = phi))
  log_piy <- switch(family,
                    "poisson" = dpois(y, lambda = mu, log = TRUE),
                    "NBlin" = dNBII(y, mu = mu, sigma = phi, log = TRUE),
                    "NBquad" = dNBI(y, mu = mu, sigma = phi, log = TRUE))
  rho <- (1 - exp(llik_zero(mu = mu, q = q, phi = phi, ZI_type = ZI_type, family = family)))/(1 - pi0)
  return(log(rho) + log_piy)
}, c("mu","q","phi","y"))

## probability function
dZI <- Vectorize(function(x, mu, q, phi, ZI_type, family, log = FALSE) {
  if(family == "poisson") {
    phi <- NA
  }
  if(x == 0) {
    px <- llik_zero(mu = mu, q = q, phi = phi, ZI_type = ZI_type, family = family)
  } else {
    px <- llik_positive(mu = mu, q = q, phi = phi, y = x, ZI_type = ZI_type, family = family)
  }
  if(!log) {
    px <- exp(px)
  }
  return(px)
}, c("x", "mu", "q", "phi"))

## combining into a single log-likelihood function
llik <- function(theta, y, X_mu, X_zero, ZI_type, family) {
  invlink <- switch(ZI_type,
                    "A" = plogis,
                    "B" = exp,
                    "C" = plogis,
                    "D" = exp)
  np_mu <- ncol(X_mu)
  np_zero <- ncol(X_zero)
  theta1 <- theta[1:np_mu]
  theta2 <- theta[(np_mu + 1):(np_mu + np_zero)]
  if(family != "poisson") {
    phi <- exp(theta[length(theta)])
  } else {
    phi <- 0
  }
  mu <- exp(X_mu %*% theta1) %>% as.numeric
  q <- invlink(X_zero %*% theta2) %>% as.numeric
  loglik_zero <- sum(llik_zero(mu = mu[y == 0], q = q[y == 0], phi = phi, ZI_type = ZI_type, family = family))
  loglik_positive <- sum(llik_positive(mu = mu[y > 0], q = q[y > 0], phi = phi, y = y[y > 0], ZI_type = ZI_type, family = family))
  return(- loglik_zero - loglik_positive)
}

## fitting wrapper
ZI <- function(fmla, ZI_fmla = ~ 1, data,
               family = c("poisson","NBlin","NBquad"), ZI_type = c("A","B","C","D"),
               init, method = "BFGS") {
  X_mu <- model.matrix(fmla, data)
  X_zero <- model.matrix(ZI_fmla, data)
  y <- model.response(model.frame(fmla, data))
  nphi <- switch(family,
                 "poisson" = 0,
                 "NBlin" = 1,
                 "NBquad" = 1)
  name_phi <- switch(family,
                     "poisson" = NULL,
                     "NBlin" = "phi",
                     "NBquad" = "phi")
  if(missing(init)) {
    init <- rep(0, ncol(X_mu) + ncol(X_zero) + nphi)
  }
  fit <- optim(par = init,
               fn = llik,
               y = y,
               X_mu = X_mu,
               X_zero = X_zero,
               ZI_type = ZI_type,
               family = family,
               method = method,
               hessian = TRUE)
  fit_coef <- fit$par
  fit_se <- sqrt(diag(solve(fit$hessian)))
  AIC <- 2 * fit$value + 2 * length(fit_coef)
  model_fit <- list(coef = fit_coef,
                    se = fit_se,
                    X_mu = X_mu,
                    X_zero = X_zero,
                    mu_fmla = fmla,
                    ZI_fmla = ZI_fmla,
                    llik = - fit$value,
                    coef_names = c(colnames(X_mu),colnames(X_zero),name_phi),
                    n = length(y),
                    npar = length(fit_coef),
                    npar_mu = ncol(X_mu),
                    npar_zero = ncol(X_zero),
                    ZI_type = ZI_type,
                    family = family,
                    AIC = AIC)
  class(model_fit) <- "ZI"
  return(model_fit)
}

## print method
print.ZI <- function(obj, ...) {
  family <- obj$family
  coefs <- obj$coef
  npar <- obj$npar
  npar_mu <- obj$npar_mu
  npar_zero <- obj$npar_zero
  n <- obj$n
  ses <- obj$se
  t_val <- coefs/ses
  llik <- obj$llik
  coef_names <- obj$coef_names
  AIC <- obj$AIC
  name_phi <- switch(family,
                     "poisson" = NULL,
                     "NBlin" = "phi",
                     "NBquad" = "phi")
  answer <- data.frame(parameter = c(rep("mu",npar_mu),rep("gamma",npar_zero),name_phi),
                       coefficient = coef_names,
                       estimate = round(coefs, 4),
                       std_error = round(ses, 4),
                       t_value = round(t_val, 4),
                       prob_wald = round(pt(abs(t_val), n - npar, lower.tail = FALSE), 4))
  cat("Family: ", family, "; ZI type ", obj$ZI_type, "\n", sep = "")
  print(answer, row.names = FALSE)
  cat("AIC =", AIC, "\n")
  return(invisible(answer))
}

## AIC method
AIC.ZI <- function(obj, ...) {
  obj$AIC
}

## model.matrix method
model.matrix.ZI <- function(obj, what = c("mu","zero"), ...) {
   X <- switch(what,
               "mu" = obj$X_mu,
               "zero" = obj$X_zero)
   return(X)
}

coef.ZI <- function(obj, what = c("mu","zero"), ...) {
  par_index <- switch(what,
                      "mu" = 1:obj$npar_mu,
                      "zero" = (obj$npar_mu + 1):(obj$npar_mu + obj$npar_zero))
  return(obj$coef[par_index])
}


## predict method
predict.ZI <- function(obj, newdata, what = c("mu","zero"), type = c("link","response"), ...) {
  if(missing(newdata)) {
    X <- model.matrix(obj = obj, what = what)
  } else {
    fmla <- switch(what,
                   "mu" = obj$mu_fmla[-2],
                   "zero" = obj$ZI_fmla)
    X <- model.matrix(fmla, data = newdata)
  }
  beta <- coef(obj = obj, what = what)
  lin_pred <- X %*% beta %>% as.numeric
  g <- exp
  if(what == "zero") {
    g <- switch(obj$ZI_type,
                "A" = plogis,
                "B" = exp,
                "C" = plogis,
                "D" = exp)
  }
  ret <- switch(type,
                "link" = lin_pred,
                "response" = g(lin_pred))
  return(ret)
}

## obtaining means (only for family Poisson)
get_mean <- function(mu, q, ZI_type) {
  rho <- (1 - dZI(x = 0, mu = mu, q = q, ZI_type = ZI_type, family = "poisson"))/
    (1 - dPO(x = 0, mu = mu))
  the_mean <- rho * mu
  return(the_mean)
}