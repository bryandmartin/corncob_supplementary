## @knitr methods

library(corncob)

B <- 1000

checkSep <- function(mod, mod_null) {
  sep1 <- brglm2::detect_separation(y = cbind(mod$dat$W, mod$dat$M - mod$dat$W),
                                    x = mod$X.mu, 
                                    family = binomial("logit"))$separation
  return(sep1)
}

mywaldchisq <- function(mod, mod_null = NULL, restrictions = NULL) {
  if (checkSep(mod = mod, mod_null = mod_null)) {
    return(1)
  } else {
    if (is.null(restrictions)) {
      restrictions <- corncob:::getRestrictionTerms(mod = mod,
                                                    mod_null = mod_null)
    }
    chi.val <- try(corncob:::waldchisq_test(mod, restrictions = restrictions$mu,
                                           restrictions.phi = restrictions$phi), 
                   silent = TRUE)
    if (class(chi.val) == "try-error") {
      return(NA)
    }
    dof.dif <- attr(chi.val, "df")
    return(stats::pchisq(chi.val, dof.dif, lower.tail = FALSE))
  }
}

mylrtest <- function(mod, mod_null) {
  dof.dif <- mod$df.model - mod_null$df.model
  chi.val <- 2 * abs(mod$logL - mod_null$logL)
  chi.val <- chi.val * (abs(chi.val) >= sqrt(.Machine$double.eps))
  pvalue <- stats::pchisq(chi.val, dof.dif, lower.tail = FALSE)
  return(pvalue)
}

mydoBoot <- function(mod, mod_null, test) {
  newW <- simulate(object = mod_null, nsim = nrow(mod_null$dat))
  newdf <- mod$dat
  newdf$W <- newW
  newout_null <- try(corncob::bbdml(formula = mod_null$formula,
                                    phi.formula = mod_null$phi.formula, 
                                    link = mod_null$link,
                                    phi.link = mod_null$phi.link, 
                                    data = newdf), silent = TRUE)
  newout_alt <- try(corncob::bbdml(formula = mod$formula,
                                   phi.formula = mod$phi.formula, 
                                   link = mod$link, phi.link = mod$phi.link,
                                   data = newdf), 
                    silent = TRUE)
  if ("try-error" %in% c(class(newout_null), class(newout_alt))) {
    return(NA)
  } else {
    if (test == "LRT") {
      test.stat <- 2 * abs(newout_alt$logL - newout_null$logL)
      test.stat <- test.stat * (abs(test.stat) >= sqrt(.Machine$double.eps))
    }
    else if (test == "Wald") {
      if (checkSep(mod = newout_alt, mod_null = newout_null)) {
        return(1)
      } 
      restrictions <- corncob:::getRestrictionTerms(mod = newout_alt, 
                                                    mod_null = newout_null)
      test.stat <- try(corncob:::waldchisq_test(mod = newout_alt,
                                                restrictions = restrictions$mu,
                                          restrictions.phi = restrictions$phi),
                       silent = TRUE)
      if (class(test.stat) == "try-error") {
        return(NA)
      }
    }
    return(test.stat)
  }
}


mypbWald <- function(mod, mod_null, B = 100) {
  if (checkSep(mod = mod, mod_null = mod_null)) {
    return(1)
  } else {
    restrictions <- corncob:::getRestrictionTerms(mod = mod,
                                                  mod_null = mod_null)
    t.observed <- try(corncob:::waldchisq_test(mod,
                                               restrictions = restrictions$mu,
                                          restrictions.phi = restrictions$phi), 
                      silent = TRUE)
    if (class(t.observed) == "try-error") {
      return(NA)
    }
    BOOT <- rep(NA, B)
    for (j in 1:B) {
      BOOT[j] <- mydoBoot(mod = mod, mod_null = mod_null, test = "Wald")
    }
    perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= 
                                           x))/(length(stats::na.omit(y)) + 1)
    p.val <- perc.rank(t.observed, BOOT)
    return(p.val)
  }
}

mypbLRT <- function(mod, mod_null, B = 100) {
  t.observed <- 2 * abs(mod$logL - mod_null$logL)
  t.observed <- t.observed * (abs(t.observed) >= sqrt(.Machine$double.eps))
  BOOT <- rep(NA, B)
  for (j in 1:B) {
    BOOT[j] <- mydoBoot(mod = mod, mod_null = mod_null, test = "LRT")
  }
  perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= 
                                         x))/(length(stats::na.omit(y)) + 1)
  p.val <- perc.rank(t.observed, BOOT)
  return(p.val)
}


# Type I error - Setting 1
setting1 <- new_method("corn1", "Corncob Setting 1",
                       method = function(model, draw) {
                         df <- data.frame(W = draw,
                                          M = model$M,
                                          X = model$X,
                                          Xstar = model$Xstar)
                         mod_null <- bbdml(formula = cbind(W, M-W)~1,
                                           phi.formula = ~1,
                                           link = "logit",
                                           phi.link = "logit",
                                           data = df,
                                           inits = rbind(c(-5.75,-5.24)))
                         mod <- bbdml(formula = cbind(W, M-W)~X1,
                                      phi.formula = ~Xst1,
                                      link = "logit",
                                      phi.link = "logit",
                                      data = df,
                                      inits = rbind(c(-5.75,0,-5.24,0)))
                         
                         waldp <- mywaldchisq(mod = mod, mod_null = mod_null)
                         lrtp <- mylrtest(mod = mod, mod_null = mod_null)
                         pbwaldp <- mypbWald(mod = mod, mod_null = mod_null,
                                             B = B)
                         pblrtp <- mypbLRT(mod = mod, mod_null = mod_null,
                                           B = B)
                         
                         list(mod = mod, mod_null = mod_null, 
                              wald = waldp, lrt = lrtp,
                              pbwald = pbwaldp, pblrt = pblrtp)
                       })

# Setting 2
setting2 <- new_method("corn2", "Corncob Setting 2",
                       method = function(model, draw) {
                         df <- data.frame(W = draw,
                                          M = model$M,
                                          X = model$X,
                                          Xstar = model$Xstar)
                         mod_null <- bbdml(formula = cbind(W, M-W)~X1,
                                           phi.formula = ~1,
                                           link = "logit",
                                           phi.link = "logit",
                                           data = df,
                                           inits = rbind(c(-5.36,-1.12,-5.69)))
                         mod <- bbdml(formula = cbind(W, M-W)~X1,
                                      phi.formula = ~Xst1,
                                      link = "logit",
                                      phi.link = "logit",
                                      data = df,
                                      inits = rbind(c(-5.36,-1.12,-5.69,0)))
                         
                         waldp <- mywaldchisq(mod = mod, mod_null = mod_null)
                         lrtp <- mylrtest(mod = mod, mod_null = mod_null)
                         pbwaldp <- mypbWald(mod = mod, mod_null = mod_null,
                                             B = B)
                         pblrtp <- mypbLRT(mod = mod, mod_null = mod_null,
                                           B = B)
                         
                         list(mod = mod, mod_null = mod_null, 
                              wald = waldp, lrt = lrtp,
                              pbwald = pbwaldp, pblrt = pblrtp)
                       })

# Setting 3
setting3 <- new_method("corn3", "Corncob Setting 3",
                       method = function(model, draw) {
                         df <- data.frame(W = draw,
                                          M = model$M,
                                          X = model$X,
                                          Xstar = model$Xstar)
                         mod_null <- bbdml(formula = cbind(W, M-W)~1,
                                           phi.formula = ~Xst1,
                                           link = "logit",
                                           phi.link = "logit",
                                           data = df,
                                           inits = rbind(c(-5.51,-5.38,-.70)))
                         mod <- bbdml(formula = cbind(W, M-W)~X1,
                                      phi.formula = ~Xst1,
                                      link = "logit",
                                      phi.link = "logit",
                                      data = df,
                                      inits = rbind(c(-5.51,0,-5.38,-.70)))
                         
                         waldp <- mywaldchisq(mod = mod, mod_null = mod_null)
                         lrtp <- mylrtest(mod = mod, mod_null = mod_null)
                         pbwaldp <- mypbWald(mod = mod, mod_null = mod_null,
                                             B = B)
                         pblrtp <- mypbLRT(mod = mod, mod_null = mod_null,
                                           B = B)
                         
                         list(mod = mod, mod_null = mod_null, 
                              wald = waldp, lrt = lrtp,
                              pbwald = pbwaldp, pblrt = pblrtp)
                       })
