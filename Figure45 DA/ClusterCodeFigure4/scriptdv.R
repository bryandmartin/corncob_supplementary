args <- commandArgs(TRUE)


devtools::install_github("bryandmartin/corncob")
require(GenSA)
require(corncob)
require(VGAM)
require(knitr)
require(ggplot2)
require(magrittr)
require(dplyr)
require(doParallel)
require(doSNOW)
require(simulator)
require(brglm2)
require(Rmpi)

mysocks <- rep(c("b07.local", "b08.local", "b09.local", "b10.local",
                 "b13.local", "b14.local", "b11.local", "b16.local",
                 "b17.local", "b18.local"), each = 10)


data(soil_phylo)
soil <- soil_phylo %>% 
  phyloseq::subset_samples(DayAmdmt %in% c(20, 21)) %>%
  phyloseq::tax_glom("Genus")


data <- soil
formula <- phi.formula <- formula_null <- ~ DayAmdmt
phi.formula_null <- ~ 1
link <- phi.link <- "logit"




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

mypbLRT <- function(mod, mod_null, B = 1000000) {
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


doSim <- function(i, data, formula, phi.formula, 
                  formula_null, phi.formula_null, 
                  link = "logit", phi.link = "logit",
                  B = 1000000) {
  data_i <- convert_phylo(data, select = phyloseq::taxa_names(data)[i])
  if (sum(data_i$W) == 0) {
    return(1)
  }
  formula_i <- stats::update(formula, cbind(W, M) ~ .)
  formula_null_i <- stats::update(formula_null, cbind(W, 
                                                      M) ~ .)
  
  fit_unr <-  try(bbdml(formula = formula_i, phi.formula = phi.formula, 
                        data = data_i, link = link, phi.link = phi.link),
                  silent = TRUE)
  fit_red <-  try(bbdml(formula = formula_null_i, 
                        phi.formula = phi.formula_null, data = data_i,
                        link = link, phi.link = phi.link), silent = TRUE)
  
  if ("try-error" %in% c(class(fit_unr), class(fit_red))) {
    return(NA)
  }
  p_val <- mypbLRT(mod = fit_unr, mod_null = fit_red, B = B)
  return(p_val)
}




cl <- parallel::makePSOCKcluster(names = mysocks)
libs <- "corncob"
parallel::clusterExport(cl, varlist = ls(envir = globalenv()))
parallel::clusterEvalQ(cl, sapply(libs,
                                  function(pkgnam) { do.call("library",
                                                             list(pkgnam))}))


corn_p_dv <- parallel::parLapply(cl, seq(length(phyloseq::taxa_names(data))),
                                 fun = doSim, data = data, formula = formula,
                                 phi.formula = phi.formula, 
                                 formula_null = formula_null,
                                 phi.formula_null = phi.formula_null, 
                                 link = "logit", phi.link = "logit",
                                 B = 1000000)
parallel::stopCluster(cl)


save(corn_p_dv, file =  "cornpdv.rda")







