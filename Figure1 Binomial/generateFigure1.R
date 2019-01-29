library(corncob)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(gridExtra)

data(soil_phylo)
soil <- soil_phylo %>% 
  phyloseq::subset_samples(DayAmdmt %in% c(11)) %>%
  tax_glom("Genus") 

myplot_glm <- function(mod, bs = 30, mytitle = NULL, 
                       xlab = "Sample ID", ylab = "Relative Abundance") {
  data0 <- mod$data
  W <- data0$Wi
  M <- data0$Ni
  mu_est <- predict(mod, type = "response")
  sims <- matrix(NA, nrow = 10000, ncol = nrow(data0))
  
  for (i in 1:10000) {
    sim <- simulate(mod)$sim_1
    refit <- glm(sim ~ 1, family = binomial)
    repred <- predict(refit, type = "response", newdata = data0)
    sims[i,] <- rbinom(length(repred), size = data0$Ni, prob = repred)
  }
  
  predint <- apply(sims, 2, quantile, c(.025, .975))
  ymin <- predint[1,]/data0$Ni
  ymax <- predint[2,]/data0$Ni

  df <- data.frame(RA = W/M,
                   E_RA = mu_est,
                   group = factor(mod$data$DayAmdmt),
                   index = 1:nrow(mod$data),
                   ymin = ymin, ymax = ymax
  )
  ggplot(df, aes(x = index, y = RA, group = group, color = group)) +
    geom_point(aes(size = 2)) + 
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .2) + 
    labs(x = xlab, y = ylab, title = mytitle) +
    theme_bw(base_size = bs) + ylim(0,.008) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
}

plot_bbdml <- function(mod, bs = 30, mytitle = NULL, 
                       xlab = "Sample ID", ylab = "Relative Abundance") {
  data0 <- mod$dat
  W <- data0$Wi
  M <- data0$Ni
  mu_est <- mod$mu.resp
  sims <- matrix(NA, nrow = 10000, ncol = nrow(data0))
  
  for (i in 1:10000) {
    sim <- simulate(mod, nsim = length(W))
    newdat <- data.frame(W = sim, M = M)
    refit <- bbdml(cbind(W, M - W)~1, phi.formula = ~1, data = newdat)
    sims[i,] <- simulate(refit, nsim = length(W))
  }
  
  predint <- apply(sims, 2, quantile, c(.025, .975))
  ymin <- predint[1,]/data0$Ni
  ymax <- predint[2,]/data0$Ni
  
  df <- data.frame(RA = W/M,
                   E_RA = mu_est,
                   group = factor(mod$dat$DayAmdmt),
                   index = 1:nrow(mod$dat),
                   ymin = ymin, ymax = ymax
  )
  ggplot(df, aes(x = index, y = RA, group = group, color = group)) +
    geom_point(aes(size = 2)) + 
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .2) + 
    labs(x = xlab, y = ylab, title = mytitle) +
    theme_bw(base_size = bs) + ylim(0,.008) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
}



Mi <- soil %>% sample_sums %>% c
Wi <- prune_taxa("OTU.82", soil) %>% otu_table %>% c
otu_to_taxonomy(OTU = "OTU.82", data = soil)
dat1 <- data.frame(Wi = Wi, Ni = Mi, DayAmdmt = 11)
binom0 <- glm(cbind(Wi, Ni-Wi)~ 1, family = binomial, data = dat1)
g1 <- myplot_glm(binom0, mytitle = "Binomial Fit")

tmp2 <- bbdml(cbind(Wi, Ni-Wi)~1, phi.formula = ~1, data = dat1)
g2 <- plot_bbdml(tmp2, mytitle = "Beta-Binomial Fit")

pdf(file = 'modelFits.pdf',width = 24, height = 8.5)
gridExtra::grid.arrange(g1, g2, layout_matrix = rbind(c(1,2)))
dev.off()
