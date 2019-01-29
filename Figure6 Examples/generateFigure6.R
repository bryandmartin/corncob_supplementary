require(corncob)
require(ggplot2)
require(phyloseq)
require(magrittr)
plot_bbdml_nolab <- function(x, AA = FALSE, color = NULL, shape = NULL,
                             facet = NULL, title = NULL, ...) {
  mod <- x
  mu_est <- mod$mu.resp
  phi_est <- mod$phi.resp
  M <- mod$M
  W <- mod$W
  
  mod$dat$DayAmdmt <- gsub("20", "None", mod$dat$DayAmdmt)
  mod$dat$DayAmdmt <- gsub("21", "Biochar", mod$dat$DayAmdmt)
  
  data0 <- mod$dat
  mu_est <- mod$mu.resp
  sims <- matrix(NA, nrow = 10000, ncol = nrow(data0))
  
  for (i in 1:10000) {
    sim <- simulate(mod, nsim = length(W))
    newdat <- data.frame(W = sim, M = M, DayAmdmt = mod$dat$DayAmdmt)
    refit <- bbdml(mod$formula, phi.formula = mod$phi.formula, data = newdat)
    sims[i,] <- simulate(refit, nsim = length(W))
  }
  
  predint <- apply(sims, 2, quantile, c(.025, .975))
  ymin <- predint[1,]/M
  ymax <- predint[2,]/M
  resp <- W/M

  samp_names <- names(W)
  dat_noNA <- mod$dat[samp_names, ]
  df <- data.frame(RA = resp, samples = samp_names, ymin = ymin, 
                   ymax = ymax)
  my_ord_str <- ""
  custom_color <- custom_shape <- FALSE
  color_name <- shape_name <- NULL
  if (!is.null(color)) {
    if (length(color) == 1) {
      df[[color]] <- factor(dat_noNA[[color]])
      color_name <- color
      my_ord_str <- paste(my_ord_str, df[[color]], sep = "_")
    }
    else if (length(color) == nrow(df)) {
      df[["color"]] <- color
      color <- color_name <- "color"
      custom_color <- TRUE
    }
    else {
      stop("color must either match a variable or be a custom vector of correct length!")
    }
  }
  my_ord_str <- paste(my_ord_str, df$samples, sep = "_")
  df$order <- factor(df$samples, levels = df$samples[order(my_ord_str)])
  ylab_tmp <- ifelse(!AA, "Relative Abundance", "Absolute Abundance")
  aes_map <- ggplot2::aes_string(x = "order", y = "RA", colour = color, 
                                 shape = shape, labs = "samples")
  my_gg <- ggplot2::ggplot(df, aes_map) + 
    ggplot2::geom_point(size = 5) + 
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), 
                           width = 0.2) + 
    ggplot2::labs(title = title, x = "Sample", y = ylab_tmp,
                  colour = "Addition", shape = shape_name) + 
    ggplot2::theme_bw(base_size = 25) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) + 
    guides(colour = FALSE)
  my_gg
}



data(soil_phylo)
soil <- soil_phylo %>% 
  phyloseq::subset_samples(DayAmdmt %in% c(20, 21)) %>%
  phyloseq::tax_glom("Genus")

i <- 49
data_i <- convert_phylo(soil, select = phyloseq::taxa_names(soil)[i])

fit_unr1 <- bbdml(formula = cbind(W, M) ~ DayAmdmt, phi.formula = ~ DayAmdmt,
                  data = data_i, nstart = 5)
g1 <- plot_bbdml_nolab(fit_unr1, color = "DayAmdmt")
g11 <- g1 + scale_color_manual(values = hcl(seq(15, 375, length = 4)[2], 
                                            l = c(45, 85), c = 100)) + 
  annotate('text', x = 18, y = 0.030, 
           label = "p[corncob]=='1.00'~x~10^{-6}~(H[0]: beta == 0)",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.028, 
           label = "p[corncob]=='1.00'~x~10^{-6}~(H[0]: beta^'*' == 0)",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.026, label = "p[DESeq2]==3.58~x~10^{-15}",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.024, label = "p[edgeR]==4.96~x~10^{-13}",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.022, label = "p[metagenomeSeq]==0.0543",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.020, label = "p[ZIB]==0.00163",
           parse = TRUE, size = 5, hjust = 0)


i <- 13
data_i <- convert_phylo(soil, select = phyloseq::taxa_names(soil)[i])
fit_unr2 <-  bbdml(formula = cbind(W, M) ~ DayAmdmt, phi.formula = ~ DayAmdmt,
                   data = data_i, nstart = 5)

g2 <- plot_bbdml_nolab(fit_unr2, color = "DayAmdmt")
g22 <- g2 + scale_color_manual(values = hcl(seq(15, 375, length = 4)[1],
                                            l = c(45, 85), c = 100)) + 
  annotate('text', x = 18, y = 0.098, 
           label = "p[corncob]==0.000735~(H[0]: beta == 0)",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.094, 
           label = "p[corncob]==0.403~(H[0]: beta^'*' == 0)",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.09, label = "p[DESeq2]==0.00387",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.086, label = "p[edgeR]==0.143",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.082, label = "p[metagenomeSeq]==0.759",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.078, label = "p[ZIB]==0.000941",
           parse = TRUE, size = 5, hjust = 0)

i <- 88
data_i <- convert_phylo(soil, select = phyloseq::taxa_names(soil)[i])
fit_unr3 <-  bbdml(formula = cbind(W, M) ~ DayAmdmt, phi.formula = ~ DayAmdmt,
                   data = data_i, nstart = 5)
g3 <- plot_bbdml_nolab(fit_unr3, color = "DayAmdmt")
g33 <- g3 + scale_color_manual(values = hcl(seq(15, 375, length = 4)[3],
                                            l = c(45, 85), c = 100)) + 
  annotate('text', x = 18, y = 0.002, 
           label = "p[corncob]==0.243~(H[0]: beta == 0)",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.0019, 
           label = "p[corncob]==0.00871~(H[0]: beta^'*' == 0)",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.0018, label = "p[DESeq2]==0.496",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.0017, label = "p[edgeR]==0.743",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.0016, label = "p[metagenomeSeq]==0.791",
           parse = TRUE, size = 5, hjust = 0) + 
  annotate('text', x = 18, y = 0.0015, label = "p[ZIB]==0.757",
           parse = TRUE, size = 5, hjust = 0)

pdf(file = 'examples.pdf',width = 24, height = 6)
gridExtra::grid.arrange(g11, g22, g33, ncol = 3)
dev.off()
