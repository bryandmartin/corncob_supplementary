load("Figure3Data.RData")

powerPlot <- function(x, title = NULL, nolab = FALSE, bs = 25, star = TRUE) {
  xlab <- expression(paste("|",beta[1],"|"))
  if (star) xlab <- expression(paste("|",beta[1]^{"*"},"|"))
  bval <- -2.46370460866063
  if (star) bval <- -3.87644320689491
  x$cons <- x$cons * bval
  x$cons <- abs(x$cons)
  df <- tidyr::gather(x, key = "Test", value = "pval", wald:pblrt)
  df$n <- factor(df$n, levels = c("100", "30", "10"))
  df$Test <- df$Test %>% gsub("wald", "Wald", .) %>%
    gsub("lrt", "LRT", .) %>%
    gsub("pbWald", "Parametric Bootstrap Wald",.) %>% 
    gsub("pbLRT", "Parametric Bootstrap LRT",.)
  df$Test <- factor(df$Test, levels = c("Wald", "LRT",
                                        "Parametric Bootstrap Wald",
                                        "Parametric Bootstrap LRT"))
  df <- df[-(which(df$n == "10" & df$Test %in% c("Wald", "LRT"))),]
  df <- aggregate(pval ~ Test + cons + n, data = df,
                  FUN = function(x) mean(x < 0.05))
  g <- ggplot2::ggplot(df) + ggplot2::theme_bw(base_size = bs) +
    ggplot2::geom_line(alpha = 1, lwd = 1.1,
                       ggplot2::aes(cons, pval,
                                    color = Test,
                                    lty = n)) + 
    ggplot2::geom_abline(intercept = 0.05, slope = 0, lwd = 1.1, lty = 2) + 
    ggplot2::labs(x = xlab, y = "Power",
                  title = title,
                  color = "Test", lty = "Sample Size") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      face = "bold"),
                   legend.position = "bottom", legend.direction = "horizontal",
                   legend.title = ggplot2::element_text(face = "bold"),
                   legend.key.width = ggplot2::unit(1.7, "cm")) +
    ggplot2::guides(colour = ggplot2::guide_legend(title.position = "top",
                                                   title.hjust = 0.5),
                    lty = ggplot2::guide_legend(title.position = "top",
                                                title.hjust = 0.5))
  if (nolab) {
    g <- g + ggplot2::guides(colour = FALSE, lty = FALSE)
  }
  g
}

get_legend <- function(myggplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

g1 <- powerPlot(tmp1, title = expression(paste(H[0]:~beta[1]==0)),
                nolab = TRUE, bs = 35, star = FALSE)
g2 <- powerPlot(tmp2, title = expression(paste(H[0]:~beta[1]^'*'==0)),
                nolab = TRUE, bs = 35, star = TRUE)
g22 <- powerPlot(tmp2, title = expression(paste(H[0]:~beta[1]^'*'==0)),
                 nolab = FALSE, bs = 35)

myleg <- get_legend(g22)

pdf(file = 'power.pdf',width = 24, height = 8.5)
gridExtra::grid.arrange(g1, g2, myleg, layout_matrix = rbind(c(1,2),
                                                             c(1,2),
                                                             c(1,2),
                                                             c(1,2),
                                                             c(1,2),
                                                             c(1,2),
                                                             c(3,3)))
dev.off()

