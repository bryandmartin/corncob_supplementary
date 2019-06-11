load("Figure2Data.RData")

unifQQ <- function(x, title = NULL, nolab = FALSE, bs = 25) {
  df <- tidyr::gather(x, key = "Test", value = "pval", wald:pblrt)
  df$n <- factor(df$n, levels = c("100", "30", "10"))
  df$Test <- df$Test %>% gsub("wald", "Wald", .) %>% gsub("lrt", "LRT", .) %>%
    gsub("pbWald", "Parametric Bootstrap Wald",.) %>%
    gsub("pbLRT", "Parametric Bootstrap LRT",.)
  df$Test <- factor(df$Test, 
                    levels = c("Wald", "LRT", "Parametric Bootstrap Wald",
                               "Parametric Bootstrap LRT"))
  g <- ggplot2::ggplot(df) + ggplot2::theme_bw(base_size = bs) +
    ggplot2::geom_qq(distribution = stats::qunif, geom = "line", alpha = 1,
                     lwd = 1.1, ggplot2::aes(sample = pval, color = Test, 
                                             lty = n)) + 
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::labs(x = "Theoretical", y = "Sample", title = title,
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

g1 <- unifQQ(tmp1, title = expression(paste(H[0]:~"(",beta[1],", ",beta[1]^'*',
                                            ")",phantom()=="(",0,", ",0,")")),
             nolab = TRUE, bs = 35) 
g2 <- unifQQ(tmp2, title = expression(paste(H[0]:~beta[1]^'*'==0)),
             nolab = TRUE, bs = 35) 
g3 <- unifQQ(tmp3, title = expression(paste(H[0]:~beta[1]==0)),
             nolab = TRUE, bs = 35) 
g33 <- unifQQ(tmp3, title = expression(paste(H[0]:~beta[1]==0)),
              nolab = FALSE, bs = 35) 

myleg <- get_legend(g33)

pdf(file = 'typeI.pdf', width = 24, height = 8.5)
gridExtra::grid.arrange(g1, g2, g3, myleg, layout_matrix = rbind(c(1,2,3),
                                                            c(1,2,3),
                                                            c(1,2,3),
                                                            c(1,2,3),
                                                            c(1,2,3),
                                                            c(1,2,3),
                                                            c(4,4,4)))
dev.off()

#####
# Code to get results for Table
#####

# for (i in 5:8) {
#   tmp <- round(quantile(tmp1[1:10000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# for (i in 5:8) {
#   tmp <- round(quantile(tmp1[10001:20000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# for (i in 5:8) {
#   tmp <- round(quantile(tmp1[20001:30000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# 
# for (i in 5:8) {
#   tmp <- round(quantile(tmp2[1:10000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# for (i in 5:8) {
#   tmp <- round(quantile(tmp2[10001:20000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# for (i in 5:8) {
#   tmp <- round(quantile(tmp2[20001:30000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# 
# for (i in 5:8) {
#   tmp <- round(quantile(tmp3[1:10000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# for (i in 5:8) {
#   tmp <- round(quantile(tmp3[10001:20000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
# for (i in 5:8) {
#   tmp <- round(quantile(tmp3[20001:30000, i],
#                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)
#   tmp <- as.character(tmp)
#   print(paste(paste(tmp, collapse = " & "), "\\"))
# }
