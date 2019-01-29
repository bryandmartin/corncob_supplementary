require(ggplot2)
require(magrittr)

load("./cornp.rda")
load("./cornpdv.rda")
corn_p <- unlist(corn_p)
corn_p_dv <- unlist(corn_p_dv)

tmp <- p.adjust(corn_p, method = "fdr")
tmp <- tmp[order(tmp)]
xx <- yy <- 1:241
for (i in xx) {
  yy[i] <- tmp[i]
}
xx <- c(0,xx); yy <- c(0,yy)
tmp <- p.adjust(corn_p_dv, method = "fdr")
tmp <- tmp[order(tmp)]
yy2 <- 1:241
for (i in xx) {
  yy2[i] <- tmp[i]
}
yy2 <- c(0, yy2)
mydf <- data.frame(x = xx, y = yy, y2 = yy2)
mydf <- reshape2::melt(mydf, id = "x")

pdf(file = 'FDR.pdf',width = 12, height = 5)
ggplot(mydf, aes(x = x, y = value, color = variable)) +
  geom_line(lwd = 1.2) +
  theme_bw(base_size = 25) + 
  scale_x_continuous(limits = c(0, 228), expand = c(0, 0)) + 
  ggplot2::labs(title = NULL, colour = "Test") +
  xlab("Number of Genera") + ylab("Estimated FDR") + 
  scale_color_discrete(labels = c(expression(paste("Differential Abundance ",
                                                   (H[0]: beta == 0))),
                                  expression(paste("Differential Variability ",
                                                   (H[0]: beta^'*' == 0)))))
dev.off()
