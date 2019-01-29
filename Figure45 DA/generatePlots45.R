require(ggplot2)
require(tidyr)

load("./myp.rda")

my_p$col <- "black"
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

my_p$col[c(13,49,88)] <- gg_color_hue(3)
my_p$size <- 1
my_p$size[my_p$col != "black"] <- 3
df <- my_p %>%
  gather(key, value, ds:zib)
df$key <- factor(df$key, levels = c("ds", "edge", "mg", "zib"), 
                 labels = c(expression(-log[10](p[DESeq2])),
                            expression(-log[10](p[edgeR])),
                            expression(-log[10](p[metagenomeSeq])),
                            expression(-log[10](p[ZIB]))))


pdf("DA_Full.pdf", width = 22, height = 5)
ggplot(df, aes(x = corn, y = value)) + 
  geom_point(color = df$col, size = df$size) +
  geom_abline(slope = 1, lwd = 1.2, col = "#999999") +
  geom_point(data = subset(df, size == 3), color = subset(df, size == 3)$col, 
             size = 6) + 
  facet_wrap(~key, strip.position = "left", scales = "free_y",
             labeller = label_parsed, nrow = 1) +
  labs(x = expression(-log[10](p[corncob])~H[0]:~beta~"="~0), y = NULL) +
  theme_bw(base_size = 35) + ylim(0,33) + xlim(0,6.1) +
  theme(strip.background = element_blank(),
        strip.placement = "outside", panel.spacing = unit(0, "lines"),
        strip.text = element_text(margin = margin(l = 10)))
dev.off()

# Remove edgeR for this plot, add corncob beta = 0
df2 <- my_p %>%
  gather(key, value, corn:zib)
df2$key <- factor(df2$key, levels = c("ds", "edge", "mg", "zib", "corn"),
                  labels = c(expression(-log[10](p[DESeq2])),
                             expression(-log[10](p[edgeR])),
                             expression(-log[10](p[metagenomeSeq])),
                             expression(-log[10](p[ZIB])),
                             expression(-log[10](p[corncob])~H[0]:~beta~"="~0)))

df2 <- df2[-c(483:(483 + 240)),]
pdf("DV_Full.pdf", width = 22, height = 5)
ggplot(df2, aes(x = corndv, y = value)) + 
  geom_point(color = df2$col, size = df2$size) +
  geom_abline(slope = 1, lwd = 1.2, col = "#999999") +
  geom_point(data = subset(df2, size == 3), color = subset(df2, size == 3)$col,
             size = 6) + 
  facet_wrap(~key, strip.position = "left", scales = "free_y",
             labeller = label_parsed, nrow = 1) +
  labs(x = expression(-log[10](p[corncob])~H[0]:~beta^"*"~"="~0), y = NULL) +
  theme_bw(base_size = 30) + ylim(0,33) + xlim(0,6.1) +
  theme(strip.background = element_blank(),
        strip.placement = "outside", panel.spacing = unit(0, "lines"),
        strip.text = element_text(margin = margin(l = 10)))
dev.off()
