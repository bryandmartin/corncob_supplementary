require(corncob)
require(DESeq2)
require(edgeR)
require(metagenomeSeq)
require(ZIBseq)
require(phyloseq)
require(PathoStat)
require(magrittr)
require(tidyr)


load("./cornp.rda")
load("./cornpdv.rda")
corn_p <- unlist(corn_p)
corn_p_dv <- unlist(corn_p_dv)

data(soil_phylo)
soil <- soil_phylo %>% 
  phyloseq::subset_samples(DayAmdmt %in% c(20, 21)) %>%
  phyloseq::tax_glom("Genus")


# DESeq2
DS_1 <- phyloseq_to_deseq2(soil, ~ DayAmdmt)
diagdds_1 <- DESeq(DS_1, test = "Wald", fitType = "parametric")
ds_p <- results(diagdds_1)$pvalue

# EdgeR
ed_2 <- PathoStat::phyloseq_to_edgeR(soil, group = "DayAmdmt")
dgList <- estimateGLMCommonDisp(ed_2)
fit <- glmFit(dgList)
lrt_2 <- glmLRT(fit, coef = 2)
tt_2 <- topTags(lrt_2, n = nrow(ed_2$table), adjust.method = "BH",
                sort.by = "none")
res_2 = tt_2@.Data[[1]]
edge_p <- res_2$PValue

# metagenomeSeq
phyloseq_to_metagenomeSeq <- function(physeq) {
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  OTU = as(otu_table(physeq), "matrix")
  ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
  TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                      row.names = taxa_names(physeq)))
  MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
  MGS = cumNorm(MGS)
  return(MGS)
}

mg_1 <- phyloseq_to_metagenomeSeq(soil)
lungData <- cumNorm(mg_1, p = 0.5)
pd <- pData(mg_1)
mod <- model.matrix(~1 + DayAmdmt, data = pd)
fit_1 <- fitFeatureModel(lungData, mod)
mg_p <- fit_1$pvalues

# ZIB
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
z_dat <- t(as.data.frame(otu_table(soil)))
z_con <- sample_data(soil)$DayAmdmt
ZIBseq2 <- function(data, outcome, transform = F, alpha = 0.05) {
  genodata = data
  Y = outcome
  useF = which(colSums(genodata) > 0 * dim(genodata)[1])
  X = genodata[, useF]
  ST = rowSums(X)
  P = dim(X)[2]
  beta = matrix(data = NA, P, 2)
  for (i in 1:P) {
    x.prop = X[, i]/ST
    if (transform == T) {
      x.prop = sqrt(x.prop)
    }
    bereg = gamlss(x.prop ~ Y, family = BEZI(sigma.link = "identity"), 
                   trace = FALSE, control = gamlss.control(n.cyc = 100))
    out = summary(bereg)
    beta[i, ] = out[2, c(1, 4)]
  }
  pvalues = beta[, 2]
  list(pvalues = pvalues, used = useF)
}


z_1 <- quiet(ZIBseq2(z_dat, z_con, transform = F, alpha = 0.05))
zib_1 <- z_1$pvalues
zib_p <- rep(NA, 241)
zib_p[z_1$used] <- zib_1

my_p <- data.frame(corn = -log10(corn_p),
                   ds = -log10(ds_p),
                   edge = -log10(edge_p),
                   mg = -log10(mg_p),
                   zib = -log10(zib_p),
                   corndv = -log10(corn_p_dv))

save(my_p, file = "myp.rda")

#####
# Getting summary statistics
#####
mysummary <- function(data) {
  data <- data[complete.cases(data),]
  cat("Spearman: ", round(cor(data[,1], data[,2], method = "spearman"), 3),
      "\n")
  cat("Percent less: ", round(sum(data[,1] < data[,2])/nrow(data), 3) * 100,
      "\n")
  cat("Median: ", round(median(data[,2]), 3), "\n")
}

round(median(corn_p), 3)
tmp <- data.frame(corn_p, ds_p); mysummary(tmp)
tmp <- data.frame(corn_p, edge_p); mysummary(tmp)
tmp <- data.frame(corn_p, mg_p); mysummary(tmp)
tmp <- data.frame(corn_p, zib_p); mysummary(tmp)

round(median(corn_p_dv), 3)
tmp <- data.frame(corn_p_dv, ds_p); mysummary(tmp)
tmp <- data.frame(corn_p_dv, edge_p); mysummary(tmp)
tmp <- data.frame(corn_p_dv, mg_p); mysummary(tmp)
tmp <- data.frame(corn_p_dv, zib_p); mysummary(tmp)
tmp <- data.frame(corn_p_dv, corn_p); mysummary(tmp)
