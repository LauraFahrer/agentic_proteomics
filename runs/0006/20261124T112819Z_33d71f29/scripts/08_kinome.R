#!/usr/bin/env Rscript
# Step 08 - kinome landscape (Fig. 6A/B), outlier kinases (Fig. S6A-C), kinome NMF (Fig. S8B/C),
#           40 differentially expressed kinases (Fig. 8B)
# Kinase universe: HGNC/Manning kinase list shipped with InKA (Kinases_HGNC2017.01.20_v2, 536 genes; data
# file only). Kinase group/family: UniProt 'Protein families' comment (rest.uniprot.org, 2026).
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressPackageStartupMessages({ library(NMF); library(limma); library(pheatmap) })
e <- new.env(); load(file.path(SRC, "inka-intensity/inka_20220330/Kinases_HGNC2017.01.20_v2.Rdata"), e)
KIN <- unique(e$Kinases$Gene.Symbol.HGNC21Apr2016)
fam <- fread(file.path(DC, "uniprot/human_kinase_families.tsv"))
setnames(fam, c("acc", "gene", "family"))
group_of <- function(f) {
  if (is.na(f) || f == "") return("Other")
  for (g in c("AGC", "CAMK", "CK1", "CMGC", "STE", "TKL")) if (grepl(paste0("\\b", g, " "), f)) return(g)
  if (grepl("Tyr protein kinase family", f)) return("TK")
  if (grepl("PI3/PI4-kinase|RIO-type|ABC1|PIKK|BCR|Alpha-type|PDK/BCKDK|FAST|TIF1|BRD", f)) return("Atypical")
  "Other"
}
tk_sub <- function(f) {
  m <- regmatches(f, regexpr("Tyr protein kinase family[.,] ?[^.;]*subfamily", f))
  if (!length(m)) return(NA_character_)
  s <- sub(".*family[.,] ?", "", sub(" subfamily", "", m))
  map <- c("Fibroblast growth factor receptor" = "FGFR", "Ephrin receptor" = "EPH", "SRC" = "SRC", "EGF receptor" = "EGFR",
           "JAK" = "JAK", "Insulin receptor" = "INSR", "ABL" = "ABL", "AXL/UFO" = "AXL", "CSF-1/PDGF receptor" = "PDGFR",
           "FAK" = "FAK", "ACK" = "ACK", "Tec" = "TEC", "SYK/ZAP-70" = "SYK", "CSK" = "CSK", "FES/FPS" = "FER",
           "Insulin receptor" = "INSR", "ROR" = "ROR", "Trk" = "TRK", "VEGF receptor" = "VEGFR")
  ifelse(s %in% names(map), map[s], s)
}
fam$group <- sapply(fam$family, group_of); fam$tk_family <- sapply(fam$family, tk_sub)
fam <- fam[!duplicated(gene)]
grp_of <- function(g) { x <- fam$group[match(g, fam$gene)]; x[is.na(x)] <- "Other"; x }
# MET/DDR are TK families in Manning; UniProt names them differently - patch common cases
tkf_of <- function(g) {
  x <- fam$tk_family[match(g, fam$gene)]
  x[g %in% c("MET", "MST1R")] <- "MET"; x[g %in% c("DDR1", "DDR2")] <- "DDR"; x[is.na(x)] <- "Other"; x
}

WP <- load_WP(); PP <- load_PP(); YP <- load_YP()
gW <- attr(WP, "annot")$Gene; gP <- attr(PP, "annot")$GeneSymbol; gY <- attr(YP, "annot")$GeneSymbol
kW <- gW %in% KIN; kP <- gP %in% KIN; kY <- gY %in% KIN
cnt <- data.table(dataset = c("WP", "PP", "YP"), kinases_mapped = c(length(unique(gW[kW])), length(unique(gP[kP])), length(unique(gY[kY]))),
                  kinase_features = c(sum(kW), sum(kP), sum(kY)), paper_kinases = c(198, 193, 140))
fwrite(cnt, file.path(TAB, "08_kinome_counts.csv")); print(cnt)

share <- function(m, genes, keyfun) {
  inten <- rowSums(2^m); k <- keyfun(genes)
  s <- tapply(inten, k, sum); round(100 * s / sum(s), 2)
}
fA_wp <- share(WP[kW, ], gW[kW], grp_of); fA_pp <- share(PP[kP, ], gP[kP], grp_of); fB_yp <- share(YP[kY, ], gY[kY], tkf_of)
paperA_wp <- c(Other = 33.02, Atypical = 22.01, STE = 4.36, TK = 4.03, AGC = 0.48, CAMK = 0.04, CK1 = 1.4e-05, TKL = 0.07, CMGC = 36.05)
paperA_pp <- c(CMGC = 20, Other = 8, Atypical = 23, STE = 11, TK = 6, AGC = 6, CAMK = 20, TKL = 2, CK1 = 3)
paperB_yp <- c(FGFR = 32.97, EPH = 32.89, SRC = 14.67, EGFR = 8.79, JAK = 3.83, MET = 3.37, INSR = 1.18, DDR = 0.77, ABL = 0.67, AXL = 0.3, Other = 0.54)
mk <- function(x, p, d) { n <- union(names(x), names(p)); data.table(dataset = d, family = n, this_run_pct = as.numeric(x[n]), paper_pct = as.numeric(p[n])) }
FS <- rbind(mk(fA_wp, paperA_wp, "WP"), mk(fA_pp, paperA_pp, "PP"), mk(fB_yp, paperB_yp, "YP"))
fwrite(FS, file.path(TAB, "08_Fig6AB_kinase_family_share.csv")); print(FS)
save_plot("08_Fig6AB_kinase_family_share", {
  par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
  for (d in c("WP", "PP", "YP")) { x <- FS[dataset == d & !is.na(this_run_pct) & this_run_pct > 0]
    pie(x$this_run_pct, labels = sprintf("%s %.1f%%", x$family, x$this_run_pct), main = paste(d, "kinome intensity share (this run)"), cex = 0.8) }
}, w = 15, h = 5)

# ---------------- outlier kinases (Fig. S6A-C) ----------------
# 'top 5% of molecular features with the highest expression/phosphorylation level across all samples'.
# Operationalised as: kinase features x samples, values in the top 5% of the dataset's kinome matrix (def A, raw
# level) or of the per-feature z-scores (def B, relative to other lines). Both reported.
outl <- function(m, genes, name) {
  ids <- rownames(m)
  zA <- m; thrA <- quantile(zA, 0.95)
  zB <- (m - rowMeans(m)) / apply(m, 1, sd); thrB <- quantile(zB, 0.95, na.rm = TRUE)
  f <- function(z, thr, def) {
    w <- which(z >= thr, arr.ind = TRUE)
    if (!nrow(w)) return(NULL)
    d <- data.table(feature = ids[w[, 1]], gene = genes[w[, 1]], sample = colnames(m)[w[, 2]], value = z[w])
    d[, line := line_of(sample)]
    d[, .(n_samples = .N, lines = paste(unique(line), collapse = ";")), by = .(feature, gene)][, `:=`(dataset = name, def = def)]
  }
  rbind(f(zA, thrA, "A_top5pct_level"), f(zB, thrB, "B_top5pct_zscore"))
}
O <- rbind(outl(WP[kW, ], gW[kW], "WP"), outl(PP[kP, ], gP[kP], "PP"), outl(YP[kY, ], gY[kY], "YP"))
fwrite(O, file.path(TAB, "08_FigS6_outlier_kinases.csv"))
ok <- O[, .(n_outlier_kinases = uniqueN(gene)), by = def]
ok$union_over_3_datasets <- TRUE; print(ok)
# recovery of the paper's named outliers (only lines where all replicates are outliers count as 'line outlier')
named <- data.table(feature_pattern = c("^ERBB2$", "^CDK12$", "^ERBB2_", "^CDK12_", "^FGFR2_", "^FGFR2_", "^CSNK1D_S384", "^CSNK1E_S389"),
                    dataset = c("WP", "WP", "PP", "PP", "PP", "YP", "PP", "PP"),
                    paper_line = c("N87", "N87", "N87", "N87", "KatoIII", "KatoIII", "MKN28", "MKN28"))
rec <- rbindlist(lapply(seq_len(nrow(named)), function(i) {
  x <- O[dataset == named$dataset[i] & grepl(named$feature_pattern[i], feature)]
  x[, .(dataset = named$dataset[i], pattern = named$feature_pattern[i], paper_line = named$paper_line[i], def,
        n_features = .N, features_outlier_in_paper_line = sum(grepl(named$paper_line[i], lines)),
        lines_flagged = paste(unique(unlist(strsplit(lines, ";"))), collapse = ";"))], by = def][, -1]
}), fill = TRUE)
fwrite(rec, file.path(TAB, "08_FigS6_named_outliers_recovery.csv")); print(rec)
# per-line view of the named features (line means, z-scores)
lm_ <- function(m) sapply(LINES, function(l) rowMeans(m[, line_of(colnames(m)) == l, drop = FALSE]))
show <- list(WP = c("ERBB2", "CDK12", "CDK4", "CDK6", "SRC"),
             PP = grep("^(ERBB2|CDK12|FGFR2|CSNK1D|CSNK1E)_", rownames(PP), value = TRUE),
             YP = grep("^(ERBB2|FGFR2|SRC|EGFR|MET)_", rownames(YP), value = TRUE))
zz <- rbind(t(scale(t(lm_(WP[show$WP, ])))), t(scale(t(lm_(PP[show$PP, ])))), t(scale(t(lm_(YP[show$YP, ])))))
fwrite(data.table(feature = rownames(zz), round(zz, 2)), file.path(TAB, "08_FigS6_named_features_line_zscores.csv"))
save_plot("08_FigS6_named_outlier_features", pheatmap(zz, cluster_cols = FALSE, cluster_rows = FALSE,
          main = "Named outlier kinases: line-mean z-scores (WP / PP / YP)", fontsize_row = 6,
          color = colorRampPalette(c("navy", "white", "firebrick"))(50), breaks = seq(-3.5, 3.5, length.out = 51)), w = 8, h = 11)
# CDK4/6 MKN28 vs YCCEL1 (Fig. S6G)
s6g <- rbindlist(lapply(c("CDK4", "CDK6"), function(g) {
  a <- WP[g, line_of(colnames(WP)) == "YCCEL1"]; b <- WP[g, line_of(colnames(WP)) == "MKN28"]
  data.table(gene = g, YCCEL1_mean = mean(a), MKN28_mean = mean(b), t_test_P = t.test(a, b)$p.value)
}))
s6g$paper_P <- c(1e-4, 2.14e-2); fwrite(s6g, file.path(TAB, "08_FigS6G_CDK46_MKN28_vs_YCCEL1.csv")); print(s6g)

# ---------------- kinome NMF (Fig. S8B/C) ----------------
INT <- fread(file.path(TAB, "03_integrative_clusters.csv")); intlab <- setNames(INT$integrative_subtype, INT$line)
km <- WP[kW, ]
fitk <- nmf(km, rank = 2, method = "brunet", nrun = 100, seed = 123456, .opt = "-p")
clk <- predict(fitk, what = "consensus"); silk <- silhouette(fitk, what = "consensus")
kres <- data.table(sample = names(clk), line = line_of(names(clk)), kinome_cluster = as.integer(clk), silhouette = silk[, "sil_width"],
                   integrative = intlab[line_of(names(clk))])
fwrite(kres, file.path(TAB, "08_FigS8BC_kinome_NMF.csv"))
kin_summary <- list(n_kinases = nrow(km), cophenetic = cophcor(fitk), mean_silhouette = mean(silk[, "sil_width"]),
                    ARI_vs_integrative = ari(kres$kinome_cluster, kres$integrative),
                    contingency = capture.output(print(table(kres$kinome_cluster, kres$integrative))))
jsonlite::write_json(kin_summary, file.path(TAB, "08_FigS8BC_kinome_NMF_summary.json"), auto_unbox = TRUE, pretty = TRUE)
print(kin_summary)
save_plot("08_FigS8B_kinome_NMF_consensus", consensusmap(fitk, annCol = data.frame(integrative = kres$integrative, row.names = kres$sample),
          main = "WP kinome NMF (k=2, 100 runs)"), w = 9, h = 8)

# ---------------- 40 differentially expressed kinases (Fig. 8B) ----------------
g <- factor(intlab[line_of(colnames(km))], levels = c("Metabolism", "EMT"))
f <- eBayes(lmFit(km, model.matrix(~ g))); tt <- topTable(f, coef = 2, number = Inf)
tt$kinase <- rownames(tt)
paper40 <- c("BAZ1B", "STK24", "CASK", "KALRN", "STK39", "COQ8A", "RPS6KA3", "TAOK3", "BRD2", "RPS6KA1", "MST1R", "MET", "ABR", "MAP2K3",
             "MAPK13", "SRC", "EPHB2", "GAK", "PAK1", "CAMK2D", "LYN", "CDK2", "CDK1", "CSNK2A1", "CDK7", "TGM2", "STK3", "BMP2K", "STK4",
             "NRBP1", "SCYL2", "MAP3K7", "CDK9", "YES1", "MAP3K20", "ILK", "BRD4", "PTK7", "CAMK1", "IGF1R")
tt$in_paper40 <- tt$kinase %in% paper40
fwrite(as.data.table(tt), file.path(TAB, "08_Fig8B_kinase_limma.csv"))
top40 <- head(tt$kinase, 40)
p40 <- data.table(kinase = paper40, in_WP = paper40 %in% rownames(km), rank_this_run = match(paper40, tt$kinase),
                  adjP = tt$adj.P.Val[match(paper40, tt$kinase)], log2FC_EMT_vs_Met = tt$logFC[match(paper40, tt$kinase)])
fwrite(p40, file.path(TAB, "08_Fig8B_paper40_check.csv"))
cat("overlap of this run's top-40 with paper's 40:", sum(top40 %in% paper40), "\n")
# clustering of samples with the paper's 40 kinases (Fig. 8B)
x40 <- km[intersect(paper40, rownames(km)), ]
hc <- cutree(hclust(dist(t(x40 - rowMeans(x40)))), 2)
fwrite(data.table(sample = names(hc), cluster = hc, integrative = intlab[line_of(names(hc))]), file.path(TAB, "08_Fig8B_paper40_hclust.csv"))
cat("ARI (hclust on paper's 40 kinases vs integrative):", ari(hc, intlab[line_of(names(hc))]), "\n")
save_plot("08_Fig8B_kinase_heatmap", {
  o <- order(intlab[line_of(colnames(x40))], colnames(x40))
  pheatmap(x40[, o], scale = "row", cluster_cols = FALSE, annotation_col = data.frame(subtype = intlab[line_of(colnames(x40))], row.names = colnames(x40)),
           color = colorRampPalette(c("navy", "black", "yellow"))(50), breaks = seq(-2, 2, length.out = 51), fontsize = 6,
           main = "Paper's 40 kinases (WP, Table S3)")
}, w = 10, h = 8)
