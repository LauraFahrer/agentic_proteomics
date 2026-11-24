#!/usr/bin/env Rscript
# Step 04 - summarise NMF grid, final subtypes, integrative ConsensusClusterPlus clustering
# (paper Fig. 2A, 4A, 4C, 5A, S2A/B, S3A, S4A/B/E/F)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressPackageStartupMessages({ library(NMF); library(ConsensusClusterPlus); library(cluster) })
ND <- file.path(DC, "nmf")

grid <- rbindlist(lapply(list.files(TAB, "^03_nmf_grid_.*_mean.csv$", full.names = TRUE), fread))
fwrite(grid[, !"clusters"], file.path(TAB, "04_nmf_grid_all.csv"))

# ---- Fig S2A / S4A / S4E: cophenetic vs k by feature fraction ----
save_plot("04_FigS2A_S4AE_cophenetic_grid", {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (d in c("WP", "PP", "YP", "WPq")) {
    g <- grid[dataset == d]
    if (!nrow(g)) next
    fr <- sort(unique(g$frac)); cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
    plot(NA, xlim = c(2, 5), ylim = range(c(g$cophenetic, 0.85, 1)), xlab = "k", ylab = "cophenetic coefficient",
         main = c(WP = "WP (Table S3)", PP = "PP (Table S3)", YP = "YP (Table S3)", WPq = "WP (quantms)")[d], xaxt = "n")
    axis(1, 2:5)
    for (i in seq_along(fr)) { s <- g[frac == fr[i]][order(k)]; lines(s$k, s$cophenetic, type = "b", col = cols[i], pch = 19) }
    legend("bottomleft", paste0(fr * 100, "% (n=", g[k == 2][order(frac)]$n_features, ")"), col = cols, pch = 19, cex = 0.7, bty = "n", title = "top MAD")
  }
}, w = 10, h = 9)

best <- grid[, .SD[which.max(cophenetic)], by = dataset]
best_k3 <- grid[k > 2, .SD[which.max(cophenetic)], by = dataset]
paper_set <- data.table(dataset = c("WP", "PP", "YP", "WPq"), frac = c(0.1, 0.2, 1.0, 0.1), k = c(3, 2, 3, 3),
                        paper_n_features = c(568, 1691, 4780, NA))
sel <- merge(paper_set, grid, by = c("dataset", "frac", "k"))
sel$rank_of_cophenetic_in_grid <- sapply(seq_len(nrow(sel)), function(i) {
  g <- grid[dataset == sel$dataset[i]]; sum(g$cophenetic > sel$cophenetic[i]) + 1 })
fwrite(sel[, !"clusters"], file.path(TAB, "04_paper_parameters_in_grid.csv"))
fwrite(best[, !"clusters"], file.path(TAB, "04_grid_best_cophenetic.csv"))
fwrite(best_k3[, !"clusters"], file.path(TAB, "04_grid_best_cophenetic_k_gt2.csv"))
print(sel[, !"clusters"]); print(best[, !"clusters"]); print(best_k3[, !"clusters"])

paper_lab <- list(WP = PAPER_WP, PP = PAPER_PP, YP = PAPER_YP, WPq = PAPER_WP)
load_fit <- function(d, fr, k) readRDS(file.path(ND, sprintf("%s_mean_%03d_k%d.rds", d, round(fr * 100), k)))

assign <- list(); comp <- list()
for (i in seq_len(nrow(sel))) {
  d <- sel$dataset[i]; fr <- sel$frac[i]; k <- sel$k[i]
  o <- load_fit(d, fr, k); fit <- o$fit
  cl <- predict(fit, what = "consensus"); sil <- silhouette(fit, what = "consensus")
  smp <- names(cl)
  pl <- paper_lab[[d]][line_of(smp)]
  tab <- table(reproduced = cl, paper = pl)
  comp[[d]] <- data.frame(dataset = d, frac = fr, k = k, n_features = length(o$features),
                          ARI_vs_paper = ari(cl, pl), cophenetic = cophcor(fit),
                          mean_silhouette = mean(sil[, "sil_width"]),
                          contingency = paste(capture.output(print(tab)), collapse = " | "))
  assign[[d]] <- data.frame(dataset = d, sample = smp, line = line_of(smp), cluster = as.integer(cl),
                            silhouette = sil[, "sil_width"], paper_cluster = pl)
  ann <- data.frame(paper = factor(pl), row.names = smp)
  save_plot(sprintf("04_consensus_%s_%dpct_k%d", d, round(fr * 100), k), {
    consensusmap(fit, annRow = ann, annCol = ann, main = sprintf("%s: top %d%% (n=%d), k=%d, coph=%.3f",
                 d, round(fr * 100), length(o$features), k, cophcor(fit)), fontsize = 7)
  }, w = 9, h = 8)
  save_plot(sprintf("04_silhouette_%s_%dpct_k%d", d, round(fr * 100), k), plot(sil, main = paste(d, "silhouette (consensus)")), w = 7, h = 7)
}
comp <- do.call(rbind, comp); fwrite(comp, file.path(TAB, "04_subtypes_vs_paper.csv")); print(comp)
A <- do.call(rbind, assign); fwrite(A, file.path(TAB, "04_subtype_assignments.csv"))

# ---- median-AD sensitivity for the paper parameter sets (computed if grid files exist) ----
med <- list.files(TAB, "^03_nmf_grid_.*_median.csv$", full.names = TRUE)
if (length(med)) fwrite(rbindlist(lapply(med, fread))[, !"clusters"], file.path(TAB, "04_nmf_grid_median_ad.csv"))

# ---- integrative clustering (Methods 2.10; Fig. 5A) ----
# binary subtype x sample table, 2 replicates per line (r1, r2) to align with YP
smp28 <- intersect(A$sample[A$dataset == "YP"], A$sample[A$dataset == "WP"])
bin <- NULL
for (d in c("WP", "PP", "YP")) {
  a <- A[A$dataset == d, ]; a <- a[match(smp28, a$sample), ]
  for (c in sort(unique(a$cluster))) bin <- rbind(bin, setNames(as.integer(a$cluster == c), smp28))
  rownames(bin)[(nrow(bin) - length(unique(a$cluster)) + 1):nrow(bin)] <- paste0(d, "subtype", sort(unique(a$cluster)))
}
fwrite(data.frame(subtype = rownames(bin), bin, check.names = FALSE), file.path(TAB, "04_integrative_binary_table.csv"))
ccp_dir <- file.path(DC, "ccp"); dir.create(ccp_dir, showWarnings = FALSE)
ccp <- ConsensusClusterPlus(bin, maxK = 4, reps = 500, pItem = 0.8, pFeature = 1, clusterAlg = "hc",
                            distance = "pearson", seed = 123456, plot = "pdf", title = ccp_dir)
ic <- ccp[[2]]$consensusClass
icl <- calcICL(ccp, title = ccp_dir, plot = "pdf")
# name integrative clusters by overlap with the WP metabolism-like cluster of the paper
met_like <- names(which.max(tapply(PAPER_INT[line_of(names(ic))] == "Metabolism", ic, mean)))
int_lab <- ifelse(ic == as.integer(met_like), "Metabolism", "EMT")
IS <- data.frame(sample = names(ic), line = line_of(names(ic)), integrative_cluster = ic, integrative_subtype = int_lab,
                 paper_integrative = PAPER_INT[line_of(names(ic))])
fwrite(IS, file.path(TAB, "04_integrative_assignments_samples.csv"))
lines_int <- IS[!duplicated(IS$line), c("line", "integrative_subtype", "paper_integrative")]
cons <- tapply(IS$integrative_subtype, IS$line, function(x) length(unique(x)) == 1)
lines_int$replicates_concordant <- cons[lines_int$line]
fwrite(lines_int, file.path(TAB, "03_integrative_clusters.csv"))
print(lines_int)
cat("integrative ARI vs paper:", ari(IS$integrative_subtype, IS$paper_integrative), "\n")
# consensus stats for k=2..4
ccs <- do.call(rbind, lapply(2:4, function(k) {
  cm <- ccp[[k]]$consensusMatrix; cl <- ccp[[k]]$consensusClass
  data.frame(k = k, mean_within_consensus = mean(sapply(unique(cl), function(c) mean(cm[cl == c, cl == c]))),
             sizes = paste(table(cl), collapse = "/"))
}))
fwrite(ccs, file.path(TAB, "04_integrative_ccp_k_stats.csv"))

save_plot("04_Fig5A_integrative_binary_heatmap", {
  o <- order(IS$integrative_subtype, IS$line)
  ann <- data.frame(Integrated = IS$integrative_subtype, Paper = IS$paper_integrative, row.names = IS$sample)
  pheatmap::pheatmap(bin[, IS$sample[o]], cluster_cols = FALSE, cluster_rows = FALSE, color = c("black", "#e8312f"),
                     annotation_col = ann, legend = FALSE, main = "Integrative clustering (ConsensusClusterPlus, 500 reps, k=2)",
                     fontsize_col = 7)
}, w = 11, h = 4.5)

# ---- WP metagene-specific proteins (Fig. S3A/B) vs the paper's 42 ----
paper42 <- c("UBE2E2", "CD55", "LGALS1", "MOV10", "SLC2A3", "PTK7", "CRABP2", "ITGA2", "ANXA3", "EPS8L2", "SH3KBP1",
             "TACSTD2", "APOBEC3B", "ISG15", "BCAM", "AKR1C3", "STARD10", "AGR2", "AGR3", "ALDH3A1", "ALDH3A2", "VIL1",
             "S100A6", "ITGB4", "CTSE", "LGALS3BP", "LYZ", "LGALS4", "FXYD3", "ANXA5", "ALDH1A1", "LZTS1", "ANKRD45",
             "PHGDH", "ANXA6", "PFN2", "MAGEA4", "TGM2", "VIM", "FSCN1", "TPD52L1", "FKBP10")
o <- load_fit("WP", 0.1, 3); feats <- extractFeatures(o$fit)
mg <- unlist(lapply(seq_along(feats), function(i) if (all(is.na(feats[[i]]))) character(0) else rownames(basis(o$fit))[feats[[i]]]))
mgdf <- data.frame(protein = mg, metagene = rep(seq_along(feats), sapply(feats, function(f) sum(!is.na(f)))))
mgdf$in_paper_42 <- mgdf$protein %in% paper42
fwrite(mgdf, file.path(TAB, "04_WP_metagene_features.csv"))
cat("metagene features:", nrow(mgdf), " overlap with paper's 42:", sum(mgdf$in_paper_42),
    " paper42 present in top-10% feature set:", sum(paper42 %in% o$features), "\n")
