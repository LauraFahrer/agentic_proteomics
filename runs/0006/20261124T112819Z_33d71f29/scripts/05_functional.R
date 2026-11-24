#!/usr/bin/env Rscript
# Step 05 - functional annotation (paper Methods 2.11; Fig. 2B, 4B, 4D, 5B, S2C/D, S4C/D/G; Table S6)
#  (a) GSEA desktop/CLI (paper 4.0.3; here 4.3.2 from bioconda - 4.0.3 not obtainable) with
#      MSigDB v7.0 gene sets (bundled in broadinstitute/ssGSEA2.0; v7.0 ~ GSEA 4.0.3 era)
#      WP: Hallmark, each cluster vs rest; PP: KEGG, cluster 1 vs 2; YP: KEGG, each cluster vs rest.
#      Phosphosite -> gene collapse by GSEA's own Max_probe mode (CHIP file).
#  (b) WGCNA co-expression modules (WP, PP), module-cluster correlation, GO BP enrichment
#  (c) one-way ANOVA (WP, YP clusters) / limma (PP clusters) + GO BP over-representation
#      (geneontology.org/PANTHER unreachable -> clusterProfiler::enricher on MSigDB v7.0 c5.bp)
#  (d) integrative EMT vs Metabolism limma (|log2FC|>1, adj.P<0.05) on WP and PP (Fig. 5B)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressPackageStartupMessages({ library(limma); library(clusterProfiler); library(WGCNA); library(pheatmap) })
allowWGCNAThreads(as.integer(Sys.getenv("NCORES", "4")))
MS <- file.path(SRC, "ssGSEA2.0/db/msigdb")
GMT <- list(H = file.path(MS, "h.all.v7.0.symbols.gmt"), KEGG = file.path(MS, "c2.cp.kegg.v7.0.symbols.gmt"),
            BP = file.path(MS, "c5.bp.v7.0.symbols.gmt"))
read_gmt2 <- function(f) { x <- clusterProfiler::read.gmt(f); colnames(x) <- c("term", "gene"); x }
T2G <- lapply(GMT, read_gmt2)
GD <- file.path(DC, "gsea"); dir.create(GD, showWarnings = FALSE)
GSEA_BIN <- Sys.getenv("GSEA_BIN", "/opt/mm/envs/repro/bin/gsea-cli.sh")
if (!file.exists(GSEA_BIN)) GSEA_BIN <- file.path(dirname(Sys.which("R")), "../share/gsea-4.3.2-0/gsea-cli.sh")

A <- fread(file.path(TAB, "04_subtype_assignments.csv"))
INT <- fread(file.path(TAB, "03_integrative_clusters.csv")); intlab <- setNames(INT$integrative_subtype, INT$line)
WP <- load_WP(); PP <- load_PP(); YP <- load_YP()
gene_of <- list(WP = attr(WP, "annot")$Gene, PP = attr(PP, "annot")$GeneSymbol, YP = attr(YP, "annot")$GeneSymbol)
names(gene_of$WP) <- rownames(WP); names(gene_of$PP) <- rownames(PP); names(gene_of$YP) <- rownames(YP)
labs <- function(d, m) { a <- A[dataset == d]; setNames(paste0("C", a$cluster), a$sample)[colnames(m)] }

# ---------------- (a) GSEA ----------------
write_gct <- function(m, f) {
  con <- file(f, "w"); writeLines(c("#1.2", paste(nrow(m), ncol(m), sep = "\t"),
    paste(c("NAME", "Description", colnames(m)), collapse = "\t")), con)
  close(con); fwrite(data.table(NAME = rownames(m), Description = "na", m), f, sep = "\t", append = TRUE, col.names = FALSE)
}
write_cls <- function(lab, f) writeLines(c(paste(length(lab), length(unique(lab)), 1, sep = " "),
  paste("#", paste(unique(lab), collapse = " ")), paste(lab, collapse = " ")), f)
write_chip <- function(ids, genes, f) fwrite(data.table(`Probe Set ID` = ids, `Gene Symbol` = genes, `Gene Title` = "na"), f, sep = "\t")
run_gsea <- function(m, lab, target, gmt, label, chip = NULL) {
  out <- file.path(GD, label); unlink(out, recursive = TRUE); dir.create(out)
  gct <- file.path(out, "data.gct"); cls <- file.path(out, "pheno.cls")
  l2 <- ifelse(lab == target, target, "REST")
  write_gct(m, gct); write_cls(l2, cls)
  args <- c("GSEA", "-res", gct, "-cls", paste0(cls, "#", target, "_versus_REST"), "-gmx", gmt,
            "-collapse", ifelse(is.null(chip), "No_Collapse", "Collapse"), "-mode", "Max_probe",
            "-norm", "meandiv", "-nperm", "1000", "-permute", "phenotype", "-rnd_type", "no_balance",
            "-scoring_scheme", "weighted", "-metric", "Signal2Noise", "-sort", "real", "-order", "descending",
            "-create_gcts", "false", "-create_svgs", "false", "-include_only_symbols", "true", "-make_sets", "true",
            "-median", "false", "-num", "100", "-plot_top_x", "0", "-rnd_seed", "149", "-save_rnd_lists", "false",
            "-set_max", "500", "-set_min", "15", "-zip_report", "false", "-out", out, "-rpt_label", label)
  if (!is.null(chip)) args <- c(args, "-chip", chip)
  system2(GSEA_BIN, args, stdout = file.path(out, "log.txt"), stderr = file.path(out, "log.txt"))
  rd <- list.dirs(out, recursive = FALSE); rd <- rd[grepl(label, basename(rd))][1]
  pos <- list.files(rd, paste0("^gsea_report_for_", target, "_.*\\.tsv$"), full.names = TRUE)
  if (!length(pos)) { warning("GSEA failed: ", label); return(NULL) }
  r <- fread(pos[1]); r <- r[, .(NAME, SIZE, ES, NES, `NOM p-val`, `FDR q-val`)]
  r$target <- target; r$analysis <- label; r
}
gsea_all <- list()
# WP: Hallmark, one-vs-rest (Fig. 2B / Table S6)
wl <- labs("WP", WP)
for (t in sort(unique(wl))) gsea_all[[paste0("WP_H_", t)]] <- run_gsea(WP, wl, t, GMT$H, paste0("WP_H_", t))
# paper labels as control for the functional step alone
wlp <- setNames(paste0("C", PAPER_WP[line_of(colnames(WP))]), colnames(WP))
for (t in sort(unique(wlp))) gsea_all[[paste0("WPpaperlab_H_", t)]] <- run_gsea(WP, wlp, t, GMT$H, paste0("WPpaperlab_H_", t))
# PP: KEGG, cluster vs other (Fig. 4B)
chip_pp <- file.path(GD, "PP.chip"); write_chip(rownames(PP), gene_of$PP, chip_pp)
pl <- labs("PP", PP)
for (t in sort(unique(pl))) gsea_all[[paste0("PP_KEGG_", t)]] <- run_gsea(PP, pl, t, GMT$KEGG, paste0("PP_KEGG_", t), chip_pp)
# YP: KEGG, one-vs-rest (Fig. 4D)
chip_yp <- file.path(GD, "YP.chip"); write_chip(rownames(YP), gene_of$YP, chip_yp)
yl <- labs("YP", YP)
for (t in sort(unique(yl))) gsea_all[[paste0("YP_KEGG_", t)]] <- run_gsea(YP, yl, t, GMT$KEGG, paste0("YP_KEGG_", t), chip_yp)
G <- rbindlist(gsea_all, fill = TRUE)
fwrite(G, file.path(TAB, "05_gsea_all.csv"))
cat("GSEA significant (FDR<0.25) per analysis:\n"); print(G[`FDR q-val` < 0.25, .N, by = analysis])

# compare WP hallmark with paper Table S6
s6 <- as.data.table(readxl::read_excel(file.path(DC, "supp/x/MOL2-18-1437-s004.xlsx"), skip = 2))
paper_s6 <- rbind(data.table(set = s6[[1]], group = "Cell prolif/EMT"), data.table(set = s6[[3]], group = "EMT"),
                  data.table(set = s6[[5]], group = "Metabolism"))[!is.na(set) & grepl("^HALLMARK", set)]
fwrite(paper_s6, file.path(TAB, "05_paper_TableS6_hallmarks.csv"))

# ---------------- (b) WGCNA ----------------
wgcna_run <- function(m, name, clus) {
  x <- t(m)
  sft <- pickSoftThreshold(x, powerVector = 1:20, verbose = 0)
  fi <- sft$fitIndices; pw <- fi$Power[which(-sign(fi$slope) * fi$SFT.R.sq >= 0.8)[1]]
  if (is.na(pw)) pw <- 6
  net <- blockwiseModules(x, power = pw, TOMType = "unsigned", minModuleSize = 30, mergeCutHeight = 0.25,
                          numericLabels = FALSE, maxBlockSize = 10000, randomSeed = 12345, verbose = 0)
  mods <- net$colors; ME <- net$MEs
  trait <- sapply(sort(unique(clus)), function(c) as.integer(clus == c)); colnames(trait) <- sort(unique(clus))
  cr <- cor(ME, trait); pv <- corPvalueStudent(cr, nrow(x))
  uni <- unique(gene_of[[name]])
  enr <- rbindlist(lapply(setdiff(unique(mods), "grey"), function(mc) {
    g <- unique(gene_of[[name]][names(mods)[mods == mc]])
    e <- tryCatch(as.data.table(enricher(g, universe = uni, TERM2GENE = T2G$BP, pvalueCutoff = 0.05, qvalueCutoff = 1)@result),
                  error = function(err) NULL)
    if (is.null(e) || !nrow(e)) return(NULL)
    e <- e[p.adjust < 0.05][order(p.adjust)][seq_len(min(.N, 5))]
    if (!nrow(e)) return(NULL)
    e$module <- mc; e[, .(module, ID, GeneRatio, p.adjust)]
  }), fill = TRUE)
  sz <- table(mods)
  summ <- data.table(module = sub("^ME", "", rownames(cr)), size = as.integer(sz[sub("^ME", "", rownames(cr))]),
                     round(cr, 3), setNames(as.data.table(signif(pv, 3)), paste0("p_", colnames(pv))))
  summ$top_GO_BP <- sapply(summ$module, function(mc) paste(head(enr[module == mc]$ID, 3), collapse = "; "))
  fwrite(summ, file.path(TAB, sprintf("05_WGCNA_%s_modules.csv", name)))
  fwrite(enr, file.path(TAB, sprintf("05_WGCNA_%s_module_GO.csv", name)))
  save_plot(sprintf("05_WGCNA_%s_module_cluster_correlation", name), {
    labeledHeatmap(cr, xLabels = colnames(cr), yLabels = rownames(cr), ySymbols = paste0(rownames(cr), " (", summ$size, ")"),
                   colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = paste0(signif(cr, 2), "\n(", signif(pv, 1), ")"),
                   setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1),
                   main = sprintf("%s WGCNA modules vs NMF clusters (power=%d)", name, pw))
  }, w = 7, h = 9)
  list(power = pw, n_modules = length(setdiff(unique(mods), "grey")), summary = summ)
}
par(mar = c(6, 10, 3, 3))
wg_wp <- wgcna_run(WP, "WP", labs("WP", WP))
wg_pp <- wgcna_run(PP, "PP", labs("PP", PP))
jsonlite::write_json(list(WP = list(power = wg_wp$power, n_modules = wg_wp$n_modules, paper_n_modules = 12),
                          PP = list(power = wg_pp$power, n_modules = wg_pp$n_modules, paper_n_modules = 5)),
                     file.path(TAB, "05_WGCNA_summary.json"), auto_unbox = TRUE, pretty = TRUE)

# ---------------- (c) ANOVA / limma + GO BP ----------------
go_by_cluster <- function(m, cl, name, sig) {
  uni <- unique(gene_of[[name]])
  means <- sapply(sort(unique(cl)), function(c) rowMeans(m[sig, cl == c, drop = FALSE]))
  top <- colnames(means)[apply(means, 1, which.max)]
  rbindlist(lapply(sort(unique(cl)), function(c) {
    g <- unique(gene_of[[name]][sig[top == c]])
    e <- tryCatch(as.data.table(enricher(g, universe = uni, TERM2GENE = T2G$BP, pvalueCutoff = 0.05, qvalueCutoff = 1)@result),
                  error = function(err) NULL)
    if (is.null(e) || !nrow(e)) return(NULL)
    e <- e[p.adjust < 0.05]; if (!nrow(e)) return(NULL)
    e$cluster <- c; e$n_features_up <- sum(top == c); e[, .(cluster, n_features_up, ID, GeneRatio, p.adjust)]
  }), fill = TRUE)
}
anova_p <- function(m, cl) apply(m, 1, function(v) summary(aov(v ~ factor(cl)))[[1]][["Pr(>F)"]][1])
wl <- labs("WP", WP); pa <- p.adjust(anova_p(WP, wl), "BH")
sig_wp <- names(pa)[pa < 0.05]
go_wp <- go_by_cluster(WP, wl, "WP", sig_wp)
fwrite(data.table(id = names(pa), adjP = pa), file.path(TAB, "05_WP_ANOVA.csv.gz"))
fwrite(go_wp, file.path(TAB, "05_WP_ANOVA_GO_BP.csv"))
yl <- labs("YP", YP); pay <- p.adjust(anova_p(YP, yl), "BH")
sig_yp <- names(pay)[pay < 0.05]
go_yp <- go_by_cluster(YP, yl, "YP", sig_yp)
fwrite(data.table(id = names(pay), adjP = pay), file.path(TAB, "05_YP_ANOVA.csv.gz"))
fwrite(go_yp, file.path(TAB, "05_YP_ANOVA_GO_BP.csv"))
pl <- labs("PP", PP); des <- model.matrix(~ 0 + factor(pl)); colnames(des) <- sort(unique(pl))
fit <- eBayes(contrasts.fit(lmFit(PP, des), makeContrasts(contrasts = paste(colnames(des)[1], "-", colnames(des)[2]), levels = des)))
tt <- topTable(fit, number = Inf, sort.by = "none")
sig_pp <- rownames(tt)[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1]
go_pp <- go_by_cluster(PP, pl, "PP", sig_pp)
fwrite(go_pp, file.path(TAB, "05_PP_limma_GO_BP.csv"))
jsonlite::write_json(list(WP_ANOVA_sig = length(sig_wp), YP_ANOVA_sig = length(sig_yp), PP_limma_sig = length(sig_pp)),
                     file.path(TAB, "05_differential_counts.json"), auto_unbox = TRUE, pretty = TRUE)

# ---------------- (d) integrative EMT vs Metabolism (Fig. 5B) ----------------
limma2 <- function(m, name) {
  g <- factor(intlab[line_of(colnames(m))], levels = c("EMT", "Metabolism"))
  des <- model.matrix(~ g)
  f <- eBayes(lmFit(m, des)); t <- topTable(f, coef = 2, number = Inf, sort.by = "none")
  data.table(id = rownames(m), gene = gene_of[[name]], log2FC_EMT_vs_Met = -t$logFC, adjP = t$adj.P.Val, dataset = name)
}
I <- rbind(limma2(WP, "WP"), limma2(PP, "PP"), limma2(YP, "YP"))
I[, sig := adjP < 0.05 & abs(log2FC_EMT_vs_Met) > 1]
fwrite(I, file.path(TAB, "05_integrative_limma_all.csv.gz"))
named <- data.table(gene = c("AAMP", "CCN1", "LGALS1", "TGFBI", "STK3", "STK4", "TGFB1I1", "CTNNB1", "YAP1", "ABL2", "APC", "PHLDB2",
                              "GLUL", "ALAD", "ACAA1", "ACAA2", "GCLC", "AK3", "G6PD", "GSTA2", "ALDH3A1", "ACSL5", "AKR1A1"),
                    paper_group = c(rep("EMT", 12), rep("Metabolism", 11)),
                    paper_level = c("protein", "protein", "protein", "protein", "protein", "protein", "phospho", "phospho", "phospho", "phospho",
                                    "phospho", "both", "protein", "protein", "protein", "protein", "protein", "protein", "phospho", "phospho",
                                    "both", "both", "both"))
chk <- rbindlist(lapply(seq_len(nrow(named)), function(i) {
  g <- named$gene[i]; r <- I[gene == g]
  w <- r[dataset == "WP"]; p <- r[dataset == "PP"]
  data.table(named[i], WP_log2FC = if (nrow(w)) round(w$log2FC_EMT_vs_Met[1], 2) else NA,
             WP_sig = if (nrow(w)) w$sig[1] else NA,
             PP_n_sites = nrow(p), PP_sites_sig_up_in_paper_group = sum(p$sig & sign(p$log2FC_EMT_vs_Met) == ifelse(named$paper_group[i] == "EMT", 1, -1)))
}))
chk[, WP_direction_agrees := sign(WP_log2FC) == ifelse(paper_group == "EMT", 1, -1)]
fwrite(chk, file.path(TAB, "05_Fig5B_named_proteins_check.csv")); print(chk)
# pathway mapping of significant WP+PP genes (KEGG + Hallmark over-representation)
ora <- rbindlist(lapply(c("EMT", "Metabolism"), function(grp) {
  s <- I[sig == TRUE & dataset %in% c("WP", "PP") & sign(log2FC_EMT_vs_Met) == ifelse(grp == "EMT", 1, -1)]
  g <- unique(s$gene); uni <- unique(I[dataset %in% c("WP", "PP")]$gene)
  rbindlist(lapply(c("KEGG", "H"), function(db) {
    e <- tryCatch(as.data.table(enricher(g, universe = uni, TERM2GENE = T2G[[db]], pvalueCutoff = 1, qvalueCutoff = 1)@result), error = function(err) NULL)
    if (is.null(e)) return(NULL); e$group <- grp; e$db <- db; e$n_genes_input <- length(g); e[p.adjust < 0.25]
  }), fill = TRUE)
}), fill = TRUE)
fwrite(ora, file.path(TAB, "05_Fig5B_integrative_pathway_ORA.csv"))
