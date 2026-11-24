#!/usr/bin/env Rscript
# Step 06 - PP kinase-activity inference, 6 activity matrices (paper Methods 2.12; Fig. 8A, Fig. S8A)
#   {raw PP, protein-normalised PP} x {mean fold difference, KSEA, MLR}, all with the original
#   PhosMap implementations (github ecnuzdd/PhosMap @3545398, functions sourced unmodified) and
#   PhosMap's bundled human kinase-substrate table (human_ksrr.csv).
#   Groups: integrative EMT vs Metabolism subclusters.
#   Differential activity: two-sided t-test per kinase, P < 0.05 and |log2FC| > 1 (paper);
#   number of matrices (0-6) in which a kinase is assigned to a subgroup is counted.
# Assumptions (not specified in the paper; documented in REPORT.md):
#   * input to all three methods = per-site log2 intensity centred on the panel mean (log2 ratio
#     to the average cell line) - KSEA requires ratios; centring cancels for mean-fold group FCs.
#   * protein normalisation = log2(PP) - log2(WP) of the same gene (sites w/o protein dropped).
#   * log2FC for mean method = log2 of ratio of group means of 2^activity; for KSEA (signed
#     -log10 P) and MLR (ridge coefficients) the scores are already on a log-like scale, so
#     log2FC is taken as the difference of group means.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
ncores <- as.integer(Sys.getenv("NCORES", "4"))

# --- original PhosMap code ---------------------------------------------------
system.file <- function(..., package = "base", mustWork = FALSE) {
  if (package == "PhosMap") file.path(SRC, "PhosMap/inst", ...) else base::system.file(..., package = package, mustWork = mustWork)
}
for (f in c("get_ka_by_mean_or_mlr.R", "get_ksea_result_list.R", "get_kses.R", "compute_kses.R"))
  source(file.path(SRC, "PhosMap/R", f))
ksrr <- read.csv(file.path(SRC, "PhosMap/inst/extdata/kinase_substrate_regulation_relationship_table/human/human_ksrr.csv"),
                 header = TRUE, stringsAsFactors = NA)

groups_file <- file.path(TAB, "03_integrative_clusters.csv")
if (file.exists(groups_file)) {
  g <- fread(groups_file); grp <- setNames(g$integrative_subtype, g$line)
  cat("using reproduced integrative clusters\n")
} else { grp <- PAPER_INT; cat("using paper integrative clusters\n") }

PP <- load_PP(); WP <- load_WP()
ann <- attr(PP, "annot")
# PhosMap ID = SYMBOL_SITE
rownames(PP) <- paste0(ann$GeneSymbol, "_", ann$Residue)
wpg <- attr(WP, "annot")$Gene
idx <- match(ann$GeneSymbol, wpg)
PPn <- PP[!is.na(idx), ] - WP[idx[!is.na(idx)], colnames(PP)]
cat("raw PP sites:", nrow(PP), " protein-normalised sites:", nrow(PPn), "\n")
centre <- function(m) m - rowMeans(m)
mats <- list(raw = centre(PP), norm = centre(PPn))

to_df <- function(m) data.frame(ID = rownames(m), m, check.names = FALSE)
run_ksea <- function(m) {
  set.seed(2024)
  res <- mclapply(colnames(m), function(s) {
    set.seed(match(s, colnames(m)))
    r <- get_ksea_result_list(m[, s], rownames(m), ksrr, 0.05)
    r$ksea_activity_i
  }, mc.cores = ncores)
  ks <- sort(unique(unlist(lapply(res, names))))
  out <- sapply(res, function(v) v[ks]); rownames(out) <- ks; colnames(out) <- colnames(m)
  data.frame(Kinase = ks, out, check.names = FALSE)
}

act <- list()
for (nm in names(mats)) {
  m <- mats[[nm]]
  t0 <- Sys.time()
  act[[paste0(nm, "_mean")]] <- get_ka_by_mean_or_mlr(to_df(m), "human", log2_label = TRUE, method = "mean")
  set.seed(2024)
  act[[paste0(nm, "_mlr")]] <- get_ka_by_mean_or_mlr(to_df(m), "human", log2_label = TRUE, method = "mlr")
  act[[paste0(nm, "_ksea")]] <- run_ksea(m)
  cat(nm, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
dir.create(file.path(TAB, "kinase_activity"), showWarnings = FALSE)
for (n in names(act)) fwrite(act[[n]], file.path(TAB, "kinase_activity", paste0("06_PP_activity_", n, ".csv")))

# --- differential activity EMT vs Metabolism ---------------------------------
diff1 <- function(a, nm) {
  k <- a$Kinase; x <- as.matrix(a[, -1]); rownames(x) <- k
  x <- x[, intersect(colnames(x), colnames(PP)), drop = FALSE]
  gl <- grp[line_of(colnames(x))]
  e <- gl == "EMT"; mt <- gl == "Metabolism"
  out <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    v <- x[i, ]; ok <- is.finite(v)
    if (sum(ok & e) < 3 || sum(ok & mt) < 3 || sd(v[ok]) == 0) return(NULL)
    p <- tryCatch(t.test(v[ok & e], v[ok & mt])$p.value, error = function(err) NA)
    lfc <- if (grepl("_mean$", nm)) log2(mean(v[ok & e]) / mean(v[ok & mt])) else mean(v[ok & e]) - mean(v[ok & mt])
    data.frame(kinase = rownames(x)[i], log2FC_EMT_vs_Met = lfc, P = p)
  }))
  out$matrix <- nm
  out$assigned <- ifelse(!is.na(out$P) & out$P < 0.05 & abs(out$log2FC_EMT_vs_Met) > 1,
                         ifelse(out$log2FC_EMT_vs_Met > 0, "EMT", "Metabolism"), "")
  out
}
D <- do.call(rbind, lapply(names(act), function(n) diff1(act[[n]], n)))
fwrite(D, file.path(TAB, "06_PP_kinase_differential_all_matrices.csv"))
S <- D[D$assigned != "", ]
cnt <- as.data.frame.matrix(table(S$kinase, S$assigned))
cnt$kinase <- rownames(cnt)
if (!"EMT" %in% names(cnt)) cnt$EMT <- 0
if (!"Metabolism" %in% names(cnt)) cnt$Metabolism <- 0
paper_emt <- c("MTOR", "CDK7", "LATS1", "PIK3CB", "PKN2", "SRPK1", "ULK1", "AKT3", "MARK2", "PDGFRB", "PIK3CA",
               "PRKACB", "RIPK3", "AKT1", "CIT", "CLK1", "CLK2", "IRAK1", "MAPK1", "MAPK3", "PAK1", "ROCK1", "ROCK2",
               "RPS6KB1", "STK17A", "STK4")
paper_met <- c("STK3", "MKNK2", "STK4", "MAP2K2", "EIF2AK1", "MAP3K15")
cnt$paper_EMT <- cnt$kinase %in% paper_emt; cnt$paper_Metabolism <- cnt$kinase %in% paper_met
cnt <- cnt[order(-(cnt$EMT + cnt$Metabolism)), c("kinase", "EMT", "Metabolism", "paper_EMT", "paper_Metabolism")]
fwrite(cnt, file.path(TAB, "06_Fig8A_assignment_counts.csv"))
print(head(cnt, 40))
# paper kinases: how many matrices / direction in this run
chk <- data.frame(kinase = unique(c(paper_emt, paper_met)))
chk$paper_group <- ifelse(chk$kinase %in% paper_met & !chk$kinase %in% paper_emt, "Metabolism",
                   ifelse(chk$kinase %in% paper_emt & !chk$kinase %in% paper_met, "EMT", "both"))
chk$n_EMT <- cnt$EMT[match(chk$kinase, cnt$kinase)]; chk$n_Met <- cnt$Metabolism[match(chk$kinase, cnt$kinase)]
chk[is.na(chk)] <- 0
chk$in_activity_matrices <- sapply(chk$kinase, function(k) sum(sapply(act, function(a) k %in% a$Kinase)))
fwrite(chk, file.path(TAB, "06_Fig8A_paper_kinases_check.csv"))
print(chk)

top <- head(cnt[(cnt$EMT + cnt$Metabolism) > 0, ], 35)
save_plot("06_Fig8A_PP_kinase_assignments", {
  par(mar = c(4, 7, 3, 1))
  y <- rev(seq_len(nrow(top)))
  plot(NA, xlim = c(-6, 6), ylim = c(0.5, nrow(top) + 0.5), yaxt = "n", xlab = "times assigned (of 6 matrices)  EMT <-  -> Metabolism", ylab = "",
       main = "Differentially active kinases (P<0.05, |log2FC|>1)", xaxt = "n")
  axis(1, at = -6:6, labels = abs(-6:6))
  rect(-top$EMT, y - 0.4, 0, y + 0.4, col = "#e8312f", border = NA)
  rect(0, y - 0.4, top$Metabolism, y + 0.4, col = "#1f3c5a", border = NA)
  lab <- ifelse(top$paper_EMT | top$paper_Metabolism, paste0(top$kinase, " *"), top$kinase)
  axis(2, at = y, labels = lab, las = 1, cex.axis = 0.7)
  abline(v = 0)
  mtext("* = also reported in paper Fig. 8A", side = 3, line = 0.2, cex = 0.7)
}, w = 7, h = 8)
