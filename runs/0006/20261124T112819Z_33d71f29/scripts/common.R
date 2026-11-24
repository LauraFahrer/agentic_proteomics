# Shared helpers for the run-0006 reproduction of Hu et al. 2024 (Mol Oncol 18:1437)
suppressPackageStartupMessages({
  library(data.table)
})

RUN <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
DC  <- file.path(RUN, "data_cache")
PROC <- file.path(DC, "processed")
FIG <- file.path(RUN, "results", "figures")
TAB <- file.path(RUN, "results", "tables")
SRC <- Sys.getenv("REPRO_SRC", "/opt/src")   # cloned GitHub tools (see 00_setup.sh)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

LINES <- c("AGS", "AZ521", "eAGS", "KatoIII", "MKN28", "N87", "NCC24", "Ocum1",
           "SNU1750", "SNU601", "SNU719", "TMK1", "YCC3", "YCCEL1")

# ---- published subtype assignments (read from the paper's figures) ----
# WP (Fig. 2A): cluster 1 = Cell prolif/EMT, 2 = EMT, 3 = Metabolism
PAPER_WP <- c(AGS = 1, eAGS = 1, MKN28 = 1, N87 = 1,
              AZ521 = 2, SNU1750 = 2, TMK1 = 2,
              KatoIII = 3, NCC24 = 3, Ocum1 = 3, SNU601 = 3, SNU719 = 3, YCC3 = 3, YCCEL1 = 3)
# PP (Fig. 4A): 1 = EMT-like, 2 = metabolism-like (identical to integrative S1/S2, Fig. 5A)
PAPER_PP <- c(AGS = 1, AZ521 = 1, eAGS = 1, MKN28 = 1, N87 = 1, SNU1750 = 1, TMK1 = 1,
              KatoIII = 2, NCC24 = 2, Ocum1 = 2, SNU601 = 2, SNU719 = 2, YCC3 = 2, YCCEL1 = 2)
PAPER_INT <- ifelse(PAPER_PP == 1, "EMT", "Metabolism")
# YP (Fig. 4C / S4g; Ocum1 replicates split, handled at sample level in paper_label())
PAPER_YP <- c(AGS = 1, AZ521 = 1, MKN28 = 1, SNU1750 = 1, TMK1 = 1, YCC3 = 1,
              KatoIII = 2, N87 = 2,
              eAGS = 3, NCC24 = 3, Ocum1 = 3, SNU601 = 3, SNU719 = 3, YCCEL1 = 3)

line_of <- function(s) sub("_r[0-9]+$", "", s)
# sample-level paper labels; Fig. 4C places Ocum1_r1 in YP cluster 1 and Ocum1_r2 in YP cluster 3
paper_label <- function(d, samples) {
  lab <- switch(d, WP = PAPER_WP, WPq = PAPER_WP, PP = PAPER_PP, YP = PAPER_YP)[line_of(samples)]
  if (d == "YP") lab[samples == "Ocum1_r1"] <- 1
  unname(lab)
}

read_mat <- function(file, idcol = "id") {
  d <- fread(file.path(PROC, file), sep = "\t", data.table = FALSE, check.names = FALSE)
  sc <- grep("_r[0-9]$", colnames(d), value = TRUE)
  m <- as.matrix(d[, sc]); rownames(m) <- d[[idcol]]
  attr(m, "annot") <- d[, setdiff(colnames(d), sc), drop = FALSE]
  m
}

load_WP <- function() read_mat("WP_S3.tsv")
load_PP <- function() read_mat("PP_S3.tsv")
load_WPq <- function() read_mat("WP_quantms_processed.tsv")
load_YP <- function() read_mat("YP_S3.tsv")   # written by 02_qc.R (gene-annotated)

# Mean absolute deviation per row (Methods: "mean absolute deviation (MAD)")
mean_ad <- function(m) rowMeans(abs(m - rowMeans(m)))
# Median absolute deviation (Fig. S2 legend: "median absolute deviations")
median_ad <- function(m) apply(m, 1, mad)

top_frac <- function(m, frac, fun = mean_ad) {
  v <- fun(m); n <- ceiling(frac * nrow(m))
  m[order(v, decreasing = TRUE)[seq_len(n)], , drop = FALSE]
}

save_plot <- function(name, expr, w = 7, h = 6) {
  ex <- substitute(expr); env <- parent.frame()
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h); eval(ex, env); invisible(dev.off())
  png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in", res = 150); eval(ex, env); invisible(dev.off())
}

# adjusted Rand index (no extra dependency)
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  s1 <- sum(choose(tab, 2)); s2 <- sum(choose(rowSums(tab), 2)); s3 <- sum(choose(colSums(tab), 2))
  e <- s2 * s3 / choose(n, 2); (s1 - e) / (0.5 * (s2 + s3) - e)
}
