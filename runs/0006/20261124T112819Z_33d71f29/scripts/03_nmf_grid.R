#!/usr/bin/env Rscript
# Step 03 - NMF subtype identification grid (paper Methods 2.10; Fig. S2A, S4A, S4E)
# Usage: Rscript 03_nmf_grid.R <WP|PP|YP|WPq> [mean|median]
#  * features = top 10/20/30/50/100 % by absolute deviation (Methods: mean absolute deviation;
#    Fig. S2 legend says median absolute deviation -> both variants available)
#  * NMF R package (Gaujoux & Seoighe), default algorithm 'brunet', nrun = 100, k = 2..5
#  * stability = cophenetic correlation coefficient of the consensus matrix
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressPackageStartupMessages(library(NMF))
args <- commandArgs(TRUE)
ds <- args[1]; madtype <- ifelse(length(args) > 1, args[2], "mean")
m <- switch(ds, WP = load_WP(), PP = load_PP(), YP = load_YP(), WPq = load_WPq())
stopifnot(min(m) >= 0)
fun <- if (madtype == "mean") mean_ad else median_ad
fracs <- c(0.1, 0.2, 0.3, 0.5, 1.0); ks <- 2:5
out <- file.path(DC, "nmf"); dir.create(out, showWarnings = FALSE)
res <- list()
for (fr in fracs) {
  x <- top_frac(m, fr, fun)
  for (k in ks) {
    t0 <- Sys.time()
    fit <- nmf(x, rank = k, method = "brunet", nrun = 100, seed = 123456, .opt = "-p")
    coph <- cophcor(fit)
    sil <- silhouette(fit, what = "consensus")
    cl <- predict(fit, what = "consensus")
    saveRDS(list(fit = fit, features = rownames(x)), file.path(out, sprintf("%s_%s_%03d_k%d.rds", ds, madtype, round(fr * 100), k)))
    res[[length(res) + 1]] <- data.frame(dataset = ds, mad = madtype, frac = fr, n_features = nrow(x), k = k,
                                         cophenetic = coph, mean_silhouette = mean(sil[, "sil_width"]),
                                         dispersion = dispersion(fit),
                                         clusters = paste(names(cl), cl, sep = ":", collapse = ";"))
    cat(ds, fr, k, "coph", round(coph, 4), "sil", round(mean(sil[, "sil_width"]), 3),
        "t", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
  }
}
res <- do.call(rbind, res)
fwrite(res, file.path(TAB, sprintf("03_nmf_grid_%s_%s.csv", ds, madtype)))
