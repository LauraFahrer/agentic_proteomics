#!/usr/bin/env Rscript
# Step 05b - GSEA sensitivity: '-permute gene_set' instead of '-permute phenotype'.
# The paper (GSEA 4.0.3) does not state the permutation type; Table S6 reports many FDR = 0 values,
# which is typical of gene-set permutation. Re-uses the GCT/CLS files written by 05_functional.R.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
GD <- file.path(DC, "gsea")
GSEA_BIN <- Sys.getenv("GSEA_BIN", "/opt/mm/envs/repro/share/gsea-4.3.2-0/gsea-cli.sh")
MS <- file.path(SRC, "ssGSEA2.0/db/msigdb")
runs <- list(list("WP_H", file.path(MS, "h.all.v7.0.symbols.gmt"), NULL),
             list("PP_KEGG", file.path(MS, "c2.cp.kegg.v7.0.symbols.gmt"), file.path(GD, "PP.chip")),
             list("YP_KEGG", file.path(MS, "c2.cp.kegg.v7.0.symbols.gmt"), file.path(GD, "YP.chip")))
res <- list()
for (r in runs) for (src in list.dirs(GD, recursive = FALSE)[grepl(paste0("^", r[[1]], "_C[0-9]$"), basename(list.dirs(GD, recursive = FALSE)))]) {
  target <- sub(".*_", "", basename(src)); label <- paste0(basename(src), "_genesetperm")
  out <- file.path(GD, label); unlink(out, recursive = TRUE); dir.create(out)
  args <- c("GSEA", "-res", file.path(src, "data.gct"), "-cls", paste0(file.path(src, "pheno.cls"), "#", target, "_versus_REST"),
            "-gmx", r[[2]], "-collapse", ifelse(is.null(r[[3]]), "No_Collapse", "Collapse"), "-mode", "Max_probe",
            "-norm", "meandiv", "-nperm", "1000", "-permute", "gene_set", "-metric", "Signal2Noise", "-set_max", "500",
            "-set_min", "15", "-rnd_seed", "149", "-plot_top_x", "0", "-zip_report", "false", "-out", out, "-rpt_label", label)
  if (!is.null(r[[3]])) args <- c(args, "-chip", r[[3]])
  system2(GSEA_BIN, args, stdout = file.path(out, "log.txt"), stderr = file.path(out, "log.txt"))
  rd <- list.dirs(out, recursive = FALSE)[1]
  f <- list.files(rd, paste0("^gsea_report_for_", target, "_.*\\.tsv$"), full.names = TRUE)
  if (length(f)) { x <- fread(f[1])[, .(NAME, SIZE, ES, NES, `NOM p-val`, `FDR q-val`)]; x$target <- target; x$analysis <- label; res[[label]] <- x }
}
G <- rbindlist(res); fwrite(G, file.path(TAB, "05b_gsea_genesetperm.csv"))
print(G[`FDR q-val` < 0.25, .N, by = analysis])
