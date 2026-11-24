#!/usr/bin/env Rscript
# Step 02 - dataset overview / QC (paper Fig. S1, Fig. 1B)
#  * YP site annotation (gene symbol, residue) -> processed/YP_S3.tsv
#  * replicate Pearson correlation (paper: mean 0.93 WP, 0.82 PP, 0.72 YP)
#  * GO cellular-component distribution of identified (phospho)proteins (Fig. 1B)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressPackageStartupMessages({ library(pheatmap) })

# ---------------- YP annotation ----------------
yp <- fread(file.path(PROC, "YP_S3_raw.tsv"), data.table = FALSE, check.names = FALSE)
e <- new.env(); load(file.path(SRC, "inka-intensity/inka_20220330/UniProt_data_08Jun2016.Rdata"), e)
up <- e$UniProt_data
g <- up$UniProt.Gene[match(yp$LeadAcc, up$Accession)]
g2 <- up$UniProt.Gene[match(sub("-.*", "", yp$LeadAcc), up$Accession)]
g[is.na(g)] <- g2[is.na(g)]
g[is.na(g)] <- yp$UniProtName[is.na(g)]
yp$GeneSymbol <- g
yp$Residue <- paste0(yp$Residue_aa, yp$LeadPos)
yp$id <- paste0(yp$GeneSymbol, "_", yp$Residue)
yp$id <- make.unique(yp$id, sep = "#")
ys <- grep("_r[0-9]$", colnames(yp), value = TRUE)
fwrite(yp[, c("id", "GeneSymbol", "Residue", "LeadAcc", "LeadPos", "SeqWindow", ys)],
       file.path(PROC, "YP_S3.tsv"), sep = "\t")

WP <- load_WP(); PP <- load_PP(); YP <- load_YP(); WPq <- load_WPq()

# ---------------- counts ----------------
pp_res <- substr(attr(PP, "annot")$Residue, 1, 1)
yp_res <- substr(attr(YP, "annot")$Residue, 1, 1)
counts <- list(
  WP_proteins_S3 = nrow(WP), WP_proteins_paper = 6091,
  PP_sites_S3 = nrow(PP), PP_residue_S3 = as.list(table(pp_res)),
  PP_unique_phosphopeptides_paper = 16685,
  YP_sites_S3 = nrow(YP), YP_residue_S3 = as.list(table(yp_res)), YP_pTyr_sites_paper = 8705,
  WP_proteins_quantms_processed = nrow(WPq))
jsonlite::write_json(counts, file.path(TAB, "02_counts.json"), auto_unbox = TRUE, pretty = TRUE)

# ---------------- replicate correlation (Fig. S1C) ----------------
corsum <- function(m, name) {
  cc <- cor(m, method = "pearson")
  ln <- line_of(colnames(m))
  same <- outer(ln, ln, "==") & upper.tri(cc)
  data.frame(dataset = name,
             mean_all_pairs = mean(cc[upper.tri(cc)]),
             mean_within_line_replicates = mean(cc[same]),
             min_within_line_replicates = min(cc[same]))
}
cs <- rbind(corsum(WP, "WP (Table S3)"), corsum(PP, "PP (Table S3)"), corsum(YP, "YP (Table S3)"),
            corsum(WPq, "WP (quantms, reprocessed here)"))
cs$paper_mean_pearson <- c(0.93, 0.82, 0.72, NA)
fwrite(cs, file.path(TAB, "02_replicate_correlation.csv"))
print(cs)

hm <- function(m, title) {
  cc <- cor(m)
  ann <- data.frame(line = line_of(colnames(m)), row.names = colnames(m))
  pheatmap(cc, cluster_rows = FALSE, cluster_cols = FALSE, main = title, fontsize = 5,
           annotation_col = ann, annotation_legend = FALSE, silent = TRUE)$gtable
}
save_plot("02_FigS1C_replicate_correlation", {
  gl <- list(hm(WP, "WP (Table S3)"), hm(PP, "PP (Table S3)"), hm(YP, "YP (Table S3)"), hm(WPq, "WP quantms (reprocessed)"))
  gridExtra::grid.arrange(grobs = gl, ncol = 2)
}, w = 14, h = 13)

# ---------------- GO cellular component (Fig. 1B) ----------------
cc_terms <- c(Nucleus = "GO:0005634", Nucleolus = "GO:0005730", `Extracellular matrix` = "GO:0031012",
              Ribosome = "GO:0005840", `Endoplasmic reticulum` = "GO:0005783", `Golgi apparatus` = "GO:0005794",
              Cytosol = "GO:0005829", Membrane = "GO:0016020", Cytoskeleton = "GO:0005856",
              Centrosome = "GO:0005813", `Extracellular exosome` = "GO:0070062", Lysosome = "GO:0005764",
              Peroxisome = "GO:0005777", Mitochondrion = "GO:0005739")
paper_1B <- data.frame(term = c("Nucleus", "Nucleolus", "Extracellular matrix", "Ribosome", "Endoplasmic reticulum",
                                "Golgi apparatus", "Cytosol", "Membrane", "Cytoskeleton", "Centrosome",
                                "Extracellular exosome", "Lysosome", "Peroxisome", "Mitochondrion"),
                       paper_proteins = c(2307, 493, 129, 138, 751, 581, 1786, 1785, 434, 172, 1523, 229, 88, 757),
                       paper_pproteins = c(2206, 464, 100, 99, 15, 350, 1768, 1627, 512, 231, 1172, 80, 64, 476))
wp_genes <- unique(attr(WP, "annot")$Gene)
pp_genes <- unique(c(attr(PP, "annot")$GeneSymbol, attr(YP, "annot")$GeneSymbol))
# direct (exact-term) Swiss-Prot annotations from QuickGO (bioconductor GO packages unreachable here)
res <- do.call(rbind, lapply(names(cc_terms), function(t) {
  f <- file.path(DC, "go", paste0(sub(":", "_", cc_terms[t]), "_exact.tsv"))
  sd <- unique(fread(f, sep = "\t")$SYMBOL)
  data.frame(term = t, GO = cc_terms[t],
             proteins_direct = sum(wp_genes %in% sd), pproteins_direct = sum(pp_genes %in% sd))
}))
res <- merge(res, paper_1B, by = "term")
fwrite(res, file.path(TAB, "02_Fig1B_GO_CC_distribution.csv"))
print(res)
cat("Spearman rho (paper vs direct, proteins):", cor(res$paper_proteins, res$proteins_direct, method = "spearman"), "\n")
cat("Spearman rho (paper vs direct, p-proteins):", cor(res$paper_pproteins, res$pproteins_direct, method = "spearman"), "\n")

save_plot("02_Fig1B_GO_CC_distribution", {
  par(mfrow = c(1, 2), mar = c(5, 11, 3, 1))
  o <- order(res$paper_proteins)
  barplot(t(as.matrix(res[o, c("paper_proteins", "proteins_direct")])), beside = TRUE,
          horiz = TRUE, names.arg = res$term[o], las = 1, col = c("grey30", "#3b7dd8"),
          main = "Proteins (WP)", xlab = "# proteins", cex.names = 0.8)
  legend("bottomright", c("paper", "this run (QuickGO, exact term)"),
         fill = c("grey30", "#3b7dd8"), cex = 0.7, bty = "n")
  barplot(t(as.matrix(res[o, c("paper_pproteins", "pproteins_direct")])), beside = TRUE,
          horiz = TRUE, names.arg = res$term[o], las = 1, col = c("grey30", "#d83b3b"),
          main = "Phosphoproteins (PP+YP)", xlab = "# phosphoproteins", cex.names = 0.8)
}, w = 12, h = 6)
