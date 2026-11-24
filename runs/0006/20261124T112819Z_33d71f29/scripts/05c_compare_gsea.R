#!/usr/bin/env Rscript
# Step 05c - compare GSEA results (phenotype and gene-set permutation) with the paper
# (Table S6 / Fig. 2B for WP Hallmarks; Fig. 4B, 4D for PP/YP KEGG)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressPackageStartupMessages(library(pheatmap))
A <- fread(file.path(TAB, "04_subtype_assignments.csv"))
# map reproduced cluster -> paper cluster (majority)
cmap <- A[, .N, by = .(dataset, cluster, paper_cluster)][order(-N)][, .SD[1], by = .(dataset, cluster)]
pname <- list(WP = c("1" = "Cell prolif/EMT", "2" = "EMT", "3" = "Metabolism"), PP = c("1" = "PP1 (EMT)", "2" = "PP2 (Metabolism)"),
              YP = c("1" = "YP1", "2" = "YP2", "3" = "YP3"))
G <- rbind(fread(file.path(TAB, "05_gsea_all.csv"))[!grepl("paperlab", analysis)][, perm := "phenotype"],
           fread(file.path(TAB, "05b_gsea_genesetperm.csv"))[, perm := "gene_set"])
G[, dataset := sub("_.*", "", analysis)]
G[, cluster := as.integer(sub("C", "", target))]
G <- merge(G, cmap[, .(dataset, cluster, paper_cluster)], by = c("dataset", "cluster"))
G[, paper_group := mapply(function(d, p) pname[[d]][as.character(p)], dataset, paper_cluster)]
G[, NAME := toupper(NAME)]
fwrite(G, file.path(TAB, "05c_gsea_with_paper_mapping.csv"))

paper <- rbind(
  fread(file.path(TAB, "05_paper_TableS6_hallmarks.csv"))[, .(dataset = "WP", set, group)],
  data.table(dataset = "PP", group = "PP2 (Metabolism)",
             set = paste0("KEGG_", c("NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "CALCIUM_SIGNALING_PATHWAY", "FC_EPSILON_RI_SIGNALING_PATHWAY",
                   "CHEMOKINE_SIGNALING_PATHWAY", "FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS", "GAP_JUNCTION", "VEGF_SIGNALING_PATHWAY",
                   "LONG_TERM_POTENTIATION", "REGULATION_OF_ACTIN_CYTOSKELETON", "GNRH_SIGNALING_PATHWAY", "ENDOCYTOSIS"))),
  data.table(dataset = "YP", group = "YP3",
             set = paste0("KEGG_", c("JAK_STAT_SIGNALING_PATHWAY", "CELL_ADHESION_MOLECULES_CAMS", "GNRH_SIGNALING_PATHWAY", "FOCAL_ADHESION",
                   "FC_EPSILON_RI_SIGNALING_PATHWAY", "ERBB_SIGNALING_PATHWAY", "T_CELL_RECEPTOR_SIGNALING_PATHWAY", "INSULIN_SIGNALING_PATHWAY",
                   "NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "CHEMOKINE_SIGNALING_PATHWAY"))),
  data.table(dataset = "YP", group = "YP1", set = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"))
fwrite(paper, file.path(TAB, "05c_paper_reported_sets.csv"))

cmp <- rbindlist(lapply(split(paper, by = c("dataset", "group")), function(p) {
  rbindlist(lapply(c("phenotype", "gene_set"), function(pm) {
    g <- G[dataset == p$dataset[1] & paper_group == p$group[1] & perm == pm]
    sig <- g[`FDR q-val` < 0.25 & NES > 0]$NAME
    data.table(dataset = p$dataset[1], paper_group = p$group[1], permutation = pm, n_paper = nrow(p),
               n_this_run_FDR25 = length(sig), n_overlap = sum(p$set %in% sig),
               n_paper_with_positive_NES = sum(p$set %in% g[NES > 0]$NAME),
               median_rank_of_paper_sets = median(match(p$set, g[order(-NES)]$NAME), na.rm = TRUE),
               n_sets_tested = nrow(g))
  }))
}))
fwrite(cmp, file.path(TAB, "05c_gsea_vs_paper_summary.csv")); print(cmp)

hm <- function(d, pm_, sets, title) {
  g <- G[dataset == d & perm == pm_ & NAME %in% sets]
  g <- rbind(g, CJ(NAME = sets, paper_group = unique(G[dataset == d]$paper_group))[, `:=`(`FDR q-val` = 1, NES = 0)], fill = TRUE)
  g <- g[!duplicated(g[, .(NAME, paper_group)])]
  w <- dcast(g, NAME ~ paper_group, value.var = "FDR q-val")
  nes <- dcast(g, NAME ~ paper_group, value.var = "NES")
  m <- -log10(pmax(as.matrix(w[, -1]), 1e-3)); rownames(m) <- w$NAME
  m[is.na(m) | as.matrix(nes[, -1]) <= 0] <- 0
  pheatmap(m[intersect(sets, rownames(m)), , drop = FALSE], cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("white", "#e8312f"))(30), breaks = seq(0, 3, length.out = 31),
           main = title, fontsize_row = 7, display_numbers = round(m[intersect(sets, rownames(m)), , drop = FALSE], 1), silent = TRUE)$gtable
}
fig2b <- paste0("HALLMARK_", c("COMPLEMENT", "GLYCOLYSIS", "UV_RESPONSE_UP", "ALLOGRAFT_REJECTION", "P53_PATHWAY", "INTERFERON_ALPHA_RESPONSE",
         "INTERFERON_GAMMA_RESPONSE", "KRAS_SIGNALING_UP", "PEROXISOME", "XENOBIOTIC_METABOLISM", "BILE_ACID_METABOLISM", "ADIPOGENESIS",
         "FATTY_ACID_METABOLISM", "OXIDATIVE_PHOSPHORYLATION", "APOPTOSIS", "PI3K_AKT_MTOR_SIGNALING", "MYC_TARGETS_V1", "UV_RESPONSE_DN",
         "MITOTIC_SPINDLE", "TGF_BETA_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "HEDGEHOG_SIGNALING", "E2F_TARGETS", "G2M_CHECKPOINT",
         "HYPOXIA", "EPITHELIAL_MESENCHYMAL_TRANSITION", "WNT_BETA_CATENIN_SIGNALING"))
for (pm in c("phenotype", "gene_set")) {
  save_plot(paste0("05_Fig2B_4BD_GSEA_", pm, "_perm"), {
    gridExtra::grid.arrange(hm("WP", pm, fig2b, paste("WP Hallmark (Fig. 2B), -log10 FDR,", pm, "perm.")),
                            hm("PP", pm, paper[dataset == "PP"]$set, "PP KEGG (Fig. 4B)"),
                            hm("YP", pm, paper[dataset == "YP"]$set, "YP KEGG (Fig. 4D)"), ncol = 3, widths = c(1.2, 1, 1))
  }, w = 20, h = 9)
}
