# Figures and tables — run 0007 (independent re-implementation)

All results were produced with self-implemented methods (`scripts/methods.py`, `scripts/networks_kinases.py`).
The paper figure is given in brackets.

## Figures (`figures/`, PNG + PDF)
| file | content |
|---|---|
| 02_FigS1C_replicate_correlation | sample Pearson heatmaps, WP/PP/YP (Table S3) and quantms WP [S1C] |
| 02_Fig1B_GO_CC_distribution | (phospho)proteins per GO-CC compartment, paper vs this run [1B] |
| 04_FigS2A_S4AE_cophenetic_grid | cophenetic vs k for 10–100 % MAD features [S2A, S4A, S4E] |
| 04_consensus_{WP,PP,YP,WPq}_* | NMF consensus matrices at the paper's settings, labelled with paper/this-run clusters [2A, 4A, 4C] |
| 04_Fig5A_integrative_binary_heatmap | binary subtype table ordered by integrative cluster [5A] |
| 05_Fig2B_4BD_GSEA_{phenotype,gene_set}_perm | −log10 FDR of paper-reported Hallmark/KEGG sets [2B, 4B, 4D] |
| 05_WGCNA_{WP,PP}_module_cluster_correlation | module eigengene vs cluster correlation [S2C, S4D] |
| 06_Fig6AB_kinase_family_share | kinome intensity shares [6A/B] |
| 06_FigS6_named_outlier_features | line-mean z-scores of the named outlier kinases [S6A–C] |
| 06_Fig8B_kinase_heatmap | paper's 40 kinases across samples [8B] |
| 06_Fig8A_PP_kinase_assignments | kinases assigned to EMT/Metabolism over the 6 activity matrices [8A] |
| 06_FigS8A_KSEA_norm_heatmap | KSEA activities, protein-normalised PP [S8A] |
| 06_Fig7A_INKA_frequency_activity | INKA top-10 frequency vs mean INKA, tyrosine kinases [7A] |
| 06_Fig7D_SRC_INKA | SRC INKA per line [7D] |
| 06_FigS5D_PTMSEA | PTM-SEA NES for the paper's signatures [S5D] |

## Tables (`tables/`)
| file | content |
|---|---|
| 01_data_summary.json, 02_* | feature counts, replicate correlation, GO-CC counts (+ Spearman vs paper) |
| 03_nmf_grid_*_mean.csv, 04_nmf_grid_all.csv | NMF grid: cophenetic, silhouette, ARI vs paper per (fraction, k) |
| 04_paper_parameters_in_grid.csv, 04_grid_best_cophenetic.csv | rank of the paper's settings; grid optimum |
| 04_subtype_assignments.csv, 04_subtypes_vs_paper.csv | per-sample clusters, silhouettes, ARI, contingency |
| 04_integrative_*, 03_integrative_clusters.csv | binary table, consensus-clustering results, k statistics |
| 04_WP_metagene_* | Kim–Park metagene proteins vs paper's 42 |
| 05_gsea_all.csv, 05_gsea_vs_paper_summary.csv, 05_paper_reported_sets.csv | GSEA (both permutation types) and recovery of paper sets |
| 05_WGCNA_* | soft-threshold fit, modules, module–cluster correlation, GO-BP |
| 05_{WP,YP}_ANOVA_GO_BP.csv, 05_PP_modt_GO_BP.csv, 05_differential_counts.json | differential features + ORA |
| 05_integrative_modt_all.csv.gz, 05_Fig5B_* | EMT vs Metabolism statistics; named proteins; pathway ORA |
| 06_kinome_counts.csv, 06_Fig6AB_kinase_family_share.csv | kinome size and family shares |
| 06_FigS6_* | outlier kinases (two definitions), named-outlier recovery, CDK4/6 test |
| 06_FigS8BC_kinome_NMF* , 06_Fig8B_* | kinome NMF; 40-kinase analysis |
| kinase_activity/06_PP_activity_*.csv | six kinase-activity matrices (raw/norm × mean/KSEA/MLR) |
| 06_PP_kinase_differential_all_matrices.csv, 06_Fig8A_* | differential activity, assignment counts, paper-kinase check |
| 06b_activity_sensitivity_* | Fig. 8A sensitivity: substring regulons, uncentred input |
| 06_INKA_scores_TK.csv, 06_Fig7A_*, 06_Fig7D_*, 06_FigS7A_top10_TK_per_sample.csv | INKA results |
| 06_PTMSEA_*.csv, 06_FigS5D_paper_signatures_check.csv | PTM-SEA NES/P per sample, EMT vs Metabolism |

Not reproduced (data unavailable): Fig. S5A (PhosphoSitePlus regulatory-site annotation), Fig. S1A/B (imputed
values not flagged in Table S3). Out of scope: Fig. 3, S3C–F (ML model), DepMap/GDSC mining (8C–E, 7E),
wet-lab experiments.
