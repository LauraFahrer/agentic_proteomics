# Figures and tables — run 0006 (reproduction of Hu et al. 2024)

Paper figure in brackets. "WPq" = quantms DIA-NN reprocessed WP; all other inputs = author-processed Table S3.

## Figures (`figures/`, PNG + PDF)
| file | content | paper |
|---|---|---|
| 02_FigS1C_replicate_correlation | Pearson sample-correlation heatmaps (WP, PP, YP, WPq) | Fig. S1C |
| 02_Fig1B_GO_CC_distribution | (phospho)proteins per GO cellular component, paper vs this run | Fig. 1B |
| 04_FigS2A_S4AE_cophenetic_grid | cophenetic coefficient vs k for top 10/20/30/50/100 % MAD | Fig. S2A, S4A, S4E |
| 04_consensus_{WP,PP,YP,WPq}_* | NMF consensus matrices at the paper's parameters, annotated with the paper's clusters | Fig. 2A, 4A, 4C |
| 04_silhouette_* | silhouette widths (consensus) | Fig. S2B, S4B, S4F |
| 04_Fig5A_integrative_binary_heatmap | binary subtype table ordered by ConsensusClusterPlus integrative cluster | Fig. 5A |
| 05_Fig2B_4BD_GSEA_{phenotype,gene_set}_perm | −log10 FDR for the paper's Hallmark (WP) and KEGG (PP, YP) sets | Fig. 2B, 4B, 4D |
| 05_WGCNA_{WP,PP}_module_cluster_correlation | module eigengene vs NMF cluster correlation | Fig. S2C, S4D |
| 08_Fig6AB_kinase_family_share | kinome intensity share by kinase group (WP, PP) and TK family (YP) | Fig. 6A, 6B |
| 08_FigS6_named_outlier_features | line-mean z-scores for ERBB2, CDK12, FGFR2, CSNK1D/E, SRC … | Fig. S6A–C |
| 08_FigS8B_kinome_NMF_consensus | NMF (k = 2) on WP kinases | Fig. S8B |
| 08_Fig8B_kinase_heatmap | paper's 40 kinases across samples, ordered by integrative subtype | Fig. 8B |

## Tables (`tables/`)
| file | content |
|---|---|
| 01_data_summary.json, 02_counts.json | feature counts (Table S3, quantms) vs paper |
| 02_replicate_correlation.csv | mean within-line / all-pair Pearson, vs paper 0.93/0.82/0.72 |
| 02_Fig1B_GO_CC_distribution.csv | GO CC counts, paper vs QuickGO |
| 03_nmf_grid_{WP,PP,YP,WPq}_mean.csv, 04_nmf_grid_all.csv | full NMF grid: frac × k → cophenetic, silhouette, dispersion, assignments |
| 04_nmf_grid_median_ad.csv | WP grid using median absolute deviation (Fig. S2 legend wording) |
| 04_paper_parameters_in_grid.csv | where the paper's chosen settings rank in the grid |
| 04_grid_best_cophenetic*.csv | grid optimum (all k; k > 2) |
| 04_subtypes_vs_paper.csv | ARI, contingency vs paper for WP/PP/YP/WPq |
| 04_subtype_assignments.csv | per-sample clusters, silhouettes, paper clusters |
| 04_integrative_* , 03_integrative_clusters.csv | binary table, ConsensusClusterPlus assignments, k statistics |
| 04_WP_metagene_features.csv | NMF extractFeatures metagene proteins vs paper's 42 |
| 05_gsea_all.csv, 05b_gsea_genesetperm.csv | GSEA results (phenotype / gene-set permutation) |
| 05c_gsea_vs_paper_summary.csv | recovery of paper-reported gene sets |
| 05_paper_TableS6_hallmarks.csv, 05c_paper_reported_sets.csv | paper's reported sets (Table S6, Fig. 4B/4D) |
| 05_WGCNA_*_modules.csv, *_module_GO.csv, 05_WGCNA_summary.json | modules, module–cluster correlation, GO-BP |
| 05_{WP,YP}_ANOVA*.csv*, 05_PP_limma_GO_BP.csv, 05_differential_counts.json | differential features and GO-BP ORA |
| 05_integrative_limma_all.csv.gz | EMT vs Metabolism limma for WP, PP, YP |
| 05_Fig5B_named_proteins_check.csv, 05_Fig5B_integrative_pathway_ORA.csv | Fig. 5B proteins and pathway ORA |
| 08_kinome_counts.csv, 08_Fig6AB_kinase_family_share.csv | kinome size and family shares vs paper |
| 08_FigS6_outlier_kinases.csv, 08_FigS6_named_outliers_recovery.csv, 08_FigS6_named_features_line_zscores.csv | outlier kinases (two definitions) |
| 08_FigS6G_CDK46_MKN28_vs_YCCEL1.csv | CDK4/6 comparison |
| 08_FigS8BC_kinome_NMF*.{csv,json} | kinome NMF |
| 08_Fig8B_kinase_limma.csv, 08_Fig8B_paper40_check.csv, 08_Fig8B_paper40_hclust.csv | differential kinases |

Not produced (tools not executed, see REPORT.md section 7): Fig. 7A/7D & S7A (InKA), Fig. 8A & S8A (PhosMap six-matrix),
Fig. S5D (PTM-SEA), Fig. S5A (PhosphoSitePlus). Excluded by scope: Fig. 3, S3C–F (ML model), 6C/D, 7B/C/E,
8C–F, S6D–J (experiments / public-database mining).
