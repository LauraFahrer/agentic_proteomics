# Figure and table index — run 20261105T121500Z_22794b7b

All figures are in `figures/` (PNG at 200 dpi + vector PDF). All tables are in `tables/`.

## Figures

| File | Content |
|---|---|
| `01_qc_ids_intensity` | Protein groups per run and log2 intensity distribution, coloured by Hu et al. subtype |
| `02_qc_missingness` | Detection breadth across the 14 lines; missingness vs mean abundance (intensity-dependent) |
| `03_qc_correlation` | Run-to-run Pearson correlation matrix (within-line r = 0.978, between-line r = 0.888) |
| `04_qc_cv` | Within-line coefficient of variation across 3 replicates |
| `05_concordance_with_published_matrix` | Per-protein agreement (r across lines) between quantms and the authors' Spectronaut matrix |
| `06_pca` | PCA of the 3,531 proteins quantified in all 42 runs (PC1–PC4) |
| `07_clustering_consensus` | Cell-line dendrogram + consensus matrix (k=2); labels coloured by published subtype |
| `08_hallmark_scores` | Hallmark protein scores per line, 20 most subtype-associated sets |
| `09_volcanoes` | Four differential contrasts: EMT vs Metabolism; natural EBV+ vs EBV−; EBV-like cluster vs rest; eAGS vs AGS |
| `10_eAGS_vs_naturalEBV` | log2FC of in-vitro EBV (eAGS/AGS) vs natural EBV+ module (Spearman ρ = 0.12) |
| `11_amplicon_positive_controls` | N87 17q12/ERBB2 and KatoIII 10q26/FGFR2 outlier amplicons (positive controls) |
| `12_lineage_TF_programme` | GI-lineage TF heatmap and lineage-TF score vs metabolic hallmark score (ρ = 0.74) |
| `13_EBV_module` | Natural-EBV module score per line (defined without TMK1) + pathway scores |
| `14_marker_heatmap` | Curated marker panel (epithelial, mesenchymal, gastric/intestinal, lineage TF, EBV-GC, drivers) |
| `15_ebv_peptide_search` | EBV-unique peptides per line from the human+EBV re-search of the DDA raw files |

## Tables

| File | Content |
|---|---|
| `sample_annotation.csv` | 42 runs with line, histology, EBV status, published subtype, key alterations |
| `qc_summary.json` | QC metrics |
| `protein_matrix_log2.csv.gz` | Protein group × 42 runs, log2 MaxLFQ, median-aligned |
| `protein_matrix_cellline_mean_log2.csv.gz` | Protein group × 14 lines, mean of present replicates |
| `concordance_quantms_vs_spectronaut.csv` | Per-protein correlation with the published matrix |
| `pca_loadings.csv` | PC1–PC3 loadings |
| `consensus_clustering_vs_paper_subtypes.csv` | Consensus-clustering ARI vs published subtypes (MAD × k grid) |
| `nmf_clustering_vs_paper_subtypes.csv` | NMF clustering vs published subtypes |
| `hallmark_scores_per_line.csv`, `hallmark_scores_vs_paper_subtype.csv` | Hallmark scores and subtype association |
| `diff_*.csv.gz` | Moderated-t differential results for each contrast |
| `geneset_*.csv` | Rank-based gene-set tests for each contrast |
| `presence_absence_proteins.csv` | On/off proteins per group |
| `protein_outlier_cytobands.csv` | Per-line protein-outlier cytoband hotspots (amplicon scan) |
| `correlation_with_lineageTF_score.csv.gz`, `geneset_correlation_with_lineageTF.csv` | Proteome-wide correlation with the lineage-TF score |
| `EBV_module_leave_one_out.csv` | Leave-one-out validation of the EBV module (incl. TMK1) |
| `eAGS_vs_naturalEBV_logFC.csv.gz` | Joint log2FC table for the concordance analysis |
| `marker_panel_log2.csv` | Marker-panel values |
| `ebv_peptide_hits.csv`, `ebv_search_summary.json` | EBV peptide identifications from the raw-file re-search |
| `*_summary.json` | Machine-readable summaries per analysis stage |
