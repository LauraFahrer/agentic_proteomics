# Evaluation

## Methods from Paper

- 3 datasets:  Whole Proteome (WP), Phospho sites (PP), tyrosine-phospho sites (YP)
- subtype classification:
    - all 3 datasets, first seperate then integrative clustering
    - selection of variable (phospho-) proteins with MAD (Median Absolute Deviation)
    - clustering with NMF (Non-negative Matrix Factorization)
    - cluster analysis: Cophenetic Correlation (Stabilität der Clusterstruktur über Wiederholungen), Consensus Matrix (Stabilität der Gruppierung der Zelllinien), Silhouette Scores (Ähnlichkeit innerhalb und Diversität zwischen Gruppen)
    - ConsensusClusterPlus (integration of all 3 datsets for evaluation of cluster consistency)
- functional annotation:
    - subcellular location: Gene Ontology Cellular Components
    - identification of protein modules (differentially co-expressed proteins, Weighted Gene Co-expression Network Analysis, WGCNA package)
    - GSEA (Gene Set Enrichment Analysis)
    - Gene Ontology Analysis
    - PTM-SEA (Post-Translational Modification Signature Enrichment Analysis) for PP dataset
- signaling pathway analysis through kinase activity analysis:
    - InKA (Integrated Kinase Activity Analysis) on YP dataset
    - PP raw data and normalised, 3 methods: Mean Fold Difference, KSEA (Kinase-Set Enrichment Analysis), MLR (Multilinear Regression) -> 6 activity matrices -> identifying kinases with activity differences between subtypes (p < 0.05, Unterschied in Aktivität um Faktor 2)
- ML-based subtype-prediction on real clinical patient data: using 42 proteins as features, 3 labels(EMT, Cell proliferation/EMT, Metabolism)

## Results from Paper

- identification of 6091 proteins, 16 685 unique phosphopeptides and 8705 phosphorylated Tyr sites across the panel
- 3 distinct clusters (WP data) -> identification of 2 subtypes when considering protein modules and functional analysis (EMT/proliferation vs metabolism)
- differential kinase activity:
    - increased SFK activity in many cell lines -> growth inhibition through eCF506 inhibitor
    - cell line dependend increased activity: ERBB2 in N87 (known), FGFR2 in KatoIII (known), CSNK1D/E in MKN28 (new)
    - subtype specific increased activity: mTOR -> EMT, MAPK (MAP2K2) -> Metabolism
    - EMT: sensitivity against mTOR inhibitors (Rapamycin), cells may activate AKT through an alternative signaling pathway that does not involve PI3K
    - Metabolism: sensitivity against MEK inhibitors (Trametinib)

## Results from Agentic Research

- no independence analysis but using paper as reference and rebuild the workflow from the paper (e.g. used same MAD thresholds as in paper)
- identification of 8152 protein groups at 1% FDR (quantms provenance: 9080 proteins; the paper’s own pipeline reported 6091), phospho layers from the paper’s supplement: 8452 PP sites + 4780 pY sites
- reproducing EMT/Metabolism subtypes for 13 of 14 celllines (except TMK1 -> agent labeled it as Metabolism)
- EMT not mesenchymal at protein level (vimentin, but specificity is controversial), differentially expressed transcription factors driving oxidative metabolism (HNF1A/HNF4A/FOXA2)
- most stable structure with four groups (one (naturally) EBV-driven -> distinct proteome module)
- differential kinase activity:
    - cell line dependend increased activity: ERBB2 in N87 (known), FGFR2 in KatoIII (known), CRKL in MKN28
    - similar subtype specific increased activity as reported in the paper
    - robust YP signaling but not for S/T after normalization
    - differential expression of RTKs (receptor tyrosine kinase) between subtypes (not reported in paper)
- mapping of known complexes and subtype-specific differences: (OXPHOS/TCA complexes up in Metabolism; transcription/ubiquitin complexes up in EMT)

## Comparison and Evaluation

| Step | Paper | Agent |
| --- | --- | --- |
| Data processing | - based on tools and R-packages<br>- e.g. MaxQuant, Spectronaut, proBatch, ... | - WP: self-written scripts using tools/packages like duckdb/numpy/...<br>- PP/YP: using processed data from paper supplementary |
| Subtype classification | - pipeline: MAD, NMF, ConsensusClusterPlus, cophenetic, silhouette R | - self-written: MAD, investigation of cluster stability (consensus matrix, hierarchical consensus clustering, bootstrap adjusted rand index)<br>- scikit-learn NMF, PCA |
| Functional annotation | - using established tools: WGCNA, GSEA, limma, PTM-SEA | - Gene Ontology MSigDB, KEGG, Reactome<br>- self-written: identification of co-expressed proteins, enrichment analysis, statistical anlysis
| Kinase-activity analysis | - INKA, KSEA, MLR, mean-fold-difference | - phosX (Kinase Library Ser/Thr Johnson 2023 + Tyr Yaron-Barir 2024 PSSMs, open-source python tool) |

- paper mainly used established tools, software applications, packages and public databases like GSEA, Gene Ontology Resource and PTMsigDB
- agentic research also has used to a large extend self-written scripts
- different processing of WP dataset -> different findings<br>using processed PP and YP data from paper -> similar findings
- accessibility issues in minor of cases: INKA, PTM-SEA on GenePattern, the GSEA desktop app’s MSigDB fetch, PhosphoSitePlus
- default bias toward coding the statistic directly (stated from Claude): writing own scripts, when method is simple enough, use tools in case of hitting a wall (e.g. phosX), avoiding installations and keeping every step auditable
