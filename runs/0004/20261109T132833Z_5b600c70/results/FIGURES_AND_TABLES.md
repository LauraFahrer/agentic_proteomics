# Figure and table index — run 0004 (protein-protein interaction / complexes)

## Figures
| File | Content |
|---|---|
| 01_coabundance_benchmark | Co-abundance vs known complexes: distribution + AUC, precision-vs-threshold (45x at bicor>=0.9), proteasome example |
| 02_module_landscape | 31 co-abundance modules (annotated by complex/process) across the 14 lines; none subtype-associated |
| 03_ppi_validation_subtype | A prior support of strong pairs (GO-CC/STRING/novel); B bicor vs STRING; C complex abundance vs subtype |

## Tables
| File | Content |
|---|---|
| gene_matrix_cellline_mean_log2.csv.gz | gene x line matrix (co-abundance input) |
| coabundance_bicor.parquet | gene x gene bicor matrix (NOT committed: ~148MB; regenerate with 02_coabundance.py) |
| coabundance_benchmark.json | AUC, precision, median correlations |
| coabundance_modules.csv | modules: size, cohesion, annotation, known complex, subtype assoc |
| modules_summary.json | module/strong-pair summary |
| strong_coabundance_pairs.csv.gz | pairs with bicor>=0.9 (+ in-known-complex flag) |
| strong_pairs_string_annotated.csv.gz | strong pairs + STRING combined/physical/experimental + support class |
| complex_subtype_association.csv | curated-complex mean abundance vs subtype (8 at q<0.1) |
| candidate_novel_gene_enrichment.csv | GO-CC compartment enrichment of candidate-novel-pair genes |
| hypotheses_summary.json | H1-H3 statistics |
