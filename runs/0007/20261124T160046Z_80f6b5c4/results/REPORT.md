# Independent reproduction of Hu et al. 2024 (Mol Oncol 18:1437) — PXD048538, quantms e60b9e9aa4c8

Run `0007` · prompt `prompt_0007.yaml` · scope: data processing, subtype classification, functional annotation,
kinase-activity analysis (the machine-learning model and public-database mining were excluded as instructed).

**Approach.** Same input data as the paper; every specialised analytical method is **re-implemented** in Python
(`scripts/methods.py`, `scripts/networks_kinases.py`). Only general-purpose libraries are used:

- numpy for linear algebra and random numbers;
- scipy for distributions, special functions, hierarchical-linkage primitives, `f_oneway` and `ttest_ind`;
- pandas and matplotlib.

No specialised package or web application (Spectronaut, NMF, ConsensusClusterPlus, WGCNA, limma, GSEA, PTM-SEA,
PhosMap, InKA, PANTHER) was executed. Run 0006 (same repository) executed several of the original tools and serves
as an implementation cross-check.

## 1. Classification of what was used

| category | items |
|---|---|
| **Self-implemented** (this run) | quantile normalisation · down-shift imputation · mean/median absolute deviation · Brunet KL-NMF with the R-NMF connectivity stopping rule · consensus matrix, cophenetic coefficient, silhouette · Kim–Park metagene selection · Monti/ConsensusClusterPlus resampling consensus clustering · limma-type moderated t (lmFit/eBayes incl. fitFDist) · BH FDR · GSEA (Signal2Noise with GSEA σ-floor, weighted KS, phenotype and gene-set permutation, NES, GSEA FDR) · GSEA Max_probe collapse · hypergeometric ORA · ssGSEA/PTM-SEA (rank normalisation, weight 0.75, area ES, directional signatures, permutation NES/P) · WGCNA-type network (scale-free fit, TOM, eigengenes, eigengene merging) with an own simplified tree cut · kinase activity by mean / PhosMap-style KSEA / ridge MLR with CV · INKA score (kinase-centric × substrate-centric) · outlier definitions · all figures |
| **General-purpose library functions** | numpy `linalg.svd/inv`, random generators; scipy `hierarchy.linkage/fcluster/cophenet`, `stats.f_oneway`, `stats.ttest_ind`, `stats.hypergeom`, `stats.t`, `special.digamma/polygamma`, `spearmanr`, `pearsonr`; pandas; matplotlib |
| **External data resources** (databases; unavoidable) | authors' processed Table S3 (WP/PP/YP); quantms DIA-NN 2.5.1 output (WP); MSigDB v7.0 Hallmark/KEGG/GO-BP gene sets; PTMsigDB v1.9.1 (flanking); PhosMap human kinase–substrate table; PhosphoSitePlus + NetworKIN (Jul 2016, KSEA-app release); InKA activation-loop table and HGNC/Manning kinase list (2016); UniProt Swiss-Prot sequences and kinase family annotations; QuickGO GO-CC annotations |
| **External analytical tools** | **none executed.** R limma was used once only to *validate* the self-implemented moderated t on simulated data: max \|Δt\| = 0.001, max \|ΔP\| = 1.3e-4 |
| **Upstream (not reproducible)** | Spectronaut 13 (commercial) and MaxQuant 1.6 searches. Downstream analyses therefore start from Table S3. The quantms DIA-NN reprocessing of WP is used as an independent-processing check (the paper's quantile normalisation + imputation re-applied with own code) |

## 2. Results versus the paper

| Paper item | This run (independent implementation) | Verdict |
|---|---|---|
| Replicate Pearson r 0.93 / 0.82 / 0.72 (Fig. S1C) | **0.931 / 0.826 / 0.721** (within-line replicate pairs) | reproduced |
| GO-CC landscape (Fig. 1B) | counts correlate with paper: ρ = 0.99 (proteins), 0.91 (phosphoproteins); annotation release differs | approximately reproduced |
| WP NMF, top-10 % MAD (568), k = 3 (Fig. 2A) | **identical clusters (ARI 1.0)**: {AGS, eAGS, MKN28, N87} / {AZ521, SNU1750, TMK1} / 7 metabolism lines; cophenetic 0.954, mean silhouette 0.73 (paper 0.87) | reproduced (silhouette lower) |
| PP NMF, top-20 %, k = 2 (Fig. 4A) | 39/42 samples match; **AGS (3 replicates) falls into the Metabolism cluster** (ARI 0.73; AGS silhouettes 0.35, −0.04, 0.05). Top-10 %, 30 % and 100 % give the paper's split exactly. Run 0006 (R `NMF`) also gave the exact split at 20 % | partially reproduced — AGS is a borderline sample |
| YP NMF, 100 %, k = 3 (Fig. 4C) | **identical (ARI 1.0)**, including the split Ocum1 replicates | reproduced |
| Integrative ConsensusClusterPlus (Fig. 5A) | **identical EMT / Metabolism split (ARI 1.0)**, both replicates concordant for all 14 lines; k = 2 consensus 0.996 | reproduced |
| Parameter choice by "most stable decomposition" (Fig. S2A, S4A/E) | the paper's settings rank 10th (WP), 4th (PP) and 11th (YP) of 20 by cophenetic coefficient. The maximum is at 10 %, k = 2 for every dataset | **not reproduced** (selection rule underdetermined) |
| 42 NMF metagene proteins (Fig. S3A) | Kim–Park selection yields 39 proteins, of which 23 are in the paper's 42; 40/42 are in the top-10 % feature set | partially |
| Independent WP processing (quantms) | Metabolism cluster identical; EMT side splits {AGS, eAGS} vs {AZ521, MKN28, N87, SNU1750, TMK1} (ARI 0.78) | EMT/Metabolism dichotomy robust; finer split processing-dependent |
| GSEA Hallmark per WP cluster (Fig. 2B, Table S6) | same leading terms; Metabolism (OXPHOS, bile-acid, fatty-acid, xenobiotic, adipogenesis, peroxisome) 4/11 Table S6 sets at FDR < 0.25 with phenotype permutation, **10/11 with gene-set permutation**. Cell-prolif/EMT 0 → 9/27, EMT 0 → 3/18 | direction reproduced; significance depends on the unreported permutation type |
| GSEA KEGG PP (Fig. 4B) / YP (Fig. 4D) | 10/11 PP and 10/10 YP paper sets have positive NES in the matching cluster; 0–2 (PP) and 0–2 (YP) reach FDR < 0.25 | direction reproduced, significance weak |
| WGCNA 12 WP / 5 PP modules (Fig. S2C, S4D) | 25 WP / 27 PP modules (own tree cut; parameters not reported in the paper; R WGCNA in run 0006 gave 32 / 44). Strongest modules track the clusters (e.g. WP module 1 r = −0.93 with Cell-prolif/EMT, module 3 r = −0.95 with EMT) | **not reproduced** (method underspecified) |
| ANOVA/limma + GO (Fig. S2D, S4C/G) | WP ANOVA 2,239 proteins (BH < 0.05); Metabolism: oxidation-reduction, organic-acid / lipid metabolism. YP ANOVA 754 sites: integrin, cell–cell adhesion, ephrin, RTK signalling; Fc-receptor signalling. PP (moderated t, 707 sites): no GO-BP term at BH < 0.05 | partially |
| Fig. 5B named proteins | **all 16** proteins named at protein (or both) level significant in the paper's direction; EMT phosphosites on TGFB1I1, CTNNB1, YAP1, ABL2, APC, PHLDB2 significant up; named metabolism phosphosites absent from Table S3 | reproduced where testable |
| PTM-SEA (Fig. S5D: CDK1/2/4/6 enriched in EMT) | **CDK6 (P = 1e-4, BH 0.002) and CDK4 (P = 0.008) significantly higher in EMT**; CDK2 (P = 0.036) and CDK1 same direction. PKCA, HIPK2 and the EGFR1 pathway higher in Metabolism. Src, CCR7, B-cell receptor and P38B signatures have < 10 sites | reproduced for the CDK signal |
| Kinome size (Fig. 6) | 178 / 158 / 155 kinases (WP / PP / YP) vs 198 / 193 / 140 (different kinase list) | approximate |
| YP TK families (Fig. 6B) | FGFR 30.5 % (33.0), EPH 25.6 (32.9), SRC 11.9 (14.7), EGFR 6.4 (8.8), JAK 3.5 (3.8), MET 2.8 (3.4), INSR 0.9 (1.2), ABL 0.7 (0.7) | reproduced closely |
| WP/PP kinase-group shares (Fig. 6A) | do not match (e.g. Atypical 22 % in paper WP) — family source not reported | not reproduced |
| Outlier kinases (Fig. S6A–C) | ERBB2 and CDK12 (N87), FGFR2 PP + YP (KatoIII only), CSNK1D S384 / CSNK1E S389 (MKN28) all recovered; total 81 or 265 outlier kinases by two operationalisations vs 124 | named findings reproduced; count definition-dependent |
| CDK4 / CDK6 YCCEL1 vs MKN28 (Fig. S6G) | Student t: **P = 1.0e-4 and 0.0214** (paper 1e-4, 0.0214) | reproduced exactly |
| Kinome NMF (Fig. S8B/C) | 179 kinases, k = 2: 41/42 samples follow the integrative split (MKN28_r3 off; ARI 0.90); silhouette 0.74 (paper 1.0) | largely reproduced |
| 40 differential kinases (Fig. 8B) | all 37 testable paper kinases (COQ8A, TGM2, MAP3K20 absent from the kinase list) are in this run's top-40; clustering on them separates the subtypes perfectly (ARI 1.0) | reproduced |
| INKA on YP: first tier SRC, ABL1, EGFR, EPHA2, PTK2, MET (Fig. 7A) | most frequent top-10 TKs (of 28 samples): EGFR 27, EPHA2 27, LYN 27, ERBB2 25, LCK 22, **SRC 21**, ABL1 19, TYK2 18, MET 15, … PTK2 12 → 5 of the 6 paper kinases in our top 9 | largely reproduced (LYN/ERBB2/LCK additionally prominent) |
| SRC activity across lines (Fig. 7D) | rank correlation with the paper's ordering ρ = 0.55 (P = 0.043); KatoIII and N87 are the two lowest in both (paper: KatoIII < N87; here N87 < KatoIII) | partially reproduced |
| PP 6-matrix kinase activity: MTOR, CDK7, LATS1 (EMT); STK3/4, MKNK2, MAP2K2, EIF2AK1 (Metabolism) (Fig. 8A) | none of these assigned in any matrix; only CLK2 (2×) and MAPK1 (1×) of the paper's 31 kinases recovered; MARK2 assigned to Metabolism (6×), the opposite of the paper. Most assignments are CDKs/MAPKs (EMT) and CDK19, PKMYT1, SRC (Metabolism). Sensitivity (`06b`): PhosMap-style substring regulon matching and uncentred input do not recover them either | **not reproduced** |

## 3. Interpretation

1. **Robust, independently reproduced:**
   - the WP and YP NMF subtypes and the integrative EMT/Metabolism classification;
   - the 40-kinase expression signature;
   - the subtype marker proteins of Fig. 5B;
   - CDK-signature (PTM-SEA) enrichment in the EMT subgroup;
   - outlier kinases of individual lines (ERBB2/CDK12 in N87, FGFR2 in KatoIII, CSNK1D/E in MKN28);
   - broad SRC-family/EGFR/EPHA2 tyrosine-kinase activity across the panel (INKA).

   These results were obtained with code written from the method descriptions, so they are not artefacts of a
   particular software package.
2. **Implementation-sensitive:** the PP 2-cluster solution at 20 % features (AGS flips between implementations),
   the silhouette values (all lower than reported), GSEA significance (permutation type), and the number of WGCNA
   modules.
3. **Not reproduced:** the subtype-specific kinase-activity calls of Fig. 8A. This is the computational basis of
   the paper's mTOR (EMT) vs MEK (Metabolism) therapeutic hypothesis. With the same data, the same kinase–substrate
   resource and the stated threshold, none of the highlighted kinases differ.

   The most likely causes are undocumented details:
   - input scaling and the replicate level used;
   - possibly line-level rather than sample-level testing;
   - PhosMap version-specific behaviour.

   Two of these were tested and excluded: substring matching and uncentred input. **The available computational
   evidence is insufficient to support Fig. 8A.** The drug-response (trametinib/rapamycin) and western-blot
   validations in the paper are experimental and outside the scope of computational reproduction.
4. **Not independent validation:** all agreement is with the same deposited data. Replicates are treated as samples
   (as in the paper), which inflates per-sample statistics.

## 4. Comparison with run 0006 (original tools on the same data)

| Item | 0006 (R/GSEA/original tools) | 0007 (self-implemented) |
|---|---|---|
| WP / PP / YP / integrative ARI at paper settings | 1 / 1 / 1 / 1 | 1 / 0.73 / 1 / 1 |
| paper settings' cophenetic rank (WP / PP / YP) | 8 / 3 / 7 | 10 / 4 / 11 |
| GSEA WP-Metabolism recovery (phenotype / gene-set permutation) | 4/11, 10/11 | 4/11, 10/11 |
| WGCNA modules WP / PP | 32 / 44 | 25 / 27 |
| kinome NMF ARI | 0.81 | 0.90 |
| Fig. 8B (37 testable kinases in top-40) | yes | yes |
| INKA, six-matrix kinase activity, PTM-SEA | not executed | executed (above) |

The two runs agree on everything that both executed, except the borderline AGS PP assignment. The core
conclusions do not depend on the specific software.

## 5. Deviations and assumptions

- *Inputs*: Table S3 (author-processed) replaces Spectronaut/MaxQuant output. Two non-human WP rows were removed
  (iRT, TRYP_PIG), leaving 5,675 proteins, and 7 duplicated PP site IDs were collapsed (8,445 sites).
- *NMF*: random-uniform initialisation, 100 runs, seed 123456. Consensus clusters come from average-linkage on
  1 − consensus. Silhouettes are computed on 1 − consensus.
- *Consensus clustering*: 500 resamplings, 80 % items, average linkage on 1 − Pearson, final average linkage on
  1 − consensus.
- *GSEA*: 1,000 permutations, Signal2Noise, gene-set size 15–500, MSigDB v7.0. Both permutation modes are
  reported because the paper does not state which it used.
- *WGCNA-type*: power = first β with scale-free signed R² ≥ 0.8 (WP 7, PP 4); unsigned TOM; own recursive tree
  cut (branch split when both children ≥ 30 features and the merge-height gap > 0.1); merging at eigengene
  correlation > 0.75. This is **not** the dynamicTreeCut algorithm, and many features stay unassigned.
- *GO enrichment*: hypergeometric ORA on MSigDB v7.0 GO-BP, replacing the PANTHER web tool
  (geneontology.org is unreachable).
- *Kinase activity*:
  - input is per-site centred log2 values;
  - protein-normalised = log2 PP − log2 WP of the gene;
  - regulons use exact symbol + site matching (PhosMap uses substring matching; tested in `06b`);
  - KSEA uses 1,000 permutations with one random ordering per sample, prefixes giving the size-n null sets;
  - MLR is ridge regression without intercept, standardised predictors, 10-fold CV over 100 λ values;
  - group test: Welch t, P < 0.05, |log2FC| > 1 (for mean: log2 ratio of group means of 2^activity; otherwise the
    difference of means).
- *INKA*: per sample on 2^log2 intensities. Kinase-centric = all pY sites on the kinase + activation-loop sites
  (InKA phomics table). Substrate-centric = PhosphoSitePlus substrates + the top-scoring NetworKIN kinase per site.
  Score = √(KC × SC). Only tyrosine kinases are ranked, as in Fig. 7A.
- *PTM-SEA*: 15-aa windows from UniProt (8,108 / 8,445 sites; 8,099 unique), centred values, 1,000 permutations,
  minimum overlap 10, PTMsigDB v1.9.1 flanking signatures. EMT vs Metabolism compared by moderated t on NES.

## 6. Network domains

Used: browse.quantms.org, www.ebi.ac.uk (Europe PMC supplements, QuickGO), rest.uniprot.org,
conda.anaconda.org (Python/R runtime), github.com (data resources inside the cloned repositories from run 0006).

No further domain access is **required** for this independent reproduction. Optional, for closer fidelity:

- www.gsea-msigdb.org / data.broadinstitute.org (MSigDB release matching GSEA 4.0.3);
- www.phosphosite.org (current kinase–substrate and regulatory-site data; Fig. S5A);
- geneontology.org / pantherdb.org (the paper's GO enrichment service);
- ftp.ebi.ac.uk (GOA files).
