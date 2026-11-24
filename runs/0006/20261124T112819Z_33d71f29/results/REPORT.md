# Reproduction of Hu et al. 2024 (Mol Oncol 18:1437) — PXD048538, quantms e60b9e9aa4c8

Run `0006` · prompt `prompt_0006.yaml` · scope: data processing, subtype classification, functional annotation,
kinase-activity analysis (the machine-learning model and public-database mining were excluded as instructed).

> **Status: partial reproduction.** Every step that could run with packaged, publicly available tools was executed.
> Three kinase/PTM steps (**InKA**, the **PhosMap mean/KSEA/MLR six-matrix analysis**, and **PTM-SEA**) were prepared
> but **not executed**. The original code is only on GitHub, and running cloned GitHub code was refused by this
> session's permission policy (section 7). Their outcomes are therefore **unknown**, not "not reproduced".

## 1. Summary

| Paper step | What was done here | Reproduction class | Outcome |
|---|---|---|---|
| Upstream search & quantification (Spectronaut 13 / MaxQuant 1.6) | Not possible (Spectronaut is commercial; raw PP/YP not reprocessed). Used (a) the authors' processed Table S3 and (b) the quantms DIA-NN reprocessing of WP | not reproduced / substituted | Table S3: 5,675 WP proteins (paper text: 6,091), 8,445 PP sites, 4,780 pY sites |
| Replicate Pearson correlation (Fig. S1C) | Recomputed on Table S3 | exact (same data, same statistic) | **0.931 / 0.826 / 0.721** vs paper 0.93 / 0.82 / 0.72 |
| GO cellular-component landscape (Fig. 1B) | QuickGO exact-term Swiss-Prot annotations (2026) | approximate (annotation release differs) | proteins: Spearman ρ = 0.99 vs paper counts; phosphoproteins ρ = 0.91 |
| NMF subtyping, MAD grid, cophenetic (Fig. 2A, 4A, 4C, S2A, S4A/E) | R `NMF` (brunet, 100 runs), MAD 10–100 %, k 2–5 | exact implementation, same input | With the paper's chosen parameters, **WP, PP and YP clusters are identical to the paper (ARI = 1.0)**. But the grid **does not** single out the paper's parameters as most stable (section 3.2) |
| Integrative clustering (Fig. 5A) | `ConsensusClusterPlus` on the binary subtype table, 500 reps | exact implementation | **Identical EMT / Metabolism split (ARI = 1.0)** |
| Subtyping on independently processed WP (quantms) | Same NMF settings | independent processing | Metabolism cluster identical; EMT-side split differs (ARI = 0.78) |
| GSEA (Fig. 2B, 4B, 4D, Table S6) | GSEA 4.3.2 CLI (paper 4.0.3), MSigDB v7.0 | approximate (minor version; permutation type not reported) | Same direction and ranking of terms; far fewer pass FDR < 0.25 with phenotype permutation. Gene-set permutation recovers more (WP Metabolism 10/11 Table S6 sets) |
| WGCNA modules (Fig. S2C, S4D) | `WGCNA` 1.73, default-type parameters | approximate (parameters not reported) | 32 WP / 44 PP modules vs paper 12 / 5. **Not reproducible from the description** |
| ANOVA / limma + GO (Fig. S2D, S4C/G) | one-way ANOVA / limma + hypergeometric ORA (MSigDB v7.0 GO-BP) | approximate (PANTHER web tool unreachable) | Metabolism: lipid/organic-acid metabolism; YP: RTK/ephrin/integrin signalling. Partly concordant |
| Integrated network proteins (Fig. 5B) | limma EMT vs Metabolism | exact statistic, same data | All 16 proteins the paper names at protein (or both) level are significant in the paper's direction; EMT phosphosites concordant; the named Metabolism phosphosites are absent from Table S3 |
| Kinome counts & families (Fig. 6A/B) | InKA/HGNC kinase list + UniProt families | approximate (kinase list and family source differ) | YP TK families agree well (FGFR 30 % vs 33 %, EPH 26 vs 33, SRC 12 vs 15, EGFR 6.4 vs 8.8). WP/PP group shares disagree |
| Outlier kinases (Fig. S6A–C) | top-5 % rule, two operationalisations | approximate (rule ambiguous) | ERBB2 & CDK12 (N87), FGFR2 (KatoIII), CSNK1D S384 / CSNK1E S389 (MKN28) **all recovered**. Total count depends on the definition (81 or 265 vs 124) |
| CDK4/6 MKN28 vs YCCEL1 (Fig. S6G) | t-test on WP | exact | CDK4 P = 1.1e-4 (paper 1e-4); CDK6 P = 0.066 (paper 0.021) |
| Kinome NMF (Fig. S8B/C) | NMF k = 2 on WP kinases | exact implementation | 40/42 samples match the integrative split (ARI 0.81); paper reports a perfect split |
| 40 differentially expressed kinases (Fig. 8B) | limma on WP kinases | exact statistic | All 37 of the paper's 40 kinases present in our kinase list rank in our top 40. Clustering on them gives ARI = 1 vs integrative subtypes |
| InKA on YP (Fig. 7A/D, S7A) | inputs built, **not run** | — | unknown |
| PP kinase activity, 6 matrices (Fig. 8A, S8A) | script written, **not run** | — | unknown |
| PTM-SEA (Fig. S5D) | flanking windows built, **not run** | — | unknown |
| PhosphoSitePlus functional sites (Fig. S5A) | not attempted (phosphosite.org unreachable; licence) | — | unknown |

**Bottom line.** Given the authors' own processed matrices and their chosen parameters, the paper's central
classification results reproduce exactly: the WP (3), PP (2), YP (3) and integrative (EMT vs Metabolism) subtypes.
The functional picture reproduces qualitatively: EMT/proliferation versus oxidative/lipid metabolism, and the named
EMT and metabolism marker proteins. Four things do **not** reproduce or cannot be verified:

1. the parameter *selection* (the cophenetic criterion does not pick the reported settings);
2. the strength of the GSEA statistics;
3. the number of WGCNA modules;
4. all kinase-activity inference results (not run).

This is **computational** reproduction on the same data; none of it constitutes experimental validation. The
western blots, inhibitor assays and siRNA experiments were not and cannot be reproduced here.

## 2. Data and processing

**Sources.**
- *Table S3* (author-processed; log2, quantile-normalised, imputed, batch-corrected): WP 5,677 rows, PP 8,452 rows
  and YP 4,780 rows. Two WP rows were removed (Biognosys iRT spike-in and porcine trypsin), giving 5,675.
  Seven duplicated PP site IDs were collapsed, giving 8,445.
- *quantms e60b9e9aa4c8*: DIA-NN 2.5.1 reprocessing of the 42 WP runs only; checksums verified against
  `provenance.json`. 8,255 protein groups pass 1 % FDR in at least one sample, and 6,092 are quantified in ≥ 50 %
  of samples. quantms has **no PP or YP layer**, so every phospho analysis necessarily uses Table S3.

**Paper processing steps vs this run.**
- *Quantile normalisation*: applied to quantms WP (missing-value-aware). Already applied in Table S3.
- *Imputation*: Spectronaut's imputation is proprietary. A Perseus-style down-shifted normal (shift 1.8 SD,
  width 0.3 SD, seeded) was used for quantms WP. This is an **approximation**.
- *proBatch*: not applied. Batch metadata are absent from the quantms SDRF, and Table S3 is already corrected.
  proBatch is not packaged for R 4.3, and installing it from GitHub was refused.
- *Imputed values*: they cannot be distinguished in Table S3, so per-line identification counts (Fig. S1A/B)
  cannot be reproduced from it.

**Counts vs paper.** The paper reports 6,091 proteins, 16,685 phosphopeptides and 8,705 pY sites, but its Table S3
contains fewer features (5,677 / 8,452 sites / 4,780). The text counts identifications, while Table S3 holds the
filtered analysis matrices.

**Replicate correlation** (`tables/02_replicate_correlation.csv`). Averaging over replicate pairs of the same line
gives 0.931 (WP), 0.826 (PP) and 0.721 (YP), matching the paper's 0.93 / 0.82 / 0.72 to two decimals. The paper's
"mean Pearson score" is therefore the within-line replicate mean. quantms WP: 0.960.

## 3. Subtype classification

### 3.1 Final subtypes (`tables/04_subtypes_vs_paper.csv`, `figures/04_consensus_*`)

The NMF settings follow the Methods: brunet algorithm, nrun = 100, seed 123456, features ranked by mean absolute
deviation. Using the paper's final settings:

| dataset | features | k | cophenetic | mean silhouette (paper) | ARI vs paper |
|---|---|---|---|---|---|
| WP (Table S3) | top 10 % = 568 (paper 568) | 3 | 0.956 | 0.74 (0.87) | **1.00** |
| PP (Table S3) | top 20 % = 1,689 (paper 1,691) | 2 | 0.954 | 0.87 (1.00) | **1.00** |
| YP (Table S3) | 100 % = 4,780 (paper 4,780) | 3 | 0.926 | 0.62 (0.92) | **1.00** |
| WP (quantms) | top 10 % = 610 | 3 | 0.921 | 0.55 | 0.78 |

- **WP.** The clusters are {AGS, eAGS, MKN28, N87} = Cell-prolif/EMT, {AZ521, SNU1750, TMK1} = EMT, and the other
  seven lines = Metabolism. All are identical to Fig. 2A.
- **YP.** The reproduction even matches the paper's oddity that the two Ocum1 replicates fall into different YP
  clusters (Fig. 4C).
- **Silhouette widths are lower than reported.** The paper does not say how silhouettes were computed; they may
  come from a different distance or the best single run rather than the consensus.
- **quantms WP.** The Metabolism cluster (7 lines) is identical. The two EMT-side clusters regroup as {AGS, eAGS}
  versus {AZ521, MKN28, N87, SNU1750, TMK1}. So the 3-cluster EMT/prolif-EMT split depends on upstream
  processing, while the EMT-vs-Metabolism dichotomy is robust to it.

### 3.2 Parameter selection is not reproduced (`tables/04_nmf_grid_all.csv`, `figures/04_FigS2A_S4AE_cophenetic_grid`)

The paper states that the final MAD fraction and k were "the best combination … that resulted in the most stable
k-cluster decomposition". In our grid of 20 combinations per dataset, the paper's settings rank:

- **WP**: 8th. The maximum is 10 %, k = 5 (cophenetic 0.970 vs 0.956).
- **PP**: 3rd. The maximum is 10 %, k = 2.
- **YP**: 7th. The maximum is 10 %, k = 2.

Cophenetic coefficients are high (> 0.85) almost everywhere, so the criterion barely discriminates. Choosing the
reported settings therefore required extra judgement that the Methods do not describe (e.g. the subjective reading
of Fig. S2A/S4A/E, silhouettes, or interpretability). This is a reproducibility gap in the *selection step*, not in
the clustering itself.

The Methods say "mean absolute deviation", but the Fig. S2 legend says "median absolute deviation". A median-AD
grid for WP was run as a sensitivity check (`tables/04_nmf_grid_median_ad.csv`). At 10 % / k = 3 it gives the
**same three WP clusters** as the paper, so the subtyping is insensitive to this ambiguity. Its cophenetic optimum
is again elsewhere (10 %, k = 2; 0.973 vs 0.963), so it does not explain the parameter choice either.

### 3.3 Integrative clustering (`tables/03_integrative_clusters.csv`, `figures/04_Fig5A_*`)

`ConsensusClusterPlus` was run on the 8 × 28 binary subtype table (2 + 3 + 3 subtypes; replicates r1/r2), with
500 reps, pItem 0.8, hierarchical clustering and Pearson distance.

- At k = 2 it returns **S1/EMT = AGS, eAGS, AZ521, N87, MKN28, SNU1750, TMK1** and **S2/Metabolism = KatoIII,
  NCC24, Ocum1, SNU601, SNU719, YCC3, YCCEL1**. This is exactly Fig. 5A, with both replicates concordant for every
  line and a consensus of 1.0.
- The integrative split equals the PP split. That is unsurprising, because WP clusters 1+2 and YP clusters 1+2 are
  nested within it apart from YCC3 and N87.

### 3.4 NMF metagenes (Fig. S3A)

`extractFeatures` (Kim–Park) gives 54 WP metagene-specific proteins, of which 23 are among the paper's 42. 40 of
the paper's 42 are in our top-10 % feature set. The paper's "42" is presumably a filtered or ranked subset; its
selection rule is not described.

## 4. Functional annotation

**GSEA** (`tables/05_gsea_all.csv`, `05b_gsea_genesetperm.csv`, `05c_gsea_vs_paper_summary.csv`,
`figures/05_Fig2B_4BD_GSEA_*`).
- *Setup*: GSEA 4.3.2 CLI, Signal2Noise, 1,000 permutations, each cluster vs the rest, MSigDB v7.0. Phosphosites
  were collapsed to genes with GSEA's own Max_probe mode.
- *Direction and ranking match the paper*:
  - Metabolism: OXPHOS, bile-acid, fatty-acid, xenobiotic, adipogenesis, peroxisome, KRAS-up.
  - Cell-prolif/EMT: G2M, E2F, MYC, mitotic spindle, TGF-β, TNFα, PI3K-AKT-mTOR.
  - EMT cluster: hypoxia, EMT, myogenesis.
  - PP metabolism cluster: NK-cytotoxicity, calcium, GnRH, gap-junction, LTP. YP: JAK-STAT, ERBB, NK.
- *Significance is much weaker*:
  - With phenotype permutation only 4 / 0 / 0 Hallmark sets pass FDR < 0.25 (Metabolism / prolif-EMT / EMT),
    against 11 / 27 / 18 in Table S6.
  - With **gene-set permutation** (not stated in the paper; Table S6's many FDR = 0 values suggest it) the
    recovery is: Metabolism 10/11, prolif-EMT 7/27, EMT 3/18, PP 5/11, YP 2/10.
  - The paper's EMT-cluster and prolif/EMT lists contain many sets (interferon, complement, allograft rejection)
    that rank only moderately here.
  - The 42 samples are 14 lines × 3 replicates. Replicate pseudo-replication inflates any per-sample statistic,
    in the paper and here.
- **Conclusion**: qualitatively reproduced, quantitatively not; the result depends on unreported GSEA settings.

**WGCNA** (`tables/05_WGCNA_*`). The paper gives no power, network type, module size or merge height. With a
scale-free-fit power (WP 7, PP 4), unsigned TOM, minModuleSize 30 and mergeCutHeight 0.25, we obtain **32 WP and
44 PP modules vs 12 and 5**.
- Module–cluster associations exist: a "brown" WP module correlates with the Metabolism cluster (r = 0.92) and a
  "red" module with the EMT cluster (r = −0.96).
- The module structure itself cannot be reproduced without the original parameters.

**ANOVA / limma + GO.**
- *Counts*: WP ANOVA gives 2,239 proteins at BH < 0.05; YP ANOVA 754 sites; PP limma 758 sites
  (|log2FC| > 1, adj.P < 0.05).
- *Enrichment*: the PANTHER web service and Bioconductor GO data were unreachable. Hypergeometric ORA on MSigDB
  v7.0 GO-BP was used instead.
- *Concordant themes*: lipid/organic-acid catabolism (Metabolism cluster); integrin/ephrin/RTK signalling and
  Fc-receptor signalling (YP clusters); actin organisation (PP).
- The specific novel terms highlighted in the paper (snRNA 3′-end processing, negative regulation of DNA-damage
  response) were not among our top terms.

**Integrated subtype markers (Fig. 5B)** (`tables/05_Fig5B_named_proteins_check.csv`).
- All 16 proteins that the paper highlights at protein level (or at both levels) are significantly changed in the
  paper's direction (|log2FC| > 1, adj.P < 0.05): AAMP, CCN1, LGALS1, TGFBI, STK3, STK4 and PHLDB2 up in EMT;
  GLUL, ALAD, ACAA1, ACAA2, GCLC, AK3, ALDH3A1, ACSL5 and AKR1A1 up in Metabolism. Proteins named only at phospho
  level (TGFB1I1, CTNNB1, YAP1, G6PD) show small protein changes in the same direction.
- EMT phosphosites on TGFB1I1, CTNNB1, YAP1, ABL2, APC and PHLDB2 include significant up-sites.
- For the metabolism phospho-proteins (G6PD, GSTA2, ALDH3A1, ACSL5, AKR1A1), **no phosphosite is present in
  Table S3's PP matrix**. That part of Fig. 5B cannot be checked from the deposited data.

## 5. Kinase analyses that were executed

- **Kinome size** (`tables/08_kinome_counts.csv`): 178 / 158 / 155 kinases in WP / PP / YP, against 198 / 193 / 140
  in the paper. The kinase list here is InKA's 2016 HGNC/Manning list; the paper used UniProt + PhosphoSitePlus.
  Three of the paper's Fig. 8B kinases (COQ8A, TGM2, MAP3K20) are absent from that list.
- **Family shares (Fig. 6A/B).**
  - YP tyrosine-kinase families reproduce well: FGFR 30.5 % (paper 33.0), EPH 25.6 (32.9), SRC 11.9 (14.7),
    EGFR 6.4 (8.8), JAK 3.5 (3.8), MET 2.8 (3.4).
  - WP/PP group shares do not. The paper reports e.g. 22 % "Atypical" in WP and 20 % CAMK in PP; we get 4.9 % and
    1.5 %. Its family assignment (and possibly the inner/outer ring of Fig. 6A) evidently differs from UniProt's,
    and the source is not given.
- **Outlier kinases (Fig. S6).** "Top 5 % of features with the highest level across all samples" was
  operationalised in two ways (raw level or per-feature z-score, top 5 % of kinase × sample values). That gives 81
  or 265 outlier kinases vs 124 in the paper; the count is definition-dependent. The paper's named examples are
  recovered:
  - ERBB2 protein and phosphosites in N87, and CDK12 protein in N87;
  - FGFR2 PP and YP sites in KatoIII, uniquely;
  - CSNK1D S384 and CSNK1E S389 in MKN28 (z-score definition).
- **Fig. S6G.** CDK4 is higher in YCCEL1 than MKN28 (P = 1.1 × 10⁻⁴, paper 1 × 10⁻⁴). CDK6: P = 0.066, against
  paper 0.021.
- **Kinome NMF (Fig. S8B/C).** 179 WP kinases, k = 2: cophenetic 0.91, silhouette 0.73. 40/42 samples follow the
  integrative split (AGS_r1 and MKN28_r3 are misassigned; ARI 0.81), whereas the paper reports silhouette 1.0. The
  EMT/Metabolism kinome segregation is reproduced in substance.
- **Fig. 8B, 40 kinases.**
  - All 37 of the paper's kinases present in our kinase list rank in our top-40 by limma (EMT vs Metabolism).
  - Hierarchical clustering on them separates the integrative subtypes perfectly (ARI 1.0).
  - Of the 37, 17 are higher in EMT and 20 higher in Metabolism.

## 6. Interpretation

1. The deposited processed data (Table S3) carry the paper's subtype structure unambiguously. Given the same matrix
   and parameters, the original R implementations return the published WP, PP, YP and integrative assignments
   exactly.
2. The **analytic degrees of freedom** sit upstream and around the clustering, not in it:
   - parameter selection (MAD/k);
   - GSEA permutation scheme;
   - WGCNA parameters;
   - the outlier and kinase-family definitions;
   - upstream quantification. Independent DIA-NN reprocessing preserves the EMT-vs-Metabolism dichotomy but not the
     finer EMT split.
3. The kinase-*activity* conclusions remain unverified in this run: SFK/SRC as a pan-GC activity (InKA) and
   mTOR/CDK7/LATS1 (EMT) vs MAP2K2/MKNK2/STK3/EIF2AK1 (Metabolism) from the six-matrix analysis. They are the
   paper's main translational claims. The supporting drug-sensitivity experiments are experimental and outside
   computational reproduction.
4. **Evidence is insufficient to call the full workflow reproduced.** It is sufficient to call the subtype
   classification reproduced given the authors' processed data and parameters.

## 7. Steps not executed and why

The original implementations are only distributed through GitHub. They were cloned at pinned commits (see
`run_info/github_tool_commits.txt`), and inputs and driver scripts were prepared. Executing the cloned code was
refused by the session's automatic permission policy, once for `R CMD INSTALL` and twice for running InKA, so:

- **InKA 1.6.0** (`digitalbiologylab/inka-intensity`): YP sites were converted to MaxQuant-style
  `Phospho (STY)Sites.txt` / `pp-Peptide-report.txt`, with one fully tryptic pseudo-peptide per pY site and
  intensities = 2^log2 (`scripts/07a_inka_inputs.py`). Runner: `scripts/07b_inka_run.R`. The paper used the InKA
  web server (inkaweb.nl), which is unreachable from here.
- **PhosMap** (`ecnuzdd/PhosMap`) mean / KSEA / MLR on raw and protein-normalised PP: `scripts/06_pp_kinase_activity.R`
  sources the unmodified PhosMap functions and applies the paper's assignment rule (P < 0.05, |log2FC| > 1,
  counting over 6 matrices).
- **PTM-SEA** (`broadinstitute/ssGSEA2.0`, PTMsigDB v1.9.1): 15-aa flanking windows were built for 8,108/8,445 PP
  sites (`scripts/04a_ptm_windows.py`). The paper used the GenePattern web module, which is unreachable.

If running these repositories is permitted, the three scripts can be executed as they are.

## 8. Deviations and assumptions (complete list)

- *Upstream*: Spectronaut 13 → not available. Table S3 was used for the downstream reproduction and quantms
  DIA-NN for the independent WP check.
- *GSEA*: 4.0.3 → 4.3.2 (bioconda). Gene sets: MSigDB v7.0. Permutation type is unspecified, so both modes were run.
- *GO enrichment*: geneontology.org (PANTHER) → hypergeometric ORA on MSigDB v7.0 GO-BP (fgsea::fora). GO CC
  counts come from QuickGO (2026 annotations).
- *WGCNA*: parameters are not in the paper; defaults as listed in section 4. WGCNA needs the Bioconductor data
  packages GO.db and GenomeInfoDbData (bioconductor.org is blocked). Empty placeholder packages were installed in a
  separate library (`/opt/stubs`) only to satisfy the namespace import; no GO.db functionality is used.
- *Kinase universe*: the InKA HGNC/Manning list (2016), not UniProt + PhosphoSitePlus. Families come from UniProt
  (2026).
- *Replicates*: treated as samples throughout, as in the paper's figures. This pseudo-replication inflates
  significance.
- *Duplicate IDs*: duplicate gene symbols in WP were made unique (GENE|ACCESSION); duplicate PP site IDs were
  collapsed.

## 9. Network domains

Reached: browse.quantms.org, www.ebi.ac.uk (Europe PMC, QuickGO), rest.uniprot.org, conda.anaconda.org, github.com
(git clone only).

Unreachable but needed for a more faithful run:
- bioconductor.org (GO.db, org.Hs.eg.db, GenomeInfoDbData);
- cloud.r-project.org (CRAN);
- www.gsea-msigdb.org / data.broadinstitute.org (GSEA 4.0.3, MSigDB);
- inkaweb.nl (InKA server);
- cloud.genepattern.org (PTM-SEA);
- www.phosphosite.org (PSP regulatory sites, Fig. S5A);
- geneontology.org / pantherdb.org (GO enrichment);
- ftp.ebi.ac.uk / ftp.pride.ebi.ac.uk (GOA, raw files);
- depmap.org, cancerrxgene.org, cbioportal.org (out of scope here).
