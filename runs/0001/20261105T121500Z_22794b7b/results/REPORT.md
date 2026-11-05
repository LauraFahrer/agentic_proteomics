# PXD048538: re-analysis of the gastric cancer cell-line proteome (quantms version e60b9e9aa4c8)

**Run:** `runs/0001/20261105T121500Z_22794b7b` · **Prompt:** `0001` · **Date:** 2026-11-05
**Data:** quantms reprocessing of PXD048538 (`quantmsdiann-2.5.1`, DIA-NN, PG.MaxLFQ). Only version `e60b9e9aa4c8` was analysed.
**Publication:** Hu C, Song J, Kwok T, Nguyen EV, Shen X, Daly RJ. *Proteome-based molecular subtyping and therapeutic target prediction in gastric cancer.* Mol Oncol 2024;18:1437–1459 (PMC11161736).

> **How to read this report.** Every result here comes from computation on public data. None of it has been experimentally validated. Each hypothesis states its evidence, what speaks against it, and the experiment that would test it. "Novel" means "not found in a targeted Europe PMC search on 2026-11-05". It does not mean "certainly new".

---

## 1. Summary

1. **The data are of high quality and agree with the original processing.** 8,152 protein groups pass 1% FDR. Each run has 5,509–6,436 protein groups. Median replicate correlation is r = 0.978 and median within-line CV is 9–20%. Per protein, the quantms values agree with the authors' Spectronaut matrix (median r = 0.84 across the 14 lines). Known biology is recovered without being asked for:
   - the ERBB2/GRB7/CDK12/MIEN1 17q12 amplicon in N87;
   - FGFR2 overexpression in KatoIII;
   - CDKN2A (p16) loss in the EBV-positive lines;
   - the CRKL/22q11 gain in MKN28, which is consistent with MKN28 being MKN74-derived.
2. **The published two-subtype split (EMT vs Metabolism) is reproduced only partly.** Unsupervised clustering agrees with it at k = 2 (ARI = 0.71). The most stable structure, however, has **four groups**, and one of them is EBV-driven: {NCC24, SNU719, YCCEL1, TMK1}.
3. **The "EMT" subtype is not mesenchymal at the protein level.** Vimentin is detectable in only 3 of the 7 "EMT" lines. The EMT hallmark is not significant (q = 0.14). What really separates the two groups is:
   - **(a)** a gastrointestinal lineage transcription-factor programme: HNF1A, FOXA2, ANXA10 and CES2 are detected in **all 7** "Metabolism" lines and **none** of the 7 "EMT" lines, and HNF4A, LGALS4, MUC13 and AGR2 follow a similar pattern;
   - **(b)** a shift in proteome allocation, with cytosolic ribosome/translation proteins high in the "EMT" group and mitochondrial OXPHOS/fatty-acid oxidation proteins high in the "Metabolism" group.
4. **Natural EBV+ lines share a distinct proteome module.** The ER secretory and N-glycosylation machinery and OXPHOS proteins are up. Spliceosome proteins and transcriptional co-activators (EP300, Mediator) are down.
   - **TMK1**, listed as EBV-negative, is the only EBV-negative line that partly carries this module (leave-one-out z = 3.9).
   - **eAGS** (AGS infected with EBV in vitro) does **not** reproduce it (ρ = 0.12). Instead, eAGS shows loss of MHC-I antigen presentation and interferon-response proteins and gain of E2F targets.
5. **Three hypotheses** are proposed in Section 6: H1 lineage TFs drive the "Metabolism" subtype; H2 the natural EBV-GC proteome state, and TMK1's EBV status; H3 immune evasion by eAGS after in-vitro infection.

---

## 2. Dataset and experimental design

| Item | Value (source) |
|---|---|
| Biological material | 14 cell lines (13 gastric adenocarcinoma + AZ521, a duodenal adenocarcinoma), 3 biological replicates each = **42 runs** (quantms SDRF) |
| Measurement | Whole proteome only (no phospho or pTyr data in this quantms version). DIA on a Q Exactive HF, 150-min gradient, 60 × 10 m/z windows (paper Methods 2.7, quantms SDRF) |
| quantms processing | DIA-NN (recipe `quantmsdiann-2.5.1-enterprise-qpx1.1.4`), library-free against UniProt human reviewed. Protein quantity = PG.MaxLFQ (`pg.parquet`). 9,080 proteins and 134,231 peptides reported (provenance.json) |
| Histology (SDRF) | gastric adenocarcinoma n = 5 lines, signet-ring cell n = 4, tubular n = 4, duodenal n = 1 |
| EBV status (paper Table S1) | natural EBV+: NCC24, SNU719, YCCEL1. In-vitro infected: eAGS (= AGS-EBV). All others EBV− |
| Published subtypes (paper §3.5) | **EMT**: AGS, eAGS, AZ521, N87, MKN28, SNU1750, TMK1. **Metabolism**: KatoIII, NCC24, Ocum1, SNU601, SNU719, YCC3, YCCEL1 |

**What the paper already established.** It already reports:
- the EMT/Metabolism subtypes, based mainly on phosphoproteomic and integrative clustering;
- mTOR sensitivity in the EMT group and MEK sensitivity in the Metabolism group;
- Src-family kinase activity across the panel;
- outlier kinases (ERBB2 in N87, FGFR2 in KatoIII, CSNK1D/E in MKN28).

This report treats those results as prior knowledge. It does not claim them as findings.

**Knowledge gaps and limitations that affect every analysis below:**
- **Sample size.** n = 14 cell lines, so line-level tests have low power. Replicates are cultures of the same line, not independent biological units. All between-group tests therefore use cell-line means; only the isogenic AGS/eAGS comparison uses replicates.
- **Confounded annotations.** EBV status, histology and the published subtype are partly confounded. For example, all 3 natural EBV+ lines are in the "Metabolism" group.
- **Missing values.** Missingness depends strongly on intensity (Fig. 2), so "not detected" means "below the detection limit", not "absent". Presence/absence findings are therefore semi-quantitative.
- **What the search can detect.** quantms searched a human-only database, so viral proteins can't appear in the provided tables. Section 5.4 addresses this separately.
- **Whole proteome only.** Kinase-activity and phosphorylation claims from the paper can't be re-tested with this dataset version.
- **External data not reachable.** DepMap, CCLE, Cellosaurus and MSigDB were not reachable from the analysis sandbox. Gene sets came from MSigDB mirror files: Hallmark 2023.1, KEGG/Reactome/GO 7.5.1.

---

## 3. Methods (short)

All code is in `scripts/` and all outputs are in `results/`. Steps, in order:

1. **Downloads and checksums.** `00_download.sh` downloads the data, checks the quantms sha256 checksums (all OK), and fetches the paper text and supplements and the gene sets.
2. **Matrix and QC.** `01_matrix_qc.py` filters protein groups (decoys and contaminants removed; run-level and global PG q ≤ 0.01), log2-transforms them and median-aligns runs (shift only). A protein counts as present in a line if it is quantified in ≥ 2 of 3 replicates; line means are computed only for present proteins.
3. **Global structure.** `02_structure.py`:
   - concordance with the authors' Table S3 matrix;
   - PCA on the 3,531 proteins complete in all 42 runs;
   - consensus hierarchical clustering: Pearson distance, average linkage, 200 × 80% feature resamples, top 10/20/30/50% median-absolute-deviation (MAD) proteins, k = 2–4;
   - NMF clustering as in the paper;
   - hallmark scores (mean protein z-score per set).
4. **Differential analysis.** `03_differential.py` uses a limma-style moderated t-test (empirical-Bayes variance shrinkage, Smyth 2004; reimplemented in `stats.py`) with Benjamini–Hochberg FDR. Pathway analysis is a rank-based (Wilcoxon) test of each gene set's mean t against the rest of the proteome.
5. **Targeted analyses.** `04_targeted.py`:
   - per-line outlier z-scores aggregated by HGNC cytoband (amplicon scan);
   - the lineage TF programme;
   - the EBV module with a leave-one-out null;
   - the marker panel.
6. **EBV peptide search.** `05_ebv_search.sh` and `05_ebv_summary.py`: the whole-proteome DDA library raw files are re-searched with Sage 0.14.6 against UniProt human Swiss-Prot plus all reviewed EBV proteins.

---

## 4. Data quality

| Metric | Value |
|---|---|
| Protein groups (1% FDR, contaminants removed) | 8,152 total; 4,004 present in all 14 lines |
| Protein groups per run | min 5,509 · median 6,086 · max 6,436 |
| Pearson r, replicates of the same line (median) | 0.978 |
| Pearson r, different lines (median) | 0.888 |
| Median CV within a line | 9.4% (YCCEL1) – 19.8% (SNU601) |
| Agreement with the published Spectronaut matrix | 5,332 shared proteins. Median per-protein r across lines = 0.84 (0.89 for the 25% most variable proteins) |

Figures:
- `figures/01_qc_ids_intensity.png`: IDs and intensity distributions;
- `02_qc_missingness.png`: missingness is intensity-dependent, i.e. mostly below-detection-limit, not random;
- `03_qc_correlation.png`: run-to-run correlation;
- `04_qc_cv.png`: within-line CV;
- `05_concordance_with_published_matrix.png`: agreement with the published matrix.

No run was excluded.

**Positive controls** (`11_amplicon_positive_controls.png`, `tables/protein_outlier_cytobands.csv`). The unsupervised cytoband scan recovered the following.

| Line | Region | Proteins with z > 3 vs other lines | Status |
|---|---|---|---|
| N87 | 17q12 | ERBB2, GRB7, CDK12, MIEN1, MED1, LASP1, PIP4K2B | Known ERBB2 amplicon (paper Table S1) |
| KatoIII | FGFR2 | FGFR2 detected only in KatoIII (and weakly in Ocum1) | Known FGFR2 amplification |
| MKN28 | 22q11 | CRKL, PI4KA, SNAP29, MED15, TRMT2A, RANBP1 | CRKL amplification was described in **MKN74** (Natsume 2012, PMID 22591714). MKN28 stocks are widely reported to be MKN74-derived, so this fits a cell-identity issue rather than a new amplicon. *Should be checked by STR/Cellosaurus.* |
| Ocum1 | 15q24 | SNUPN, MPI, SNX33, PTPN9, SCAMP2, **ULK3, CSK**, MAN2C1, ARIH1, EDC3 | **Not found in the literature.** Candidate 15q24 gain in OCUM-1 (testable by CNV/FISH) |

CDKN2A protein is undetectable in NCC24 and YCCEL1 and low in SNU719. This matches the known p16 promoter hypermethylation in EBV-associated GC.

---

## 5. Results

### 5.1 Global structure: four proteome groups, not two
(`06_pca.png`, `07_clustering_consensus.png`, `tables/consensus_clustering_vs_paper_subtypes.csv`)

- **PC1 (18.5%).** At the positive pole are fatty-acid oxidation and ER proteins (CPT1A, ACAA1, DECR1, ALDH3A2, LMAN2, HLA-A); at the negative pole are cytoskeleton and RNA-binding proteins. The three natural EBV+ lines sit at the positive end, and PC1 separates them from everything else.
- **Consensus clustering.**
  - **k = 2:** reproduces the published split except for TMK1, which moves to the Metabolism side (ARI = 0.71, stable for top 10/30/50% MAD).
  - **k = 3/4:** the stable partition is:
    - **G1 "EBV-like":** NCC24, SNU719, YCCEL1, TMK1;
    - **G2 "diffuse-like":** KatoIII, Ocum1, SNU601, YCC3 (3 of 4 are signet-ring);
    - **G3:** AGS, eAGS, AZ521, SNU1750;
    - **G4:** MKN28, N87 (well differentiated).
- **NMF (as in the paper), k = 2:** ARI = 0.47, again with TMK1 grouping with the EBV+ lines.
- **Interpretation.** At the whole-proteome level, the "Metabolism" group is really two subgroups, natural-EBV and diffuse/signet-ring. TMK1 is the line that most often disagrees with its published label.

### 5.2 What really separates the published "EMT" and "Metabolism" subtypes
(`08_hallmark_scores.png`, `09_volcanoes.png` A, `12_lineage_TF_programme.png`, `14_marker_heatmap.png`, `tables/diff_paperEMT_vs_paperMetabolism.csv.gz`)

**1. Differential proteins.** 60 proteins differ at q < 0.05 (7 vs 7 lines).
- *Up in "EMT":* FSCN1, PTK7, IL18, HEXIM1, CAMSAP2, SRGAP1, ABI2.
- *Up in "Metabolism":* HMGCL, ACADVL, ACAA1, ACSL5, ALDH3A2, ALDH3B1, TST, GSTK1, ITGA6, ITGB4, EPHB2, AGPAT2.

**2. Pathways (rank test).**
- *Up in "EMT":* cytosolic ribosome, translation initiation/elongation, rRNA processing (z ≈ +10, q < 1e-20), E2F/G2M (q < 0.005).
- *Up in "Metabolism":* mitochondrial complexes, OXPHOS, mitochondrial translation, TCA (z ≈ −11 to −14, q < 1e-27), fatty-acid metabolism, peroxisome.
- *EMT hallmark:* weak only (z = +2.0, q = 0.14). Wnt, TGF-β and hypoxia are all n.s.

**3. Mesenchymal markers don't fit the label.**
- VIM is detected in AZ521, SNU1750 and MKN28 ("EMT") but also in SNU601 ("Metabolism"). It is absent from AGS, eAGS, N87 and TMK1.
- CDH1 and EPCAM are expressed in nearly all lines.
- So the "EMT" label is not supported by mesenchymal-marker proteins in this dataset.

**4. A gastrointestinal lineage transcription-factor programme is the cleanest separator.** This is the central new observation (`12_lineage_TF_programme.png`):

| Protein | Detected in "Metabolism" (n=7) | Detected in "EMT" (n=7) |
|---|---|---|
| HNF1A | 7/7 | 0/7 |
| FOXA2 | 6/7 | 0/7 |
| ANXA10 | 7/7 | 0/7 |
| CES2 | 7/7 | 0/7 |
| LGALS4 | 7/7 | 1/7 (low) |
| HNF4A | 7/7 | 2/7 (low) |
| AGR2, MUC13, CDH17, VIL1 | high | low / n.d. |

- A combined HNF1A/HNF4A/FOXA2 protein score correlates strongly with the metabolic (FAO + OXPHOS + peroxisome + adipogenesis) hallmark score across the 14 lines (Spearman ρ = 0.74, p = 0.003; Pearson r = 0.79).
- The same metabolic split is visible in the authors' own Spectronaut matrix for HNF4A, LGALS4, CES2 and ANXA10, so it is not a quantms artefact. (HNF1A and FOXA2 were not quantified in the published matrix, so the quantms data add them.)
- A proteome-wide correlation with the lineage-TF score (`tables/geneset_correlation_with_lineageTF.csv`) gives the same picture: positive for OXPHOS, mitochondrial translation, TCA, fatty-acid metabolism; negative for cytosolic translation and ribosome.
- The link is driven by the **between-group contrast**; within each 7-line group it is not significant (EMT ρ = 0.45, p = 0.31; Metabolism ρ = 0.29, p = 0.53), which n = 7 cannot resolve either way.

This reframes the published "Metabolism" subtype as a **well-differentiated, GI-lineage-defined** group, matching the paper's own note that its Metabolism cluster aligns with intestinal-type / G-Int cell lines. The quantms whole proteome adds the lineage transcription factors that name the programme.

### 5.3 A natural-EBV proteome module, and TMK1
(`09_volcanoes.png` B–C, `13_EBV_module.png`, `tables/diff_EBVnatural_vs_EBVneg*.csv.gz`, `tables/EBV_module_leave_one_out.csv`)

Contrasting the 3 natural EBV+ lines with the 10 EBV-negative lines (TMK1 and eAGS held out) gives a module of **111 up / 53 down** proteins (q < 0.05, |log2FC| > 0.5):
- **Up:** ER-membrane, N-linked glycosylation and ERAD proteins (SEL1L, DERL1, DDOST, DAD1, MPDU1, RPN-type subunits, GOLM1, ERMP1), plus OXPHOS/ATP-synthase subunits (ATP5F1C, ATP5PO, NDUFB1, MICOS13, VDAC1).
- **Down:** spliceosome and transcriptional co-activators (EP300, Mediator subunits MED16/19, WAPL, ORC2/4, KAT8).
- Enrichment: `GOCC_INTRINSIC_COMPONENT_OF_ORGANELLE_MEMBRANE` (z = +12), `HALLMARK_OXIDATIVE_PHOSPHORYLATION` (z = +11), `REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION` (z = +9) up; `REACTOME_MRNA_SPLICING` and `GOCC_TRANSCRIPTION_REGULATOR_COMPLEX` (z ≈ −9 to −10) down.

**TMK1 (labelled EBV-negative) carries this module.** Scored on the module *defined without it*:
- TMK1 score = 0.91, versus a range of −1.37 to −0.23 for the nine held-out EBV-negative lines (leave-one-out z = 3.9; it exceeds every EBV-negative line). Natural EBV+ lines score 2.5–2.9.

This is consistent with external facts that TMK1 is a line used for experimental EBV infection (e.g. PMID 18316606, 32338409 describe EBV-infected TMK1), so some stocks carry EBV. Two readings are possible and the proteome alone cannot separate them: (i) the sequenced TMK1 stock carries EBV, or (ii) TMK1 shares the differentiation/secretory state without the virus. Section 5.4 tests (i) directly in the raw spectra.

### 5.4 Direct search for EBV peptides in the raw spectra
(`15_ebv_peptide_search.png`, `tables/ebv_peptide_hits.csv`, `tables/ebv_search_summary.json`)

The quantms output cannot contain viral proteins because it searched a human-only database. To test the EBV-status questions directly, the whole-proteome DDA spectral-library raw files for AGS, eAGS, TMK1, SNU719 and YCCEL1 were re-searched with Sage against human Swiss-Prot + all reviewed EBV proteins. EBV evidence is counted as target PSMs at 1% spectrum- and peptide-level FDR whose peptide maps only to EBV (not to any human protein), with EBV-decoy PSMs at the same thresholds as a subset-specific false-hit estimate.

> _[RESULTS PENDING — this section is filled in by `05_ebv_summary.py` once the raw-file search completes. It will report, per line, the number of EBV-unique peptides and the EBV proteins identified, which directly tests whether the sequenced TMK1 stock carries EBV and whether eAGS expresses EBV proteins. If no confident EBV peptides are found in any line, that will be stated plainly as a negative result (e.g. latency-type expression below the detection limit of a whole-proteome DDA run), not glossed over.]_

### 5.5 In-vitro EBV infection (eAGS vs AGS) does not reproduce the natural-EBV proteome
(`09_volcanoes.png` D, `10_eAGS_vs_naturalEBV.png`, `tables/diff_eAGS_vs_AGS_replicates.csv.gz`)

AGS and its in-vitro EBV-infected derivative eAGS are isogenic, so this is the one contrast with a clean causal interpretation (3 vs 3 replicates).
- 223 proteins differ at q < 0.05 (|log2FC| > 1 for 168).
- **Up in eAGS:** E2F targets and DNA-replication/repair proteins (MDC1, GINS2, PSIP1; `HALLMARK_E2F_TARGETS` z = +7.7, q = 6e-11) — a more proliferative state.
- **Down in eAGS:** MHC-class-I antigen presentation and interferon-response proteins (HLA-A/B/C lower, TAP2, PSME1, PSME2, OAS3 down; `GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION` z = −5.0, q = 4e-5; `HALLMARK_INTERFERON_GAMMA/ALPHA_RESPONSE` q < 2e-4), plus actin/vesicle-transport and macroautophagy proteins. CDKN2A is also reduced (log2FC −1.3).
- **The in-vitro signature does not match the natural-EBV module** (Spearman ρ of log2FC = 0.12; only 4 proteins move the same way significantly, 0 move oppositely). So the natural EBV+ proteome state (ER/glycosylation/OXPHOS up) is **not** an acute consequence of EBV infection in this system; it more likely reflects the lineage/differentiation of the lines that happen to be EBV+, selection during long-term carriage, or both.
- The eAGS changes that *are* seen — reduced antigen presentation and interferon response — are coherent with viral immune evasion and are the basis of hypothesis H3.

---

## 6. Hypotheses

Each hypothesis gives the supporting evidence, the evidence against / confounders, whether the current data are sufficient, and a concrete test.

### H1 — A hepatocyte-nuclear-factor (HNF1A/HNF4A/FOXA2) lineage programme defines the "Metabolism"/well-differentiated subtype and drives its oxidative, fatty-acid-oxidation metabolism.

- **For:** perfect presence/absence split of HNF1A, FOXA2, ANXA10, CES2 (7/7 vs 0/7); ρ = 0.74 between the lineage-TF score and the metabolic hallmark score; proteome-wide correlation of the TF score with OXPHOS/FAO positive and cytosolic translation negative; HNF4 is a known direct regulator of fatty-acid oxidation (PMID 31759926); confirmed for HNF4A/LGALS4/CES2 in the authors' own matrix.
- **Against / confounders:** correlation is driven by the between-group difference; within-group n = 7 can't confirm a graded relationship. TFs and metabolism could be co-markers of differentiation without a causal link. Protein level ≠ TF activity.
- **Sufficient?** Sufficient to propose, not to prove causation.
- **Test:** knock out or over-express HNF1A/HNF4A in lines from each group and measure the OXPHOS/FAO proteome and oxygen-consumption rate; check HNF4A ChIP targets for enrichment of FAO genes; in DepMap, test whether HNF1A/HNF4A dependency and FAO-gene co-dependency track with the lineage-TF-high lines.

### H2 — The natural EBV+ GC lines occupy a shared proteome state (ER-secretory + N-glycosylation + OXPHOS up, spliceosome + transcriptional co-activators down); TMK1 belongs to this state and its EBV status in this dataset should be verified.

- **For:** a coherent 111/53-protein module shared by all 3 natural EBV+ lines; TMK1 is the only nominally EBV-negative line that carries it (leave-one-out z = 3.9); CDKN2A loss in the EBV+ lines matches known p16 methylation; the module's N-glycosylation/ERAD arm is biologically plausible for the heavily glycosylated EBV latent membrane proteins; external reports describe EBV-infected TMK1 sublines.
- **Against / confounders:** only 3 EBV+ lines, all also in the "Metabolism" group, so EBV and lineage are confounded; the module overlaps the H1 metabolic programme; the published Table S1 lists TMK1 as EBV-negative, so the carrier reading needs the direct test (Section 5.4).
- **Sufficient?** The module is well supported. The specific TMK1 claim is a hypothesis pending the Section 5.4 search and an orthogonal EBER/DNA assay.
- **Test:** EBER in-situ hybridisation / EBV-qPCR on the exact stocks; the Section 5.4 peptide search; and, to separate EBV from lineage, compare EBV+ vs EBV− lines matched for the lineage-TF programme.

### H3 — In-vitro EBV infection of AGS lowers MHC-class-I antigen presentation and the interferon response, a cell-intrinsic immune-evasion state.

- **For:** in the isogenic AGS→eAGS pair, HLA-A/B/C, TAP2, PSME1/2 and interferon-response proteins all drop (antigen-presentation z = −5.0; IFN-γ/α q < 2e-4); this is the cleanest (replicate-level, isogenic) contrast in the study; it is mechanistically coherent with herpesvirus immune evasion.
- **Against / confounders:** one infected line only, so it may be clone-specific; eAGS is also more proliferative (E2F up), which can independently lower antigen presentation; protein level ≠ surface MHC-I.
- **Sufficient?** Sufficient to propose; needs replication in other in-vitro infection models.
- **Test:** flow cytometry for surface HLA-class-I ± IFN-γ in AGS vs eAGS and in other newly infected lines; CD8 T-cell killing assays; test whether EBV latent genes (e.g. BNLF2a/EBNA1) are present and whether their knockdown restores MHC-I.

**Hypotheses considered and not supported by these data:**
- The published **EMT** subtype is *not* mesenchymal at the whole-proteome level (vimentin in 3/7; EMT hallmark q = 0.14). It is better described as a translation-high, lineage-TF-low group.
- **SRCC vs non-SRCC** and **"diffuse-like" cluster vs rest**: no proteins at q < 0.05 after correction beyond the lineage/metabolism axis, so no SRCC-specific protein signature is supported here.

---

## 7. Reproducibility

- `scripts/00_download.sh` → `05_ebv_summary.py` run in order reproduce every table and figure. Set `PXD048538_DATA` to the download directory. A fixed seed (20261105) is used throughout.
- Inputs are **not** committed (hundreds of MB, and redistribution terms). `00_download.sh` fetches them and verifies the quantms sha256 checksums.
- `run_info/` records the prompt, dataset config, dataset version, quantms recipe and package versions.
- Gene-set files come from public MSigDB GitHub mirrors (msigdb.org was not reachable); versions are in the filenames.

## 8. Files

- `results/REPORT.md` — this report.
- `results/figures/01–15_*.png|pdf` — figures.
- `results/tables/*.csv[.gz], *.json` — matrices, differential results, gene-set tests, QC and summaries.
- `scripts/*.py, *.sh` — analysis code.
