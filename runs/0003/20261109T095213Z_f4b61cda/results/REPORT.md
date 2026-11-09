# PXD048538: signaling-pathway and PTM analysis of the gastric-cancer cell-line (phospho)proteome

**Prompt:** `0003` · **Run:** `runs/0003/20261109T095213Z_f4b61cda` · **Date:** 2026-11-09
**Publication:** Hu C et al., *Proteome-based molecular subtyping and therapeutic target prediction in gastric cancer*, Mol Oncol 2024;18:1437–1459 (PMC11161736).

> **Read this first — where the PTM data come from.** The quantms release of PXD048538 that this project targets (**version `e60b9e9aa4c8`**) contains **only the whole proteome** — it has **no phosphorylation or other PTM layer** (confirmed in runs 0001–0002: the `.pg`/`.feature` files are label-free protein/peptide quant with carbamidomethyl/oxidation/acetyl only). A signaling/PTM analysis therefore cannot be done on that release alone. The phosphorylation data used here are the **author-processed phosphoproteome (PP: pS/pT/pY) and tyrosine-phosphorylation (YP: pY by immunoaffinity) matrices from the study's Supplementary Table S3** (`MOL2-18-1437-s005.xls`; log2, quantile-normalised, imputed, batch-corrected). This is the same study, re-analysed independently. Every signaling claim below rests on that supplementary PTM data, and the report says so at each step. All results are computational and **not experimentally validated**.

---

## 1. Summary

1. **The subtype difference in the bulk Ser/Thr phosphoproteome is abundance-driven.** 20 S/T/Y phosphosites differ between the EMT and Metabolism proteome subtypes at q<0.05, but after **normalising each site to its parent protein, 0 sites remain significant** (Fig. 2A). So most "phospho differences" between subtypes reflect protein-abundance remodeling, not a change in phosphorylation stoichiometry. This is a limitation to respect, not a finding to force.
2. **The robust, abundance-independent signaling difference is at the tyrosine level.** 25 tyrosine sites differ (q<0.05), and **all but one are higher in the Metabolism subtype** (Fig. 2B). They concentrate on a **receptor-tyrosine-kinase → adaptor → RAS-ERK module**: EPHB2 pY780, SHC1 pY213, GAREM1 pY453, SHB pY96, RASA1 pY438, ITGB4 pY1207/pY920, PEBP1(RKIP) pY181 (Fig. 4A).
3. **Motif-based kinase-activity inference (phosX / Kinase Library) gives a clean subtype dichotomy** that independently recovers — and mechanistically explains — the paper's drug sensitivities (Fig. 3):
   - **Metabolism subtype** = high **receptor-tyrosine-kinase** activity (MET, FGFR2, FGFR3, MST1R/RON, JAK3, DDR2 at FDR<0.1; also Eph receptors, HER2) and high **MAP2K/MAP3K** Ser/Thr activity (MEK5, MEK2, ZAK, MEKK2, KHS1/2). → an **RTK→RAS-ERK** program (the paper reports this subgroup is MEK-inhibitor sensitive).
   - **EMT subtype** = high **AGC-kinase** Ser/Thr activity (AKT1/2/3, SGK1/3, P70S6K, RSK3) → a **PI3K-AKT-mTOR/SGK** program (the paper reports this subgroup is mTOR-inhibitor sensitive).
4. **Three testable signaling hypotheses** (Section 5): **H1** the Metabolism subtype is driven by a convergent multi-RTK→RAS-ERK network (MEK dependency + RTK redundancy); **H2** the EMT subtype depends on PI3K-AKT-mTOR with an under-appreciated **SGK** arm; **H3** the actionable signaling difference between subtypes is carried by tyrosine (RTK) signaling, not the bulk S/T phosphoproteome.

---

## 2. Design, prior knowledge and knowledge gaps

**Phospho data (Supplementary Table S3).** Global phosphoproteome **PP**: 8,452 sites (7,406 pS, 887 pT, 159 pY) across 14 lines × 3 replicates; tyrosine **YP** (anti-pY immunoaffinity): 4,780 pY sites × 14 lines × 2 replicates. Both are author-processed. The whole proteome (**WP**, same workbook) is used for protein-normalisation. Cell-line means are the independent unit (n = 14).

**What the paper established (prior knowledge, not findings here):** the EMT/Metabolism subtypes; a kinome landscape; mTOR sensitivity in the EMT subgroup and MEK sensitivity in the Metabolism subgroup (validated there with rapamycin/trametinib); Src-family-kinase activity; outlier kinases ERBB2 (N87), FGFR2 (KatoIII), CSNK1D/E (MKN28). This report asks what an **independent, motif-based** re-analysis of the phospho layers says about the *kinases and signaling wiring* behind the subtypes, and whether site-level evidence supports new, testable mechanisms.

**Knowledge gaps / limitations that bound every result:**
- **Not the target release.** The PTM layer is not in quantms `e60b9e9aa4c8`; it is supplementary data from the paper. Reprocessing the raw phospho runs (DIA PP, DDA pY) through a uniform pipeline was out of scope (DIA phospho needs a spectral library; very heavy). Results therefore inherit the authors' processing choices (normalisation, imputation, batch correction).
- **Abundance confounding.** Phosphosite intensity = occupancy × protein abundance. Protein-normalisation removes the abundance term but costs power (two noisy measurements) and coverage (5,731 of 8,452 S/T sites have a matched protein). The 0 protein-normalised significant sites means "no abundance-independent S/T difference is **detectable** at this n/FDR", not "none exists".
- **Kinase inference is motif-based.** phosX scores the ±5 flanking sequence of each site against the Ser/Thr (Johnson 2023) and Tyr (Yaron-Barir 2024) Kinase-Library PSSMs and tests rank enrichment of the differential statistic. It infers *which kinases' sequence-motifs are coordinately up/down*; it does not measure kinase activity directly, cannot distinguish kinases with near-identical motifs, and used 200 permutations (p floored at ~0.005; kinases are therefore ranked by the continuous KS enrichment statistic, with FDR shown).
- **Tyr seqrnk uses raw pY fold-changes** (IP-based YP is only weakly abundance-coupled); S/T seqrnk uses raw PP fold-changes, so the S/T kinase inference carries the same abundance caveat as the S/T differential (hence H3).
- **n = 14, confounded annotations** (subtype/lineage/EBV), as in runs 0001–0002. No DepMap/CCLE access to cross-check dependencies.

---

## 3. Methods (brief)

Code in `scripts/`. `phospho_common.py` loads the PP/YP/WP matrices, maps sample names to lines, and builds the 10-mer flanking window for each site (phospho-residue at position 6) from the UniProt human Swiss-Prot FASTA; it also fixes the Excel gene→date corruption (SEPT/MARCH). `01_phospho_build.py` builds line-mean PP/YP matrices and protein-normalised PP (site − parent protein). `02_diff_phospho.py` runs moderated-t differential analysis (empirical Bayes) per contrast on PP (raw + protein-normalised) and YP, and writes `seqrnk` files. `03_kinase_activity.py` parses the phosX outputs (one run per contrast × {Ser/Thr, Tyr}, 200 permutations, built-in Kinase-Library PSSMs). `04_hypotheses.py` assembles the signaling-module figure. Contrasts: **EMT vs Metabolism** (primary) and **natural-EBV+ vs EBV−** (secondary).

---

## 4. Results

### 4.1 Coverage (`figures/01_phospho_coverage.png`)
Phosphosite coverage is even across lines (no line is an outlier); the global PP layer is 88% pSer, 10% pThr, 2% pTyr, and the dedicated pY-IP layer adds 4,780 tyrosine sites — the two tyrosine sources are the signaling-rich part of the data.

### 4.2 Differential phosphorylation: abundance-driven at S/T, real at tyrosine (`figures/02_diff_phospho_volcano.png`; `tables/diff_{PP,YP}_*.csv.gz`)
EMT vs Metabolism: **PP raw** 20 sites at q<0.05; **PP protein-normalised 0 sites** (Fig. 2A is a symmetric null cloud). **YP** 25 sites at q<0.05, 24 of them higher in Metabolism (Fig. 2B). The significant pY sites map to a coherent signaling set — **EPHB2, SHC1, GAREM1, SHB, RASA1, ITGB4, SASH1, SDCBP, PEBP1/RKIP** — i.e. receptor tyrosine kinases and the Grb2/SHC-class adaptors and RasGAP that couple them to RAS-ERK (Fig. 4A). PEBP1(RKIP) pY181, a Src site that de-represses RAF, being higher in Metabolism fits an active RAF-MEK-ERK axis there.

### 4.3 Kinase-activity landscape (`figures/03_kinase_activity_subtype.png`, `04_signaling_modules.png`; `tables/kinase/*`)
phosX on the subtype contrast gives a clean, two-sided dichotomy:

| Layer | Higher in **EMT** | Higher in **Metabolism** |
|---|---|---|
| **Tyr kinome** | LCK, WEE1 (weak) | **MET, FGFR2, FGFR3, MST1R, JAK3, DDR2** (FDR<0.1); EPHB1/3, EPHA3/8, HER2, TRKA/B, CSF1R |
| **Ser/Thr kinome** | **AKT1/2/3, SGK1/3, P70S6K, RSK3** (AGC; AKT3/CAMK1D/RSK3 FDR<0.1) | **MEK5, MEK2, ZAK, MEKK2, KHS1/2, GCK** (MAP2K/MAP3K; GCK/ZAK/KHS1-2 FDR<0.1) |

Two points. First, the Metabolism subtype shows **both** elevated RTK (tyrosine) and elevated MAP2K/MAP3K (Ser/Thr) activity — a complete RTK→RAS→MEK-ERK cascade, matching the site-level adaptor pY in 4.2. Second, this independently reproduces the paper's pharmacology (EMT→mTOR via AKT/SGK/P70S6K; Metabolism→MEK via the MAP2K module) from a fully independent, motif-based method, and extends it upstream to the receptor layer. The S/T inference carries the abundance caveat (H3); the tyrosine inference does not, and is where the two subgroups differ most sharply.

### 4.4 EBV+ vs EBV− (secondary, underpowered)
Only 6 PP and 5 YP sites differ at q<0.05 between the 3 natural-EBV+ and the EBV− lines, so kinase inference for this contrast (EBV+ leaning to EPHB1/FER/FES/TEC tyrosine kinases) rests on very few differential sites and is **not** advanced as a hypothesis. Tables are provided (`diff_*_EBV_vs_nonEBV`, `kinase/EBV_*`) but the evidence is insufficient.

---

## 5. Hypotheses

### H1 — The Metabolism/GI-lineage subtype is driven by a convergent, redundant receptor-tyrosine-kinase network feeding RAS-ERK; MEK is the actionable node.
- **For:** 24/25 differential pY sites are higher in Metabolism and sit on an RTK→adaptor→RasGAP module (EPHB2, SHC1 Y213, GAREM1 Y453, SHB Y96, RASA1 Y438, ITGB4; Fig. 4A); phosX infers **multiple RTKs co-active** in Metabolism (MET, FGFR2/3, MST1R, JAK3, DDR2 at FDR<0.1, plus Eph/HER2/TRK) and elevated **MAP2K/MAP3K** (MEK5/2, MEKK2, KHS1/2); the paper independently found this subgroup MEK-inhibitor (trametinib/selumetinib) sensitive — H1 supplies the upstream mechanism.
- **Against / confounders:** phosX cannot separate kinases with similar motifs (the "RTK network" is a motif-family signal, not proof each receptor is active); subtype and GI lineage are confounded; no direct pERK measurement here.
- **Sufficient?** Sufficient to propose a convergent-RTK/MEK-node model and the redundancy prediction; not to prove any single RTK drives it.
- **Test:** phospho-ERK/MEK immunoblot across the panel; sensitivity to MEK inhibition (should recapitulate the paper) vs single-RTK inhibitors (predict single-RTK resistance due to redundancy, MEK/ERK-node or multi-RTK combinations more effective); knock-down of SHC1/GAREM1 (the convergent adaptors); DepMap MET/FGFR2/MAP2K1 dependency for the Metabolism lines.

### H2 — The EMT subtype runs a PI3K-AKT-mTOR program with a co-active, under-appreciated SGK arm.
- **For:** phosX Ser/Thr activity higher in EMT for AKT1/2/3 **and** SGK1/SGK3, plus the downstream P70S6K and RSK3 (Fig. 3, 4C); AKT and SGK share the basophilic RxRxxS/T motif and both feed mTORC1 — the paper's mTOR (rapamycin) sensitivity of this subgroup is consistent, and the SGK arm predicts that AKT inhibition alone may be bypassed (the paper indeed found a pan-AKT inhibitor *ineffective* in the EMT subgroup despite PI3K-AKT enrichment).
- **Against / confounders:** motif overlap means "AKT vs SGK" cannot be fully resolved by sequence; the S/T inference is abundance-caveated (H3); n = 14.
- **Sufficient?** Sufficient to nominate SGK as a parallel node explaining AKT-inhibitor insensitivity; needs direct test.
- **Test:** pAKT-S473 / pNDRG1 (SGK substrate) immunoblot in EMT vs Metabolism lines; compare mTOR, pan-AKT, and SGK (e.g. GSK650394) inhibitor sensitivity in the EMT subgroup; the prediction is that SGK or mTOR-node inhibition outperforms AKT inhibition.

### H3 — The actionable signaling difference between GC proteome subtypes is tyrosine-level; the bulk Ser/Thr phosphoproteome difference is abundance-coupled.
- **For:** 0 of 8,452 S/T/Y sites are differentially phosphorylated between subtypes after protein-normalisation, whereas 25 tyrosine (IP) sites are, with a coherent RTK signature; the tyrosine layer (weakly abundance-coupled) carries the robust signal, the S/T layer does not.
- **Against / confounders:** protein-normalisation loses power and coverage, so a true abundance-independent S/T difference could be present but undetected at this n; the claim is about detectability, not non-existence.
- **Sufficient?** Sufficient as a methodological caution and a prioritisation (focus signaling-based stratification/targeting of these subtypes on the tyrosine/RTK axis); not a biological mechanism on its own.
- **Test:** deeper, protein-normalised phosphoproteomics (more replicates) to set a real detection limit for abundance-independent S/T differences; targeted pY validation of the RTK-adaptor module by PRM.

**Considered but not supported:** an EBV-specific signaling program (Section 4.4 — too few differential sites, insufficient evidence).

---

## 6. Reproducibility & files
`scripts/00_download.sh` fetches the quantms data, the paper supplement (Table S3), gene sets and the UniProt FASTA; `pip install phosx` provides the Kinase-Library PSSMs. Run `01_phospho_build → 02_diff_phospho → run_phosx (phosX per contrast×layer) → 03_kinase_activity → 04_hypotheses`. phosX runs are the slow step (~15 min per Ser/Thr contrast). Inputs not committed.
- `results/REPORT.md`, `results/FIGURES_AND_TABLES.md`.
- `results/figures/01–04_*.png|pdf`; `results/tables/*` (phospho matrices, differential, `kinase/` activity tables, seqrnk files); `run_info/`.
