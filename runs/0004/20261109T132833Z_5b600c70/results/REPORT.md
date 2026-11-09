# PXD048538: protein co-abundance, interaction networks and complexes in gastric-cancer cell lines

**Prompt:** `0004` · **Run:** `runs/0004/20261109T132833Z_5b600c70` · **Date:** 2026-11-09
**Data:** quantms reprocessing of PXD048538, version `e60b9e9aa4c8` (whole-proteome DIA; 14 gastric-cancer cell lines × 3 replicates). **Publication:** Hu C et al., Mol Oncol 2024;18:1437–1459 (PMC11161736).

> **Scope.** prompt_0004 asks to infer protein–protein interactions, networks and complexes from **protein-abundance patterns** across the cell lines, and to relate them to the molecular subtypes. This is an abundance-based (guilt-by-co-regulation) analysis on the whole proteome — the layer the quantms `e60b9e9aa4c8` release *does* contain. All results are computational and **not experimentally validated**; co-abundance indicates co-regulation, **not necessarily physical interaction**. "Novel" means "no prior evidence in GO-CC complexes or STRING (any channel)" on 2026-11-09.

---

## 1. Summary

1. **Co-abundance across the 14 lines reliably recovers known complexes at high confidence.** Biweight-midcorrelation of protein abundance separates same-complex pairs from random pairs (median bicor 0.16 vs 0.02, Mann-Whitney p≈0), but the *average* signal is modest (AUC 0.61) — expected for n=14. Crucially, the **strongest edges are highly reliable**: at bicor ≥ 0.9, pairs are **45× enriched** for known-complex membership (precision 0.22 vs 0.005 baseline; Fig. 1). Textbook complexes are recovered cleanly (ATP synthase, condensin SMC2–SMC4, MICOS CHCHD3–IMMT, eIF3, OST DAD1–DDOST, proteasome).
2. **The co-abundance network decomposes into 31 modules** that map to real complexes/processes (respiratory chain, NADH dehydrogenase, condensin, ESCRT-I, RNA-Pol-II, exosome, DNA-replication pre-initiation; Fig. 2) — but these modules vary **line-to-line, not by subtype** (0 modules subtype-associated at q<0.1).
3. **Complex *abundance* is, however, coordinated with the subtype** (Fig. 3C): OXPHOS/TCA complexes (respiratory chain IV, cytochrome complex, TCA-cycle enzyme complex) are coherently higher in **Metabolism** lines, while gene-expression / protein-turnover complexes (CDC73–PAF1, transcription-elongation factor, SCF ubiquitin-ligase, eIF3, exon-junction, SUMO-ligase) are higher in **EMT** lines (8 complexes q<0.1). This connects the complex layer to the subtypes found in runs 0001–0003.
4. **62% of the strongest co-abundance pairs are already known** (GO-CC complex or STRING ≥0.4); the remaining **332 "candidate novel" pairs** have no prior interaction evidence — a set that **cannot be taken as physical interactions** at n=14, but whose top members are biologically coherent (several ER/mitochondrial membrane proteins: MTCH2–VAPA, ERP29–TMED2, ROMO1–TMEM33, AGPAT5–MLEC).
5. **Three hypotheses** (Section 5): **H1** candidate novel co-regulation at the ER–mitochondria membrane interface; **H2** subtype-coordinated complex abundance (OXPHOS complexes in Metabolism, gene-expression complexes in EMT); **H3** (methodological) co-abundance here is a reliable *known-complex map* but underpowered for de-novo interaction discovery.

---

## 2. Design, prior knowledge and knowledge gaps

**Input.** Whole-proteome matrix (8,152 protein groups; 4,004 quantified in all 14 lines), collapsed to 4,002 genes present in every line (the unit for co-abundance). **Method.** Protein co-abundance = biweight midcorrelation (robust to outlier lines) across the 14 cell-line means. Modules = hierarchical clustering on 1−bicor of the 2,000 most-variable genes. **References** (all reachable without the blocked bulk/FTP hosts): GO cellular-component "complex" terms (benchmark), the **STRING** REST API (PPI evidence channels), and the **ComplexPortal** API (EBI) for spot checks.

**What the paper established (prior knowledge, not findings here):** the EMT/Metabolism subtypes and their molecular features. The paper did **not** build a co-abundance interaction network; this report does, and asks what protein-complex structure the abundance data encode and how it relates to the subtypes.

**Knowledge gaps / limitations that bound every result:**
- **n = 14 is small for correlation-based interaction inference.** A correlation needs |r|≳0.66 for p<0.01 at n=14, and co-abundance has many false positives on average (AUC 0.61). The analysis is therefore restricted to the **strongest edges** (bicor ≥ 0.9), where enrichment for real complexes is 45×. Even there, co-abundance ≠ physical interaction.
- **Co-regulation is not interaction.** Two proteins can co-vary because they are in one complex, in one pathway, or co-induced by a shared transcriptional/stability program. "Candidate novel" pairs (no STRING/GO-CC evidence) may be indirect.
- **Reference incompleteness cuts both ways.** Benchmarking on GO-CC "complex" terms alone over-counts "novel" (e.g. ribosome/spectrin pairs); STRING annotation corrects most of this (ribosomal, prohibitin, etc. reclassified as known). Remaining novel pairs are limited by STRING's own coverage.
- **Replicates are not independent;** line means are used, so n = 14 biological units. No DepMap/CCLE access to cross-check.

---

## 3. Methods (brief)
`01_build_matrix.py` builds the protein matrix (as runs 0002–0003). `02_coabundance.py` computes the gene×gene biweight-midcorrelation and benchmarks same-GO-CC-complex pairs vs random. `03_modules.py` clusters the top-variable genes into modules, annotates them (hypergeometric enrichment vs GO-BP/CC, Reactome, KEGG), computes a module eigengene and its subtype association, and lists strong (bicor≥0.9) pairs. `04_hypotheses.py` annotates strong pairs with STRING combined/physical/experimental scores (one network call over all 836 genes; input symbols mapped to STRING preferred names to handle aliases), tests complex-level subtype association, and assembles the figures.

---

## 4. Results

### 4.1 Co-abundance recovers complexes — but only at the top edges (`figures/01_coabundance_benchmark.png`; `tables/coabundance_benchmark.json`)
Same-complex pairs are shifted to higher co-abundance than random (medians 0.16 vs 0.02), yet the distributions overlap heavily (AUC 0.61). The precision of "edge is within a known complex" rises steeply with the threshold, reaching **0.22 at bicor ≥ 0.9 (45× the 0.005 baseline)**. The proteasome, ATP synthase and other obligate complexes appear as dense co-abundance blocks. **Takeaway:** only the strongest edges are trustworthy here — which is how the rest of the analysis is scoped.

### 4.2 Modules are line-specific, not subtype-specific (`figures/02_module_landscape.png`; `tables/coabundance_modules.csv`)
31 modules (1,986 genes) map to recognisable complexes/processes (respiratory chain, NADH dehydrogenase, condensin, ESCRT-I, RNA-Pol-II, exosome/exoribonuclease, DNA-replication pre-initiation, histone-deacetylase, SCF/ubiquitin-ligase). Their eigengenes vary strongly **between individual lines** (e.g. condensin module high in MKN28/N87) but **no module separates the EMT/Metabolism subtypes** (q>0.1 for all). So whole co-regulated modules are not the substrate of the subtype distinction.

### 4.3 Complex abundance *is* coordinated with subtype (`figures/03_ppi_validation_subtype.png` C; `tables/complex_subtype_association.csv`)
Testing the mean abundance of each curated complex against subtype gives **8 complexes at q<0.1**:
- **Higher in Metabolism:** TCA-cycle enzyme complex, respiratory-chain complex IV, cytochrome complex, mitochondrial TCA complex — i.e. the **OXPHOS/TCA machinery as coherent units** (within-complex cohesion bicor 0.46–0.52).
- **Higher in EMT:** CDC73–PAF1, transcription-elongation-factor complex, SCF ubiquitin-ligase, exon–exon-junction complex, eIF3, SUMO-ligase — **gene-expression and protein-turnover complexes**.

This recapitulates, at the level of whole complexes, the proteome-allocation difference seen in runs 0002–0003 (Metabolism = oxidative metabolism; EMT = translation/transcription), and shows the subunits move **together**, consistent with coordinated complex assembly rather than scattered single-protein changes.

### 4.4 Most strong pairs are known; a residual novel set (`figures/03_ppi_validation_subtype.png` A–B; `tables/strong_pairs_string_annotated.csv.gz`)
Of 877 strong pairs (bicor ≥ 0.9): **283 share a GO-CC complex, 262 have STRING support (≥0.4), 332 are candidate novel** (STRING <0.15 and not GO-CC). 62% are thus already known — a strong positive control for the method. The candidate-novel top (all STRING score 0, verified by direct query) includes **MTCH2–VAPA, ERP29–TMED2, ROMO1–TMEM33, AGPAT5–MLEC, BPNT2–ROMO1, IMPDH1–PTGES2, GIGYF2–TJP1**.

---

## 5. Hypotheses

### H1 — The strongest co-regulated protein pairs with no prior interaction evidence are enriched for ER/secretory-membrane and ER–Golgi trafficking proteins, marking a candidate, previously unannotated secretory-membrane co-regulation set.
- **For:** the candidate-novel-pair genes (bicor ≥ 0.9, no GO-CC/STRING evidence; STRING score 0 confirmed by direct query) are significantly enriched, **relative to all strong-pair genes**, for the GO-CC *nuclear-outer-membrane/ER membrane network* (odds ratio 2.1, p = 0.002), the *ER–Golgi intermediate-compartment membrane* (p = 0.004) and *vesicle / ER-to-Golgi transport-vesicle membrane* (p = 0.03–0.05) (`tables/candidate_novel_gene_enrichment.csv`). Concrete top pairs: ERP29–TMED2 (ER chaperone + ER–Golgi cargo receptor), AGPAT5–MLEC (lipid acyltransferase + ER malectin), and MTCH2–VAPA (bicor 0.97) at the mitochondria–ER interface.
- **Against / confounders:** co-abundance ≠ physical interaction — these pairs may reflect a shared secretory-membrane biogenesis program rather than complexes; the compartment enrichment is nominal (not FDR-significant across all GO-CC terms); STRING/GO coverage gaps make "novel" provisional; n = 14.
- **Sufficient?** Not sufficient to assert interactions — these are **candidate co-regulated pairs**, explicitly flagged for validation.
- **Test:** co-immunoprecipitation / proximity labelling (BioID, APEX) for ERP29–TMED2 and MTCH2–VAPA; co-fractionation (size-exclusion / BN-PAGE co-elution) to test physical association; confirm the co-regulation in an independent cell-line proteomic panel.

### H2 — The EMT/Metabolism subtypes differ in the abundance of whole protein complexes: OXPHOS/TCA complexes are coordinately high in Metabolism, gene-expression/turnover complexes high in EMT.
- **For:** 8 curated complexes differ by subtype at q<0.1, with high within-complex co-abundance cohesion (the subunits move together). The direction — respiratory-chain IV / cytochrome / TCA-enzyme complexes up in Metabolism; CDC73–PAF1, transcription-elongation, SCF, eIF3 up in EMT — matches the pathway/altitude differences from runs 0002–0003, now shown at complex resolution.
- **Against / confounders:** complex "abundance" here is the mean of subunit abundances, not measured assembled-complex levels; subtype and lineage are confounded; n=14.
- **Sufficient?** Sufficient to propose coordinated complex-abundance differences; not to prove assembled-complex stoichiometry.
- **Test:** native/BN-PAGE or complexome profiling of a Metabolism vs EMT line pair to confirm assembled OXPHOS-complex levels; test whether Metabolism lines have higher respiratory capacity (Seahorse) tracking the complex-abundance signal.

### H3 — (methodological) Protein co-abundance over this 14-line panel is a reliable map of *known* complexes but is underpowered for de-novo interaction discovery.
- **For:** strong edges are 45× enriched for known complexes and 62% of them are already annotated; the AUC is only 0.61 overall; "novel" pairs cannot be separated from indirect co-regulation at n=14.
- **Against / confounders:** more sensitive estimators (partial correlation, more conditions) could extract more; absence of evidence (STRING) is not evidence of absence.
- **Sufficient?** As a caution and a design recommendation, yes.
- **Test:** repeat on a larger GC cell-line panel (dozens–hundreds of conditions, e.g. CCLE/ProteomeHD-scale) where co-abundance resolves novel complexes; use orthogonal co-fractionation MS to call physical complexes directly.

**Considered but not supported:** that co-regulated protein *modules* define the molecular subtypes (Section 4.2 — no module is subtype-associated); the subtype signal is in complex *abundance*, not module co-regulation structure.

---

## 6. Reproducibility & files
`scripts/00_download.sh` fetches the data and gene sets; run `01_build_matrix → 02_coabundance → 03_modules → 04_hypotheses` with `PXD048538_DATA` set. `04` queries the STRING API (string-db.org). Inputs not committed.
- `results/REPORT.md`, `results/FIGURES_AND_TABLES.md`; `results/figures/01–03_*.png|pdf`; `results/tables/*` (co-abundance matrix, modules, strong pairs with STRING annotation, complex–subtype table, summaries).
