# PXD048538: biological-process and pathway analysis of the gastric-cancer cell-line proteome

**Prompt:** `0002` · **Run:** `runs/0002/20261107T120709Z_1d6e8ba6` · **Date:** 2026-11-07
**Data:** quantms reprocessing of PXD048538, version `e60b9e9aa4c8` (DIA-NN, PG.MaxLFQ; 14 gastric-cancer cell lines × 3 replicates = 42 whole-proteome DIA runs).
**Publication:** Hu C et al., *Proteome-based molecular subtyping and therapeutic target prediction in gastric cancer*, Mol Oncol 2024;18:1437–1459 (PMC11161736).

> **Scope.** prompt_0002 asks specifically for **biological processes and molecular pathways** affected by proteomic changes, and up to three **novel, testable** hypotheses. This report is pathway-centric. It is a companion to run 0001 (which covered subtyping, cell identity and an EBV analysis); where it overlaps it does so only as the entry point (the subtype contrast) and then goes to pathway structure and three new pathway-level hypotheses. Everything here is computational and **not experimentally validated**; "novel" = not found in a targeted Europe PMC search on 2026-11-07.

---

## 1. Summary

1. **Pathway activity across the 14 lines organises into four co-regulated modules** (Fig. 1): (1) a **proliferation/proteostasis** module (MYC targets, E2F, G2M, unfolded-protein response, DNA repair); (2) a **TGF-β / EMT / adhesion–secretory** module; (3) an **interferon / hypoxia** module; and (4) a broad **oxidative-lipid metabolism** module (OXPHOS, fatty-acid/bile-acid/xenobiotic metabolism, peroxisome, glutathione/ROS, glycolysis, estrogen response).
2. **Only the metabolism module (module 4) separates the published EMT vs Metabolism subtypes** (Welch p = 0.008); the proliferation/proteostasis module varies *independently* of subtype (Fig. 1B). So the dominant pathway axis between subtypes is metabolic, and proliferation is a second, orthogonal axis.
3. **Pathways associated with the subtype difference** (preranked GSEA, Fig. 2–3): the "EMT" subtype is high in **cytoplasmic translation / ribosome / NMD / mRNA-splicing / SRP-dependent targeting**; the "Metabolism" subtype is high in **oxidative phosphorylation, mitochondrial translation, TCA cycle, respiratory electron transport and fatty-acid metabolism**. This is a translation-versus-oxidative-metabolism proteome allocation difference.
4. **Three novel, testable pathway hypotheses** (Section 5), each supported and with the experiment that would test it:
   - **H1** — the Metabolism-subtype lines run a coordinated **glutathione / antioxidant / xenobiotic-detoxification programme** (GCLC/GCLM/GSS, GPX1/4, GSTs, ALDH/AKR, G6PD, ABCC1), co-regulated with OXPHOS (r = 0.77; Metab > EMT p = 0.001), confirmed by STRING enrichment (glutathione metabolism FDR 1e-22) → predicted chemo-/ferroptosis-resistance and a **GPX4 / GCLC dependency**.
   - **H2** — **MYC-target protein load is tightly coupled to the proteostasis (UPR) machinery** across the panel (r = 0.87, p = 5e-5) and to DNA repair (r = 0.67), independent of subtype → high-MYC GC lines carry a proteostatic load and a predicted **proteasome / HSP90 / UPR vulnerability**.
   - **H3** — the panel is **not organised along a Warburg glycolysis↔OXPHOS trade-off** (the two correlate positively, r = 0.50); instead a subset (NCC24, YCCEL1) is **OXPHOS-high** → candidate responders to OXPHOS/complex-I inhibition. (Most tentative; n = 14.)

---

## 2. Experimental design, prior knowledge and knowledge gaps

**Design (quantms SDRF + paper).** 14 lines (13 gastric adenocarcinoma + AZ521, duodenal), 3 biological replicates each, label-free DIA on a Q Exactive HF, quantified by DIA-NN (PG.MaxLFQ). This quantms version is **whole-proteome only** — no phospho/pTyr layers — so pathway inference here is abundance-based, not signalling-activity-based.

**What the paper established (treated as prior knowledge, not findings):** the EMT vs Metabolism subtypes; mTOR (EMT) and MEK (Metabolism) sensitivities; Src-family-kinase activity; outlier kinases ERBB2 (N87), FGFR2 (KatoIII), CSNK1D/E (MKN28). The paper's pathway statements come largely from the phospho/kinase layers; this report asks what the **whole-proteome** says about biological processes and pathways, and what new, testable pathway hypotheses follow.

**Knowledge gaps / limitations that bound every result below:**
- **n = 14 lines.** Cell-line means are the independent unit (replicates are cultures of one line). Pathway correlations use n = 14 and are descriptive, not powered for small effects.
- **Pathway activity = mean protein z-score** of set members; it is an abundance proxy, not a measured flux or signalling activity. Abundance of a metabolic enzyme need not equal pathway flux.
- **GSEA with proteome-wide programmes:** because the two subtypes differ along broad translation/metabolism axes, many gene sets reach significance (≈800 of 3,579 at q<0.1). The NES ranking and leading edge, not the raw count, carry the interpretation.
- **Confounding:** subtype, GI lineage and (for 3 lines) EBV status are partly confounded (see run 0001). Metabolic and lineage programmes cannot be fully separated here.
- **No external validation reachable:** DepMap/CCLE drug-sensitivity and dependency data (which would directly test the vulnerability predictions) and reactome.org were not reachable from the analysis environment; MSigDB/STRING were used from mirrors/API. The hypotheses are therefore *predictions*, explicitly awaiting those orthogonal datasets or wet-lab tests.

---

## 3. Methods (brief)

Code in `scripts/`, outputs in `results/`; fixed seed 20261105. Protein matrix built as in run 0001 (decoys/contaminants removed, run- and global-PG q ≤ 0.01, log2, median-aligned; present = ≥2/3 replicates; 8,152 protein groups, 4,004 in all 14 lines). Pathway activity = mean, across lines, of the z-scored protein abundances of each gene set's members (proteins quantified in ≥12 lines; 4,996 genes). Pathway **co-regulation modules** = hierarchical clustering (correlation distance) of the 34 most variable hallmark pathways, cut into 4 blocks. **Differential abundance** between subtypes = moderated t-test (empirical Bayes; `stats.py`). **Pathway enrichment** = preranked GSEA on the gene-level moderated-t statistic (2,000 permutations; Hallmark + KEGG + Reactome + GO-BP) and STRING functional enrichment for the detox gene set. Gene sets: MSigDB Hallmark 2023.1, KEGG/Reactome/GO 7.5.1 (GitHub mirrors).

---

## 4. Results: the pathway landscape

### 4.1 Four co-regulated pathway modules (`figures/01_pathway_coregulation.png`, `02_pathway_landscape.png`; `tables/pathway_modules_vs_subtype.csv`)

Correlating hallmark pathway activity across the 14 lines gives four modules:

| Module | Biology | Member hallmarks (abbrev.) | Separates EMT vs Metab? |
|---|---|---|---|
| 1 | **Proliferation / proteostasis** | MYC_TARGETS_V1/V2, E2F, G2M, UNFOLDED_PROTEIN_RESPONSE, DNA_REPAIR | No (Welch p = 0.24) |
| 2 | **TGF-β / EMT / adhesion–secretory** | TGF_BETA, EMT, KRAS_UP, PI3K_AKT_MTOR, COAGULATION, COMPLEMENT, APICAL_JUNCTION, APOPTOSIS | No (p = 0.60) |
| 3 | **Interferon / hypoxia** | IFN-α, IFN-γ, HYPOXIA, ALLOGRAFT_REJECTION | No (p = 0.93) |
| 4 | **Oxidative-lipid metabolism** | OXPHOS, FATTY_ACID, ADIPOGENESIS, BILE_ACID, XENOBIOTIC, PEROXISOME, CHOLESTEROL, GLYCOLYSIS, ROS, ESTROGEN, P53 | **Yes (p = 0.008)** |

The single most variable hallmark pathways across the panel are MYC_TARGETS_V2, OXIDATIVE_PHOSPHORYLATION, TGF_BETA, E2F, G2M and FATTY_ACID_METABOLISM (`tables/hallmark_variability.csv`). The least variable (near-uniform across lines) are inflammatory/TNFα/IL2-STAT5 — i.e. intrinsic inflammatory signalling is relatively constant.

The key structural point (Fig. 1B): **module 4 (metabolism) is the subtype axis, and module 1 (proliferation/proteostasis) is orthogonal to it** (module 1–module 4 r = −0.12, n.s.). The published "EMT/Metabolism" split is, at the whole-proteome level, chiefly a metabolic contrast; proliferation state varies within both subtypes.

### 4.2 Pathways associated with the subtype difference (`figures/03_gsea_subtype.png`, `04_gsea_leadingedge.png`; `tables/gsea_EMT_vs_Metabolism.csv`)

60 proteins differ at q < 0.05 (9 up in EMT, 51 up in Metabolism). Preranked GSEA on the full ranking gives a clean, coherent picture:

- **Higher in "EMT":** cytoplasmic translation, ribosome (KEGG), eukaryotic translation initiation/elongation, nonsense-mediated decay, rRNA processing, SRP-dependent co-translational targeting, mRNA splicing — a **protein-synthesis / translational-apparatus** programme (NES up to +2.9, q = 0.013).
- **Higher in "Metabolism":** oxidative phosphorylation, mitochondrial translation, aerobic/cellular respiration, respiratory electron transport + ATP synthesis, TCA cycle, fatty-acid metabolism, small-molecule catabolism — a **mitochondrial oxidative-metabolism** programme (NES down to −3.0, q = 0.013).

The leading-edge plot (Fig. 3) shows OXPHOS proteins concentrated at the Metabolism end and MYC targets (weakly) at the EMT end. So the two proteome subtypes differ primarily in how proteome mass is allocated between the **cytoplasmic translation machinery** and the **mitochondrial oxidative-metabolism machinery**.

---

## 5. Hypotheses

Each gives supporting evidence, confounders / what speaks against it, sufficiency, and a concrete test.

### H1 — The Metabolism-subtype lines run a coordinated glutathione/antioxidant + xenobiotic-detoxification programme, predicting redox-based chemoresistance and a GPX4/GCLC (ferroptosis) dependency.
(`figures/05_H1_detox_programme.png`; `tables/H1_string_enrichment_detox.csv`, `hypotheses_summary.json`)

- **For:** an 18-protein detox panel — glutathione synthesis (GCLC, GCLM, GSS), GSH peroxidases/transferases (GPX1, GPX4, GSTP1/K1/M3/O1, MGST1/2), aldehyde/ketone reduction (ALDH3A2, ALDH4A1, AKR1A1), NADPH/antioxidant (G6PD, TXNRD1, PRDX6) and the efflux pump ABCC1 — is coordinately **higher in the Metabolism subtype** (Welch p = 0.001) and **tightly co-regulated with OXPHOS** (r = 0.77, p = 0.001). STRING functional enrichment of this exact set returns glutathione metabolism (FDR 1.2e-22), glutathione conjugation (2.6e-21), the ROS pathway (9.6e-16) and cellular detoxification (1.5e-15) — an orthogonal confirmation that the co-abundant proteins form a bona-fide detox network. The presence of **GPX4 + GCLC/GCLM** (the glutathione–GPX4 axis) and **G6PD/TXNRD1** (NADPH supply) specifically marks a ferroptosis-defence configuration.
- **Against / confounders:** detox and OXPHOS co-vary with the GI-lineage/metabolism programme, so this may be one facet of a differentiated metabolic state rather than an independent axis; abundance ≠ enzymatic flux or actual GSH level; n = 14.
- **Sufficient?** Sufficient to propose a specific, mechanistic vulnerability; not sufficient to prove the dependency.
- **Test:** in Metabolism- vs EMT-subtype lines, measure sensitivity to the GPX4 inhibitor RSL3/ML210 and to erastin (system xc−), and to GCLC knockdown/BSO (glutathione depletion); predict Metabolism-subtype lines are **more resistant to ROS-inducing/platinum chemotherapy** but **more dependent on GPX4/GCLC** (collateral ferroptosis sensitivity). Cross-check GPX4/GCLC/SLC7A11 dependency scores for these lines in DepMap (needs depmap.org access).

### H2 — MYC-target protein abundance is tightly coupled to the proteostasis (UPR) and DNA-repair machinery across the panel, independent of subtype, predicting a proteostasis vulnerability in high-MYC lines.
(`figures/06_H2_H3_metabolic_coupling.png` A–B)

- **For:** across the 14 lines, the mean MYC-target (V1+V2) protein score correlates with the unfolded-protein-response score at **r = 0.87 (p = 5e-5)** and with DNA repair at r = 0.67 (p = 0.008). MYC activity varies independently of the metabolism subtype (high-MYC lines include both EMT — eAGS, AGS — and Metabolism — NCC24, KatoIII). This recapitulates, at the proteome level, the known biology that MYC-driven biosynthesis raises proteotoxic stress and engages the UPR/chaperone system.
- **Against / confounders:** MYC-target and UPR hallmark sets share some ribosome/chaperone members, which can inflate the correlation (partly definitional); high MYC-target abundance is a proxy for MYC activity, not a measurement of it; n = 14.
- **Sufficient?** Sufficient to nominate a vulnerability; the shared-member caveat must be addressed before causal claims.
- **Test:** after removing overlapping genes between the MYC and UPR sets, re-test the correlation; then test whether MYC-high lines (NCC24, KatoIII, eAGS, AGS) are preferentially sensitive to proteasome inhibition (bortezomib), HSP90 inhibition, or IRE1/PERK (UPR) inhibition; MYC dependency in DepMap for the same lines.

### H3 — The GC panel is not organised along a Warburg glycolysis↔OXPHOS trade-off; an OXPHOS-high subset (NCC24, YCCEL1) is a candidate for OXPHOS-targeted therapy.
(`figures/06_H2_H3_metabolic_coupling.png` C)

- **For:** glycolysis and OXPHOS pathway activity are **positively**, not reciprocally, correlated across the panel (r = 0.50, p = 0.07) — there is no anti-correlated Warburg axis at the protein level. OXPHOS activity instead has a high-tail: NCC24 (z = 0.9) and YCCEL1 (z = 0.5) stand out as OXPHOS-high.
- **Against / confounders:** this is the weakest of the three (the glycolysis–OXPHOS correlation is only borderline, p = 0.07, and could be driven by overall "metabolic-high" lineage lines; OXPHOS-protein abundance ≠ respiration rate); n = 14.
- **Sufficient?** Insufficient on its own — it is a lead, not a conclusion. State: the evidence is suggestive but not statistically firm.
- **Test:** Seahorse OCR/ECAR to confirm NCC24/YCCEL1 are oxidative; then test complex-I inhibitor (IACS-010759) sensitivity; DepMap OXPHOS-gene dependency for these lines.

**Considered but not supported:** a proliferation↔oxidative-metabolism trade-off (modules 1 and 4 are statistically independent, r = −0.12, n.s.), so the data do **not** support a simple "proliferative vs metabolic" dichotomy across these lines.

---

## 6. Reproducibility & files

Run `scripts/00_download.sh` (from run 0001, same data) then `01_build_matrix.py → 04_hypotheses.py` with `PXD048538_DATA` set to the data directory. `03_enrichment.py` runs 2,000-permutation GSEA (~9 min). Inputs are not committed; they are fetched and checksum-verified by the download script. `run_info/` holds the prompt, dataset config and package versions.

- `results/REPORT.md` — this report; `results/FIGURES_AND_TABLES.md` — index.
- `results/figures/01–06_*.png|pdf`; `results/tables/*.csv[.gz], *.json`.
- `scripts/*.py` (+ `common.py`, `stats.py` with added preranked GSEA).
