# PXD048538: integrated synthesis of four proteomic analyses of gastric-cancer cell lines

**Prompt:** `0005` (synthesis) · **Run:** `runs/0005/20261109T145940Z_51591012` · **Date:** 2026-11-09
**Integrates:** run 0001 (exploration / subtyping / EBV), run 0002 (pathways), run 0003 (signaling / PTM), run 0004 (protein–protein interaction / complexes). **Dataset:** quantms PXD048538 `e60b9e9aa4c8` (14 GC cell lines). **Publication:** Hu C et al., Mol Oncol 2024;18:1437–1459.

> **The single most important caveat, stated first.** All four analyses were run on **one dataset (14 cell lines)**. When several analyses "agree", that is **internal consistency of correlated features computed from the same protein abundances — not independent experimental validation.** Many of the signatures below are literally re-derived from the same matrix (oxidative-metabolism pathways, OXPHOS complexes and the metabolic-subtype kinase ranking all come from the same ~8,000 protein abundances). The quantitative "convergence" (Fig. 1B) should be read as *one finding seen through several lenses*, and every hypothesis still needs orthogonal, wet-lab validation. Nothing here is experimentally validated.

---

## 1. Executive summary

Across the four analyses, the gastric-cancer cell-line proteome is organised by **two largely orthogonal axes**:

1. **A dominant GI-lineage / oxidative-metabolism axis = the paper's EMT vs Metabolism subtype.** The "Metabolism" lines (KatoIII, NCC24, Ocum1, SNU601, SNU719, YCC3, YCCEL1) are a **well-differentiated, HNF-transcription-factor-driven, oxidatively-metabolic** state; the "EMT" lines (AGS, eAGS, AZ521, MKN28, N87, SNU1750, TMK1) are **translation/biosynthesis-high and lineage-TF-low**. Four independent analytical angles land on this one axis (Fig. 1): GI-lineage TFs HNF1A/HNF4A/FOXA2 (run 0001), OXPHOS/fatty-acid pathways and the glutathione detox programme (run 0002), the RTK→MEK vs AKT-mTOR signalling split (run 0003), and OXPHOS-complex abundance (run 0004). These four metabolic signatures co-vary across the 14 lines at **mean r = 0.66**.
2. **An orthogonal MYC-target / proteostasis / proliferation axis** (run 0002/0003) that varies **within** both subtypes (r = −0.24 with the metabolic axis) — a second, independent source of heterogeneity.

The strongest, most actionable integrated ideas are: the Metabolism subtype's **RTK→MEK dependency** and **ferroptosis (GPX4/GCLC) vulnerability**; the EMT subtype's **PI3K-AKT-mTOR/SGK dependency**; and, as a clean mechanistic aside, **EBV-induced MHC-I/interferon suppression** in the isogenic AGS→eAGS pair. One earlier claim was **rejected on re-analysis** (TMK1 as a cryptic EBV carrier — no EBV peptides found), illustrating that the project did discriminate, not just accumulate, positive results.

---

## 2. The integrated biological picture

### 2.1 Axis 1 — GI-lineage / oxidative metabolism (the subtype) (`figures/01_cross_run_scorecard.png`)
Per cell line, eight derived signatures were recomputed on the shared protein matrix (`tables/per_line_axis_scores.csv`). Everything that marks the **Metabolism** subtype points the same way and is statistically associated with it:

| Signature (source run) | EMT − Metab | Welch p |
|---|---|---|
| GI-lineage TF HNF1A/HNF4A/FOXA2 (r1) | −1.67 | <0.001 |
| Oxidative metabolism OXPHOS+FAO (r2,r4) | −0.58 | 0.002 |
| Glutathione/antioxidant detox (r2) | −0.74 | 0.006 |
| OXPHOS complexes (resp chain/cytochrome) (r4) | −0.55 | 0.042 |
| RTK-adaptor pY module EPHB2/SHC1/RASA1 (r3) | −0.46 | 0.071 |
| Cytoplasmic translation / ribosome (r2,r3) | **+0.75** | 0.015 |

The first four form a tight block (mean r = 0.66; Fig. 1B): a coherent **differentiated, HNF-driven, mitochondrial-oxidative** programme that also carries antioxidant/detox capacity and receptor-tyrosine-kinase input. The EMT subtype is its mirror: high cytoplasmic-translation/biosynthetic machinery, low lineage TFs. This matches — and sharpens with molecular detail — the authors' observation that their Metabolism cluster aligns with intestinal-type / G-Int cell lines.

### 2.2 Axis 2 — MYC / proteostasis / proliferation (orthogonal)
MYC-target abundance, the unfolded-protein-response machinery and E2F/DNA-repair move together (run 0002: MYC–UPR r = 0.87) but are **independent of the metabolic subtype** (r = −0.24; Fig. 1B). High-MYC lines span both subtypes (NCC24, KatoIII, eAGS, AGS). This is a genuinely separate axis of heterogeneity and a separate therapeutic angle (proteostasis stress).

### 2.3 Signalling wiring on the subtype axis (run 0003)
The subtype is not just metabolic — it has a signalling identity. Motif-based kinase-activity inference independently reproduced the paper's pharmacology: **Metabolism = RTK (MET/FGFR2-3/RON/Eph/HER2) → RAS-ERK → MEK** (the paper: MEK-inhibitor sensitive); **EMT = PI3K-AKT-mTOR with an SGK arm** (the paper: mTOR-sensitive, pan-AKT ineffective). Critically, run 0003 also showed the **bulk Ser/Thr phosphoproteome difference between subtypes is abundance-driven** (0 sites after protein-normalisation) — so the robust, abundance-independent signalling difference is tyrosine-level.

### 2.4 Complex-level and interaction view (run 0004)
Co-abundance reliably recovers known complexes at high confidence (45× enrichment at bicor ≥ 0.9). Whole co-regulated **modules do not define the subtypes**, but **complex abundance does**: OXPHOS/TCA complexes are coordinately high in Metabolism, gene-expression/turnover complexes (CDC73–PAF1, eIF3, SCF) high in EMT — the same Axis 1, now at complex resolution, with subunits moving together.

---

## 3. Overlapping, complementary and conflicting findings

- **Overlapping (same axis, multiple methods):** oxidative metabolism appears as a pathway (r2), a set of assembled complexes (r4), and the backdrop to the RTK/MEK signalling state (r3); the GI-lineage identity (r1) is the likely upstream driver. These are **not independent** — they are different measurements of one programme.
- **Complementary (add genuinely new layers):** r1 adds the lineage-TF *driver*; r3 adds the *signalling actionability* (RTK/MEK vs AKT-mTOR/SGK) and the key methodological correction (tyrosine, not bulk S/T); r2 adds a specific *metabolic vulnerability* (glutathione/ferroptosis) and the *orthogonal* MYC-proteostasis axis; r4 adds *complex-level* coordination.
- **Conflicting / corrected:** run 0001 proposed TMK1 might be a cryptic EBV carrier (it scored on the EBV proteome module); the run 0001 direct EBV peptide search **found no EBV peptides in TMK1**, so the carrier reading was rejected and the module membership re-attributed to shared cell state. Run 0002's "no Warburg trade-off" (H2.3) is weak (p = 0.07) and should not be over-read. No hard contradictions otherwise.

---

## 4. Critical evaluation of all hypotheses (`tables/hypothesis_register.csv`, `figures/02_hypothesis_prioritisation.png`)

All 12 prior hypotheses are registered with a documented rating of evidence strength, robustness, novelty, status and key limitation. Headlines:
- **Best-supported:** 1.H1 (lineage TFs), 3.H1 (RTK/MEK), 3.H2 (AKT-mTOR/SGK), 2.H1 (glutathione/ferroptosis), 1.H3 (EBV immune evasion), 4.H2 (complex abundance), plus the two methodological ones (3.H3 tyrosine-robust, 4.H3 co-abundance maps known complexes).
- **Weak / insufficient:** 2.H3 (no-Warburg, p = 0.07) and 4.H1 (novel ER/secretory co-regulation — co-abundance only, n = 14, nominal enrichment): candidate-level, explicitly not asserted.
- **Partly rejected:** 1.H2 (TMK1 EBV-carrier).
- **Recurring confounder across nearly all:** n = 14, and subtype/GI-lineage/EBV are mutually confounded, so "driven by metabolism" vs "driven by lineage/differentiation" cannot be separated from this dataset.

---

## 5. Prioritised hypotheses (top 5) (`tables/prioritised_hypotheses.csv`)

Ranked by integrated evidence and actionability; the "supporting analyses" count is cross-method agreement on **one** dataset, not independent validation.

**P1 — The "Metabolism" subtype is a well-differentiated, HNF-driven oxidative-metabolic state.** *Integrates 1.H1 + 2.H1 + 4.H2 + metabolic arm of 3.H1 (4 analyses).* Evidence: lineage TFs, OXPHOS/FAO pathways, glutathione detox and OXPHOS complexes all co-vary (r = 0.66) and split the subtypes; consistent with the authors' intestinal/G-Int mapping. Limitation: convergence is single-dataset; lineage-vs-metabolism causality unresolved. **Validation:** HNF4A/HNF1A knock-down / over-expression → OXPHOS proteome + oxygen-consumption (Seahorse); HNF4A ChIP for FAO/OXPHOS targets; replicate in an independent GC cohort.

**P2 — Metabolism subtype depends on RTK→RAS-ERK (MEK node); EMT subtype on PI3K-AKT-mTOR/SGK.** *Integrates 3.H1 + 3.H2.* Evidence: motif-based kinase activity independently reproduces the paper's MEK (Metabolism) vs mTOR (EMT) sensitivities and adds the upstream RTK layer and the SGK arm (explaining pan-AKT failure). Limitation: PTM data are from the paper supplement (not quantms); kinase inference is motif-based; the paper's drug data use the same lines. **Validation:** pERK vs pAKT/pNDRG1 by subtype; MEK vs mTOR vs SGK inhibitor panels; single- vs multi-RTK inhibition to test RTK redundancy.

**P3 — The Metabolism subtype carries a glutathione/antioxidant programme → ferroptosis (GPX4/GCLC) vulnerability.** *2.H1.* Evidence: GCLC/GCLM/GPX4/GSTs/G6PD coordinately higher (p = 0.001), co-regulated with OXPHOS; STRING confirms a glutathione/ROS network (FDR 1e-22). Specific, actionable, previously un-nominated. Limitation: single analysis; abundance ≠ flux; vulnerability untested. **Validation:** GPX4i (RSL3/ML210) + erastin + BSO sensitivity, Metabolism vs EMT; DepMap GPX4/GCLC/SLC7A11 dependency; lipid-ROS assay.

**P4 — A MYC-target / proteostasis (UPR) proliferation axis varies orthogonally to the metabolic subtype.** *2.H2.* Evidence: MYC–UPR r = 0.87, independent of subtype (r = −0.24) → a second heterogeneity axis. Limitation: MYC/UPR gene-set overlap inflates r; single dataset. **Validation:** re-test after removing set overlap; proteasome/HSP90/IRE1 inhibitor sensitivity in MYC-high lines.

**P5 — In-vitro EBV infection imposes a cell-intrinsic immune-evasion state (MHC-I / IFN down).** *1.H3.* Evidence: the cleanest contrast in the project (isogenic, replicate-level AGS→eAGS): EBV lowers MHC-I antigen presentation and the IFN response. Limitation: one infected line; eAGS also more proliferative. **Validation:** surface HLA-I ± IFN-γ flow cytometry (AGS vs eAGS and other freshly infected lines); CD8 killing; EBV latent-gene (BNLF2a) dependence.

---

## 6. Knowledge gaps, alternative explanations, confounders

- **No independent validation.** Four analyses, one dataset — the central limitation. The most valuable next step is *orthogonal* data, not more re-analysis.
- **Causality (lineage vs metabolism).** P1 cannot say whether HNF-lineage drives the metabolism or vice-versa, or whether both are co-markers of differentiation. HNF perturbation is the discriminating experiment.
- **n = 14 and confounding.** Subtype, GI lineage and EBV status are entangled (all 3 natural-EBV+ lines are "Metabolism"); correlation-based inference (co-abundance, kinase activity, pathway scores) is underpowered and abundance-based (abundance ≠ activity/flux/assembled-complex).
- **Data provenance.** The signalling layer (run 0003) is the authors' supplementary phospho data, not the quantms release; the quantms release is whole-proteome only. Kinase and complex inferences are motif-/correlation-based predictions.
- **Alternative explanations.** Metabolic/translation differences could reflect proliferation rate or culture-media adaptation rather than fixed subtype biology; the RTK "network" is a motif-family signal, not proof each receptor is active; "novel" co-abundance pairs may be indirect co-regulation.

---

## 7. Most promising future directions

1. **One orthogonal experiment settles the most.** Across P1–P3 the Metabolism subtype is predicted to be **oxidative, antioxidant-dependent and MEK-driven**: a single functional screen (Seahorse OCR, RSL3/erastin, trametinib, plus HNF4A dependency) in a Metabolism-vs-EMT line panel would test the integrated model at once.
2. **Cross-validate in independent cohorts** (CPTAC/other GC proteomes, DepMap/CCLE dependency & drug data) — turns single-dataset convergence into real validation and would confirm/refute the RTK, ferroptosis and SGK predictions (needs DepMap/CCLE network access, currently blocked).
3. **Resolve lineage vs metabolism causality** by HNF1A/HNF4A perturbation with proteomic + metabolic read-out.
4. **Scale the interaction map** — co-abundance over many more conditions, or co-fractionation MS, to move run 0004's candidate pairs from co-regulation to physical interaction.

---

## 8. Files
- `results/REPORT.md` (this synthesis); `results/FIGURES_AND_TABLES.md`.
- `figures/01_cross_run_scorecard.png` (per-line axis scorecard + axis convergence), `02_hypothesis_prioritisation.png`.
- `tables/per_line_axis_scores.csv`, `axis_subtype_association.csv`, `hypothesis_register.csv` (all 12 prior hypotheses evaluated), `prioritised_hypotheses.csv` (top 5), `integration_summary.json`, `register_summary.json`.
- Source material: the full reports, figures and tables of runs 0001–0004 under `runs/0001`–`runs/0004`.
