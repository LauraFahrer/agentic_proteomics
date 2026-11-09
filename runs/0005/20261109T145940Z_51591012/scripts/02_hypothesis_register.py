"""Register and critical evaluation of every hypothesis from runs 0001-0004, and the
prioritised top-5 integrated hypotheses. The scores are the synthesis author's critical
judgement (documented here so they are auditable), not a statistical computation.

Scoring (1-5): evidence = strength/number of supporting observations; convergence = how many
independent ANALYSES (not datasets) support it; novelty = beyond the paper's own conclusions;
stat/robustness = statistical strength & freedom from the key confounder. All four runs share
ONE dataset, so convergence is internal consistency, NOT independent experimental validation.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CAT, INK2, SUBTYPE_COLORS, TAB, savefig, style

style()

# --- register of all prior hypotheses -------------------------------------
REG = [
    # run, id, short, supporting analyses, evidence(1-5), robustness(1-5), novelty(1-5), status, key_limitation
    ("0001", "1.H1", "GI-lineage TFs (HNF1A/HNF4A/FOXA2) define the Metabolism subtype",
     "presence 7/7 vs 0/7; rho=0.74 with metabolic score; confirmed in authors' Spectronaut matrix", 5, 4, 3, "supported",
     "lineage/subtype/metabolism confounded; TF level != activity"),
    ("0001", "1.H2", "Natural-EBV proteome module; TMK1 a cryptic EBV carrier",
     "111/53-protein module in 3 EBV+ lines; TMK1 scores on it", 3, 2, 3, "partly rejected",
     "TMK1-carrier claim REJECTED by direct EBV peptide search; module confounded with lineage; n=3 EBV+"),
    ("0001", "1.H3", "In-vitro EBV infection (eAGS) lowers MHC-I antigen presentation & IFN response",
     "isogenic AGS->eAGS: HLA/TAP/PSME/IFN down (antigen-presentation q=4e-5)", 4, 4, 3, "supported",
     "single infected line; eAGS also more proliferative (could lower MHC-I indirectly)"),
    ("0002", "2.H1", "Metabolism subtype runs a glutathione/antioxidant detox programme -> GPX4/GCLC (ferroptosis) vulnerability",
     "detox proteins Metab>EMT p=0.001; detox-OXPHOS r=0.77; STRING glutathione FDR 1e-22", 4, 4, 4, "supported",
     "abundance != flux/GSH level; prediction of vulnerability untested; n=14"),
    ("0002", "2.H2", "MYC-target load is coupled to the proteostasis/UPR machinery (proliferation axis)",
     "MYC vs UPR r=0.87 p=5e-5; vs DNA-repair r=0.67; orthogonal to subtype", 4, 3, 3, "supported",
     "MYC & UPR gene sets share members (partly definitional); MYC level != activity"),
    ("0002", "2.H3", "No Warburg glycolysis<->OXPHOS trade-off; an OXPHOS-high subset exists",
     "glycolysis-OXPHOS r=+0.50 p=0.07; NCC24/YCCEL1 OXPHOS-high", 2, 2, 2, "weak/insufficient",
     "borderline (p=0.07); abundance != respiration; n=14"),
    ("0003", "3.H1", "Metabolism subtype = convergent RTK (MET/FGFR/RON/Eph) -> RAS-ERK -> MEK dependency",
     "24/25 differential pY on RTK-adaptor module; phosX RTK+MAP2K activity; matches paper's MEK sensitivity", 5, 4, 4, "supported",
     "phosX can't resolve same-motif kinases; PTM data from paper supplement not quantms; n=14"),
    ("0003", "3.H2", "EMT subtype = PI3K-AKT-mTOR with a co-active SGK arm -> mTOR/SGK dependency",
     "phosX AKT1/2/3+SGK1/3+P70S6K activity; matches paper's mTOR sensitivity & pan-AKT failure", 4, 4, 4, "supported",
     "S/T kinase inference abundance-confounded; motif overlap AKT/SGK; n=14"),
    ("0003", "3.H3", "The robust subtype signalling difference is tyrosine-level; bulk S/T phospho is abundance-coupled",
     "0/8452 S/T sites differ after protein-normalisation; 25 pY sites do", 4, 4, 3, "supported (methodological)",
     "protein-normalisation loses power; detectability not non-existence"),
    ("0004", "4.H1", "Candidate novel ER/secretory-membrane co-regulation (ERP29-TMED2, MTCH2-VAPA)",
     "strong co-abundance (bicor>=0.95) with 0 STRING/GO evidence; nominal ER-membrane enrichment", 2, 2, 4, "weak/candidate",
     "co-abundance != interaction; enrichment not FDR-significant; n=14 underpowered for novel PPI"),
    ("0004", "4.H2", "Subtype-coordinated complex abundance: OXPHOS complexes high in Metabolism, gene-expression complexes in EMT",
     "8 complexes q<0.1; OXPHOS/TCA complexes Metab-high move together (cohesion 0.46-0.52)", 4, 3, 2, "supported",
     "mean-of-subunits != assembled complex; same axis as 2.H1/4 (not independent)"),
    ("0004", "4.H3", "Co-abundance here is a reliable known-complex map but underpowered for de-novo PPI discovery",
     "strong edges 45x enriched for complexes; AUC 0.61; 62% of strong pairs already known", 4, 4, 2, "supported (methodological)",
     "n=14; absence of STRING evidence != absence of interaction"),
]
reg = pd.DataFrame(REG, columns=["run", "id", "hypothesis", "supporting_analyses", "evidence",
                                 "robustness", "novelty", "status", "key_limitation"])
reg.to_csv(TAB / "hypothesis_register.csv", index=False)

# --- prioritised top-5 integrated hypotheses ------------------------------
PRIOR = [
    ("P1", "The 'Metabolism' subtype is a well-differentiated, HNF-driven OXIDATIVE-METABOLIC state",
     "1.H1 + 2.H1 + 4.H2 (+ metabolic arm of 3.H1)", 4,
     "Lineage TFs, OXPHOS/FAO pathways, glutathione detox, and OXPHOS complexes all co-vary (mean r=0.66) and all split EMT/Metab; recovered by 4 analyses and consistent with the authors' G-Int/intestinal mapping.",
     "Convergence is across methods on ONE dataset, not independent validation; lineage vs metabolism causally unresolved.",
     "HNF4A/HNF1A knockdown or over-expression -> measure OXPHOS proteome + oxygen consumption (Seahorse); HNF4A ChIP for FAO/OXPHOS targets; test in an independent GC proteomic/transcriptomic cohort."),
    ("P2", "Metabolism subtype depends on RTK->RAS-ERK signalling (MEK node); EMT subtype on PI3K-AKT-mTOR/SGK",
     "3.H1 + 3.H2 (cross-validated by the paper's own drug data)", 2,
     "Motif-based kinase activity (phosX) independently reproduces the paper's MEK (Metabolism) vs mTOR (EMT) sensitivities and adds the upstream RTK layer (Metabolism) and the SGK arm (EMT, explaining pan-AKT failure).",
     "PTM data are from the paper supplement, not quantms; kinase inference is motif-based; the paper's drug data are not independent of the same cell lines.",
     "pERK vs pAKT/pNDRG1 immunoblot by subtype; MEK vs mTOR vs SGK inhibitor sensitivity; multi-RTK vs single-RTK inhibition (test RTK redundancy)."),
    ("P3", "The Metabolism subtype carries a glutathione/antioxidant programme -> ferroptosis (GPX4/GCLC) vulnerability",
     "2.H1 (STRING-validated gene network)", 1,
     "GCLC/GCLM/GPX4/GSTs/G6PD coordinately higher in Metabolism (p=0.001), co-regulated with OXPHOS; STRING confirms a glutathione/ROS network (FDR 1e-22). A specific, actionable, previously un-nominated vulnerability.",
     "Single analysis; abundance not flux; the vulnerability itself is a prediction, untested.",
     "GPX4 inhibitor (RSL3/ML210) + erastin + BSO sensitivity in Metabolism vs EMT lines; DepMap GPX4/GCLC/SLC7A11 dependency; lipid-ROS measurement."),
    ("P4", "A MYC-target / proteostasis (UPR) proliferation axis varies ORTHOGONALLY to the metabolic subtype",
     "2.H2 (+ run0003 proliferation signal)", 1,
     "MYC-target abundance tracks the UPR machinery (r=0.87) and DNA repair; this axis is statistically independent of the metabolic subtype axis (r=-0.24), so it is a second source of heterogeneity.",
     "MYC/UPR gene-set member overlap inflates r; single-dataset; proliferation state also co-varies with culture conditions.",
     "After removing set overlap, re-test MYC-UPR coupling; proteasome/HSP90/IRE1 inhibitor sensitivity in MYC-high lines (NCC24, KatoIII, eAGS, AGS)."),
    ("P5", "In-vitro EBV infection imposes a cell-intrinsic immune-evasion state (MHC-I / IFN down)",
     "1.H3 (isogenic AGS->eAGS)", 1,
     "The cleanest contrast in the project (isogenic, replicate-level): EBV infection lowers MHC-I antigen presentation and the interferon response. Mechanistically coherent with herpesvirus immune evasion.",
     "One infected line; eAGS is also more proliferative; not linked to the main subtype axis.",
     "Flow cytometry for surface HLA-I +/- IFN-gamma in AGS vs eAGS and other freshly infected lines; CD8 T-cell killing; test EBV latent-gene (BNLF2a) dependence."),
]
prio = pd.DataFrame(PRIOR, columns=["priority", "hypothesis", "integrates", "n_supporting_analyses",
                                    "evidence_summary", "key_limitation", "validation"])
prio.to_csv(TAB / "prioritised_hypotheses.csv", index=False)
json.dump({"n_hypotheses_registered": len(reg), "n_prioritised": len(prio),
           "status_counts": reg.status.value_counts().to_dict()},
          open(TAB / "register_summary.json", "w"), indent=2)
print(reg[["id", "status", "evidence", "robustness", "novelty"]].to_string())
print("\nPrioritised:\n", prio[["priority", "hypothesis", "n_supporting_analyses"]].to_string())

# --- figure: register overview + priorities -------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"wspace": 0.05, "width_ratios": [1.15, 1]})
ax = axes[0]
reg2 = reg.copy()
y = np.arange(len(reg2))[::-1]
runcol = {"0001": CAT[0], "0002": CAT[2], "0003": CAT[3], "0004": CAT[4]}
ax.barh(y, reg2.evidence, color=[runcol[r] for r in reg2.run], height=0.6)
for yi, (_, r) in zip(y, reg2.iterrows()):
    ax.text(0.05, yi, f"{r.id}  {r.hypothesis[:52]}", va="center", fontsize=6.6, color="white" if r.evidence > 1.5 else INK2)
    mark = {"supported": "", "partly rejected": "  [carrier rejected]", "weak/insufficient": "  [weak]",
            "weak/candidate": "  [candidate]"}.get(r.status, "")
    if mark:
        ax.text(r.evidence + 0.05, yi, mark, va="center", fontsize=6, color=INK2)
ax.set_yticks([]); ax.set_xlabel("evidence strength (author rating 1-5)")
ax.set_xlim(0, 6.5)
ax.set_title("A  All 12 hypotheses from runs 0001-0004\n(bar colour = source run)", fontsize=9.5)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=runcol[k], label=f"run {k[-1]}") for k in runcol], fontsize=7, loc="lower right")
ax = axes[1]
y = np.arange(len(prio))[::-1]
ax.barh(y, prio.n_supporting_analyses, color=CAT[0], height=0.55)
for yi, (_, r) in zip(y, prio.iterrows()):
    ax.text(0.05, yi + 0.0, f"{r.priority}: {r.hypothesis[:46]}", va="center", fontsize=6.8, color="white" if r.n_supporting_analyses >= 2 else INK2)
ax.set_yticks([]); ax.set_xlabel("# supporting analyses (of 4 runs)")
ax.set_xlim(0, 4.3)
ax.set_title("B  Prioritised integrated hypotheses\n(convergence = same data, not independent validation)", fontsize=9.5)
savefig(fig, "02_hypothesis_prioritisation")
plt.close(fig)
