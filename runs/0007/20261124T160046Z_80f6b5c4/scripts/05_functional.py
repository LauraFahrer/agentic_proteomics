#!/usr/bin/env python3
"""Step 05 - functional annotation with self-implemented methods (Fig. 2B, 4B, 4D, 5B, S2C/D, S4C/D/G; Table S6).
 (a) GSEA (methods.gsea): S2N, weighted KS, 1000 permutations (phenotype and gene-set), NES, FDR; MSigDB v7.0;
     phosphosites collapsed to genes by per-sample maximum (GSEA 'Max_probe').
 (b) WGCNA-type modules (networks_kinases.wgcna_modules), module-cluster correlation, GO-BP ORA.
 (c) one-way ANOVA (scipy.stats.f_oneway) / moderated t (methods.moderated_t) + GO-BP ORA (methods.ora).
 (d) integrative EMT vs Metabolism moderated t (|log2FC| > 1, adj.P < 0.05) and Fig. 5B named-protein check."""
import json
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from common import *
import methods as M
import networks_kinases as NK

step = sys.argv[1] if len(sys.argv) > 1 else "all"
H = M.read_gmt(os.path.join(DB, "h.all.v7.0.symbols.gmt"))
KEGG = M.read_gmt(os.path.join(DB, "c2.cp.kegg.v7.0.symbols.gmt"))
BP = M.read_gmt(os.path.join(DB, "c5.bp.v7.0.symbols.gmt"))
A = pd.read_csv(os.path.join(TAB, "04_subtype_assignments.csv"))
INT = pd.read_csv(os.path.join(TAB, "03_integrative_clusters.csv")); intlab = dict(zip(INT.line, INT.integrative_subtype))
WP, wa = load_WP(); PP, pa = load_PP(); YP, ya = load_YP()
gene = {"WP": wa["Gene"], "PP": pa["GeneSymbol"], "YP": ya["GeneSymbol"]}
mats = {"WP": WP, "PP": PP, "YP": YP}


def labels(d):
    a = A[A.dataset == d].set_index("sample"); return a.loc[mats[d].columns, "cluster"].values


# map reproduced cluster -> paper cluster (majority) for naming
cmap = A.groupby(["dataset", "cluster"]).paper_cluster.agg(lambda x: x.value_counts().index[0]).to_dict()

if step in ("all", "gsea"):
    res = []
    for d, sets in [("WP", H), ("PP", KEGG), ("YP", KEGG)]:
        X = mats[d] if d == "WP" else M.collapse_max_probe(mats[d], gene[d].values)
        if d == "WP":
            X = X.copy(); X.index = gene["WP"].values; X = X.groupby(level=0).max()
        lab = labels(d)
        for c in sorted(set(lab)):
            for perm in ["phenotype", "gene_set"]:
                r = M.gsea(X.values, X.index.values, lab == c, sets, nperm=1000, permute=perm, seed=149)
                r["dataset"], r["cluster"], r["paper_cluster"], r["permutation"] = d, c, cmap[(d, c)], perm
                res.append(r)
                print(d, c, perm, "FDR<0.25 & NES>0:", int(((r.FDR_q < 0.25) & (r.NES > 0)).sum()), flush=True)
    G = pd.concat(res); G.to_csv(os.path.join(TAB, "05_gsea_all.csv"), index=False)

    # comparison with paper (Table S6, Fig. 4B, 4D)
    s6 = pd.read_excel(os.path.join(DC, "supp/x/MOL2-18-1437-s004.xlsx"), skiprows=2)
    ps = []
    for col, grp in [(0, 1), (2, 2), (4, 3)]:
        ps += [("WP", grp, s) for s in s6.iloc[:, col].dropna() if str(s).startswith("HALLMARK")]
    ps += [("PP", 2, "KEGG_" + s) for s in ["NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "CALCIUM_SIGNALING_PATHWAY",
           "FC_EPSILON_RI_SIGNALING_PATHWAY", "CHEMOKINE_SIGNALING_PATHWAY", "FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS", "GAP_JUNCTION",
           "VEGF_SIGNALING_PATHWAY", "LONG_TERM_POTENTIATION", "REGULATION_OF_ACTIN_CYTOSKELETON", "GNRH_SIGNALING_PATHWAY", "ENDOCYTOSIS"]]
    ps += [("YP", 3, "KEGG_" + s) for s in ["JAK_STAT_SIGNALING_PATHWAY", "CELL_ADHESION_MOLECULES_CAMS", "GNRH_SIGNALING_PATHWAY",
           "FOCAL_ADHESION", "FC_EPSILON_RI_SIGNALING_PATHWAY", "ERBB_SIGNALING_PATHWAY", "T_CELL_RECEPTOR_SIGNALING_PATHWAY",
           "INSULIN_SIGNALING_PATHWAY", "NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "CHEMOKINE_SIGNALING_PATHWAY"]]
    ps += [("YP", 1, "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION")]
    P = pd.DataFrame(ps, columns=["dataset", "paper_cluster", "set"]); P.to_csv(os.path.join(TAB, "05_paper_reported_sets.csv"), index=False)
    cmp = []
    for (d, pc), p in P.groupby(["dataset", "paper_cluster"]):
        for perm in ["phenotype", "gene_set"]:
            g = G[(G.dataset == d) & (G.paper_cluster == pc) & (G.permutation == perm)]
            sig = set(g[(g.FDR_q < 0.25) & (g.NES > 0)].set)
            rk = g.sort_values("NES", ascending=False).set.tolist()
            cmp.append(dict(dataset=d, paper_cluster=pc, permutation=perm, n_paper=len(p), n_this_run_FDR25=len(sig),
                            n_overlap=len(sig & set(p.set)), n_paper_with_positive_NES=int(g[g.set.isin(p.set)].NES.gt(0).sum()),
                            median_rank_of_paper_sets=np.median([rk.index(s) + 1 for s in p.set if s in rk]) if any(s in rk for s in p.set) else np.nan,
                            n_sets_tested=len(g)))
    cmp = pd.DataFrame(cmp); cmp.to_csv(os.path.join(TAB, "05_gsea_vs_paper_summary.csv"), index=False); print(cmp.to_string())
    fig2b = ["HALLMARK_" + s for s in ["COMPLEMENT", "GLYCOLYSIS", "UV_RESPONSE_UP", "ALLOGRAFT_REJECTION", "P53_PATHWAY",
             "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "KRAS_SIGNALING_UP", "PEROXISOME", "XENOBIOTIC_METABOLISM",
             "BILE_ACID_METABOLISM", "ADIPOGENESIS", "FATTY_ACID_METABOLISM", "OXIDATIVE_PHOSPHORYLATION", "APOPTOSIS",
             "PI3K_AKT_MTOR_SIGNALING", "MYC_TARGETS_V1", "UV_RESPONSE_DN", "MITOTIC_SPINDLE", "TGF_BETA_SIGNALING",
             "TNFA_SIGNALING_VIA_NFKB", "HEDGEHOG_SIGNALING", "E2F_TARGETS", "G2M_CHECKPOINT", "HYPOXIA",
             "EPITHELIAL_MESENCHYMAL_TRANSITION", "WNT_BETA_CATENIN_SIGNALING"]]
    for perm in ["phenotype", "gene_set"]:
        fig, axs = plt.subplots(1, 3, figsize=(18, 9), gridspec_kw=dict(width_ratios=[3, 1.3, 2]))
        for ax, (d, sets, names) in zip(axs, [("WP", fig2b, WP_NAMES), ("PP", P[P.dataset == "PP"].set.tolist(), {1: "PP1 EMT", 2: "PP2 Met"}),
                                              ("YP", P[P.dataset == "YP"].set.tolist(), {1: "YP1", 2: "YP2", 3: "YP3"})]):
            g = G[(G.dataset == d) & (G.permutation == perm)]
            cols = sorted(g.paper_cluster.unique())
            mat = np.zeros((len(sets), len(cols)))
            for j, c in enumerate(cols):
                gg = g[g.paper_cluster == c].set_index("set")
                for i, s in enumerate(sets):
                    if s in gg.index and gg.loc[s, "NES"] > 0:
                        mat[i, j] = -np.log10(max(gg.loc[s, "FDR_q"], 1e-3))
            ax.imshow(mat, cmap="Reds", vmin=0, vmax=3, aspect="auto")
            for i in range(mat.shape[0]):
                for j in range(mat.shape[1]):
                    ax.text(j, i, f"{mat[i,j]:.1f}", ha="center", va="center", fontsize=6)
            ax.set_yticks(range(len(sets))); ax.set_yticklabels([s.split("_", 1)[1] for s in sets], fontsize=7)
            ax.set_xticks(range(len(cols))); ax.set_xticklabels([names[c] for c in cols], rotation=45, fontsize=8)
            ax.set_title(f"{d}: -log10 FDR ({perm} perm.)", fontsize=9)
        savefig(fig, f"05_Fig2B_4BD_GSEA_{perm}_perm")

if step in ("all", "wgcna"):
    summ = {}
    for d in ["WP", "PP"]:
        X = mats[d].values.T
        sft = NK.scale_free_fit(X); sft.to_csv(os.path.join(TAB, f"05_WGCNA_{d}_soft_threshold.csv"), index=False)
        ok = sft[sft.signed_R2 >= 0.8]; beta = int(ok.power.iloc[0]) if len(ok) else 6
        lab = NK.wgcna_modules(X, beta, min_size=30, merge_cut=0.25, gap=0.1)
        ME = NK.eigengenes(X, lab)
        cl = labels(d)
        rows = []
        uni = set(gene[d])
        for m in ME.columns:
            r = dict(module=int(m), size=int((lab == m).sum()))
            for c in sorted(set(cl)):
                rr, pv = stats.pearsonr(ME[m], (cl == c).astype(float)); r[f"r_C{c}"] = rr; r[f"p_C{c}"] = pv
            o = M.ora(set(gene[d][lab == m]), uni, BP)
            r["top_GO_BP"] = "; ".join(o[o.adjP < 0.05].set.head(3))
            rows.append(r)
        S = pd.DataFrame(rows); S.to_csv(os.path.join(TAB, f"05_WGCNA_{d}_modules.csv"), index=False)
        summ[d] = dict(power=beta, n_modules=int(len(ME.columns)), n_unassigned=int((lab == 0).sum()),
                       paper_n_modules={"WP": 12, "PP": 5}[d])
        print(d, summ[d], flush=True)
        rc = [c for c in S.columns if c.startswith("r_C")]
        fig, ax = plt.subplots(figsize=(5, max(4, 0.25 * len(S))))
        ax.imshow(S[rc].values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_yticks(range(len(S))); ax.set_yticklabels([f"M{m} ({s})" for m, s in zip(S.module, S["size"])], fontsize=6)
        ax.set_xticks(range(len(rc))); ax.set_xticklabels(rc); ax.set_title(f"{d} modules vs NMF clusters (beta={beta})")
        savefig(fig, f"05_WGCNA_{d}_module_cluster_correlation")
    json.dump(summ, open(os.path.join(TAB, "05_WGCNA_summary.json"), "w"), indent=1)

if step in ("all", "diff"):
    out = {}
    def go_by_cluster(d, sig_idx, lab):
        X = mats[d].values; uni = set(gene[d]); res = []
        means = np.array([X[sig_idx][:, lab == c].mean(axis=1) for c in sorted(set(lab))]).T
        top = np.array(sorted(set(lab)))[means.argmax(axis=1)]
        for c in sorted(set(lab)):
            o = M.ora(set(gene[d].values[sig_idx][top == c]), uni, BP); o = o[o.adjP < 0.05].copy()
            o["cluster"], o["paper_cluster"], o["n_features_up"] = c, cmap[(d, c)], int((top == c).sum()); res.append(o)
        return pd.concat(res)
    for d in ["WP", "YP"]:
        lab = labels(d); p = M.oneway_anova(mats[d].values, lab); q = M.bh(p)
        sig = np.where(q < 0.05)[0]; out[f"{d}_ANOVA_sig"] = len(sig)
        go_by_cluster(d, sig, lab).to_csv(os.path.join(TAB, f"05_{d}_ANOVA_GO_BP.csv"), index=False)
    lab = labels("PP"); X = np.column_stack([np.ones(len(lab)), (lab == 1).astype(float)])
    r, _ = M.moderated_t(PP.values, X, [0, 1]); sig = np.where((r.adjP < 0.05) & (r.logFC.abs() > 1))[0]
    out["PP_moderated_t_sig"] = len(sig)
    go_by_cluster("PP", sig, lab).to_csv(os.path.join(TAB, "05_PP_modt_GO_BP.csv"), index=False)
    # integrative EMT vs Metabolism
    allr = []
    for d in ["WP", "PP", "YP"]:
        g = np.array([intlab[line_of(s)] == "EMT" for s in mats[d].columns], float)
        r, prior = M.moderated_t(mats[d].values, np.column_stack([np.ones(len(g)), g]), [0, 1])
        r.insert(0, "id", mats[d].index); r.insert(1, "gene", gene[d].values); r["dataset"] = d
        allr.append(r.rename(columns={"logFC": "log2FC_EMT_vs_Met"}))
    I = pd.concat(allr); I["sig"] = (I.adjP < 0.05) & (I.log2FC_EMT_vs_Met.abs() > 1)
    I.to_csv(os.path.join(TAB, "05_integrative_modt_all.csv.gz"), index=False)
    named = [("AAMP", "EMT", "protein"), ("CCN1", "EMT", "protein"), ("LGALS1", "EMT", "protein"), ("TGFBI", "EMT", "protein"),
             ("STK3", "EMT", "protein"), ("STK4", "EMT", "protein"), ("TGFB1I1", "EMT", "phospho"), ("CTNNB1", "EMT", "phospho"),
             ("YAP1", "EMT", "phospho"), ("ABL2", "EMT", "phospho"), ("APC", "EMT", "phospho"), ("PHLDB2", "EMT", "both"),
             ("GLUL", "Metabolism", "protein"), ("ALAD", "Metabolism", "protein"), ("ACAA1", "Metabolism", "protein"),
             ("ACAA2", "Metabolism", "protein"), ("GCLC", "Metabolism", "protein"), ("AK3", "Metabolism", "protein"),
             ("G6PD", "Metabolism", "phospho"), ("GSTA2", "Metabolism", "phospho"), ("ALDH3A1", "Metabolism", "both"),
             ("ACSL5", "Metabolism", "both"), ("AKR1A1", "Metabolism", "both")]
    chk = []
    for gname, grp, lev in named:
        sgn = 1 if grp == "EMT" else -1
        w = I[(I.dataset == "WP") & (I.gene == gname)]; p = I[(I.dataset == "PP") & (I.gene == gname)]
        chk.append(dict(gene=gname, paper_group=grp, paper_level=lev,
                        WP_log2FC=round(w.log2FC_EMT_vs_Met.iloc[0], 2) if len(w) else np.nan,
                        WP_sig_in_paper_direction=bool(len(w) and w.sig.iloc[0] and np.sign(w.log2FC_EMT_vs_Met.iloc[0]) == sgn),
                        PP_n_sites=len(p), PP_sites_sig_in_paper_direction=int((p.sig & (np.sign(p.log2FC_EMT_vs_Met) == sgn)).sum())))
    chk = pd.DataFrame(chk); chk.to_csv(os.path.join(TAB, "05_Fig5B_named_proteins_check.csv"), index=False); print(chk.to_string())
    oras = []
    for grp, sgn in [("EMT", 1), ("Metabolism", -1)]:
        s = I[I.sig & I.dataset.isin(["WP", "PP"]) & (np.sign(I.log2FC_EMT_vs_Met) == sgn)]
        uni = set(I[I.dataset.isin(["WP", "PP"])].gene)
        for db, sets in [("KEGG", KEGG), ("H", H)]:
            o = M.ora(set(s.gene), uni, sets); o["group"], o["db"] = grp, db; oras.append(o[o.adjP < 0.25])
    pd.concat(oras).to_csv(os.path.join(TAB, "05_Fig5B_integrative_pathway_ORA.csv"), index=False)
    json.dump(out, open(os.path.join(TAB, "05_differential_counts.json"), "w"), indent=1)
