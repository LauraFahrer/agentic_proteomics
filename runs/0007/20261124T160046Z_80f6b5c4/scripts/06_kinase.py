#!/usr/bin/env python3
"""Step 06 - kinase analyses with self-implemented methods.
  kinome   : kinome counts / family shares (Fig. 6A/B), outlier kinases (Fig. S6A-C, S6G), kinome NMF (Fig. S8B/C),
             40 differentially expressed kinases (Fig. 8B)
  activity : PP kinase activity, 6 matrices {raw, protein-normalised} x {mean, KSEA, MLR}; EMT vs Metabolism
             assignment counting (P < 0.05, |log2FC| > 1) (Fig. 8A, S8A)
  inka     : INKA scores on YP (Fig. 7A, 7D, S7A)
  ptmsea   : PTM-SEA (ssGSEA on PTMsigDB v1.9.1 flanking signatures) on PP (Fig. S5D)
Databases (external, used as data): HGNC/Manning kinase list (InKA resource, 2016), UniProt families (2026),
PhosMap human kinase-substrate table, PhosphoSitePlus+NetworKIN Jul-2016 (KSEA app release), InKA activation-loop
table, PTMsigDB v1.9.1, UniProt Swiss-Prot sequences."""
import json
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from common import *
import methods as M
import networks_kinases as NK

step = sys.argv[1] if len(sys.argv) > 1 else "all"
INT = pd.read_csv(os.path.join(TAB, "03_integrative_clusters.csv")); intlab = dict(zip(INT.line, INT.integrative_subtype))
WP, wa = load_WP(); PP, pa = load_PP(); YP, ya = load_YP()
KIN = set(pd.read_csv(os.path.join(DB, "Kinases_HGNC2017.01.20_v2.csv"))["Gene.Symbol.HGNC21Apr2016"])
fam = pd.read_csv(os.path.join(DC, "uniprot", "human_kinase_families.tsv"), sep="\t", keep_default_na=False)
fam.columns = ["acc", "gene", "family"]; fam = fam.drop_duplicates("gene").set_index("gene")


def group_of(g):
    f = fam.family.get(g, "")
    for k in ["AGC", "CAMK", "CK1", "CMGC", "STE", "TKL"]:
        if re.search(rf"\b{k} ", f):
            return k
    if "Tyr protein kinase family" in f:
        return "TK"
    if re.search("PI3/PI4-kinase|RIO-type|ABC1|PDK/BCKDK|FAST|TIF1|BRD", f):
        return "Atypical"
    return "Other"


TKMAP = {"Fibroblast growth factor receptor": "FGFR", "Ephrin receptor": "EPH", "SRC": "SRC", "EGF receptor": "EGFR",
         "JAK": "JAK", "Insulin receptor": "INSR", "ABL": "ABL", "AXL/UFO": "AXL", "CSF-1/PDGF receptor": "PDGFR",
         "FAK": "FAK", "ACK": "ACK", "Tec": "TEC", "SYK/ZAP-70": "SYK", "CSK": "CSK", "FES/FPS": "FER", "ROR": "ROR", "Trk": "TRK"}


def tk_family(g):
    if g in ("MET", "MST1R"): return "MET"
    if g in ("DDR1", "DDR2"): return "DDR"
    m = re.search(r"Tyr protein kinase family[.,] ?([^.;]*) subfamily", fam.family.get(g, ""))
    return TKMAP.get(m.group(1), m.group(1)) if m else "Other"


is_emt = lambda cols: np.array([intlab[line_of(c)] == "EMT" for c in cols])

# ====================================================================== kinome
if step in ("all", "kinome"):
    gW, gP, gY = wa["Gene"].values, pa["GeneSymbol"].values, ya["GeneSymbol"].values
    kW, kP, kY = np.isin(gW, list(KIN)), np.isin(gP, list(KIN)), np.isin(gY, list(KIN))
    cnt = pd.DataFrame(dict(dataset=["WP", "PP", "YP"], kinases=[len(set(gW[kW])), len(set(gP[kP])), len(set(gY[kY]))],
                            features=[kW.sum(), kP.sum(), kY.sum()], paper=[198, 193, 140]))
    cnt.to_csv(os.path.join(TAB, "06_kinome_counts.csv"), index=False); print(cnt)

    def share(m, genes, keyf):
        s = pd.Series((2 ** m).sum(axis=1).values, index=[keyf(g) for g in genes]).groupby(level=0).sum()
        return (100 * s / s.sum()).round(2)
    kYt = kY & np.array([group_of(g) == "TK" for g in gY])
    sh = {"WP": share(WP[kW], gW[kW], group_of), "PP": share(PP[kP], gP[kP], group_of), "YP": share(YP[kYt], gY[kYt], tk_family)}
    paper = {"WP": dict(Other=33.02, Atypical=22.01, STE=4.36, TK=4.03, AGC=0.48, CAMK=0.04, CK1=1.4e-05, TKL=0.07, CMGC=36.05),
             "PP": dict(CMGC=20, Other=8, Atypical=23, STE=11, TK=6, AGC=6, CAMK=20, TKL=2, CK1=3),
             "YP": dict(FGFR=32.97, EPH=32.89, SRC=14.67, EGFR=8.79, JAK=3.83, MET=3.37, INSR=1.18, DDR=0.77, ABL=0.67, AXL=0.3, Other=0.54)}
    FS = pd.concat([pd.DataFrame(dict(dataset=d, family=sorted(set(sh[d].index) | set(paper[d])),
                                      this_run_pct=[sh[d].get(f, np.nan) for f in sorted(set(sh[d].index) | set(paper[d]))],
                                      paper_pct=[paper[d].get(f, np.nan) for f in sorted(set(sh[d].index) | set(paper[d]))])) for d in sh])
    FS.to_csv(os.path.join(TAB, "06_Fig6AB_kinase_family_share.csv"), index=False); print(FS[FS.dataset == "YP"].to_string())
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.5))
    for ax, d in zip(axs, ["WP", "PP", "YP"]):
        x = FS[(FS.dataset == d) & (FS.this_run_pct > 0)].dropna(subset=["this_run_pct"])
        ax.pie(x.this_run_pct, labels=[f"{f} {v:.1f}%" for f, v in zip(x.family, x.this_run_pct)], textprops=dict(fontsize=7))
        ax.set_title(f"{d} kinome intensity share (this run)")
    savefig(fig, "06_Fig6AB_kinase_family_share")

    # outlier kinases
    def outl(m, genes, d):
        V = m.values; Zs = (V - V.mean(1, keepdims=True)) / V.std(1, ddof=1, keepdims=True)
        res = []
        for nm, Q in [("A_top5pct_level", V), ("B_top5pct_zscore", Zs)]:
            thr = np.nanquantile(Q, 0.95); ii, jj = np.where(Q >= thr)
            t = pd.DataFrame(dict(feature=m.index[ii], gene=genes[ii], line=[line_of(c) for c in m.columns[jj]]))
            t = t.groupby(["feature", "gene"]).line.agg(lambda x: ";".join(dict.fromkeys(x))).reset_index()
            t["n_lines"] = t.line.str.count(";") + 1; t["def"] = nm; t["dataset"] = d; res.append(t)
        return pd.concat(res)
    O = pd.concat([outl(WP[kW], gW[kW], "WP"), outl(PP[kP], gP[kP], "PP"), outl(YP[kY], gY[kY], "YP")])
    O.to_csv(os.path.join(TAB, "06_FigS6_outlier_kinases.csv"), index=False)
    print(O.groupby("def").gene.nunique())
    named = [("WP", r"^ERBB2$", "N87"), ("WP", r"^CDK12$", "N87"), ("PP", r"^ERBB2_", "N87"), ("PP", r"^CDK12_", "N87"),
             ("PP", r"^FGFR2_", "KatoIII"), ("YP", r"^FGFR2_", "KatoIII"), ("PP", r"^CSNK1D_S384", "MKN28"), ("PP", r"^CSNK1E_S389", "MKN28")]
    rec = []
    for d, pat, ln in named:
        for df_, g in O[(O.dataset == d) & O.feature.str.contains(pat)].groupby("def"):
            rec.append(dict(dataset=d, pattern=pat, paper_line=ln, definition=df_, n_features=len(g),
                            n_outlier_in_paper_line=int(g.line.str.contains(ln).sum()),
                            lines_flagged=";".join(dict.fromkeys(";".join(g.line).split(";")))))
    pd.DataFrame(rec).to_csv(os.path.join(TAB, "06_FigS6_named_outliers_recovery.csv"), index=False); print(pd.DataFrame(rec).to_string())
    s6g = []
    for g, pp in [("CDK4", 1e-4), ("CDK6", 0.0214)]:
        a = WP.loc[g, [c for c in WP.columns if line_of(c) == "YCCEL1"]]; b = WP.loc[g, [c for c in WP.columns if line_of(c) == "MKN28"]]
        s6g.append(dict(gene=g, YCCEL1=a.mean(), MKN28=b.mean(), student_t_P=stats.ttest_ind(a, b).pvalue, paper_P=pp))
    pd.DataFrame(s6g).to_csv(os.path.join(TAB, "06_FigS6G_CDK46.csv"), index=False); print(pd.DataFrame(s6g))
    lmeans = lambda m: pd.DataFrame({l: m[[c for c in m.columns if line_of(c) == l]].mean(axis=1) for l in LINES})
    feats = pd.concat([lmeans(WP.loc[["ERBB2", "CDK12", "CDK4", "CDK6", "SRC"]]),
                       lmeans(PP[PP.index.str.match(r"^(ERBB2|CDK12|FGFR2|CSNK1D|CSNK1E)_")]),
                       lmeans(YP[YP.index.str.match(r"^(ERBB2|FGFR2|SRC|EGFR|MET)_")])])
    Zl = feats.sub(feats.mean(1), axis=0).div(feats.std(1), axis=0)
    fig, ax = plt.subplots(figsize=(8, 0.18 * len(Zl) + 2))
    ax.imshow(Zl.values, cmap="RdBu_r", vmin=-3.5, vmax=3.5, aspect="auto")
    ax.set_yticks(range(len(Zl))); ax.set_yticklabels(Zl.index, fontsize=5); ax.set_xticks(range(14)); ax.set_xticklabels(LINES, rotation=90)
    ax.set_title("Named outlier kinases: line-mean z-scores (WP / PP / YP)")
    savefig(fig, "06_FigS6_named_outlier_features")

    # kinome NMF (Fig. S8B/C)
    km = WP[kW]
    C, (W_, H_, _) = M.nmf_consensus(km.values, 2, nrun=100, seed=123456)
    Z, cl = M.consensus_hclust(C, 2); D = 1 - C; np.fill_diagonal(D, 0); sil = M.silhouette(D, cl)
    integ = [intlab[line_of(c)] for c in km.columns]
    kres = pd.DataFrame(dict(sample=km.columns, kinome_cluster=cl, silhouette=sil, integrative=integ))
    kres.to_csv(os.path.join(TAB, "06_FigS8BC_kinome_NMF.csv"), index=False)
    ks = dict(n_kinases=len(km), cophenetic=M.cophenetic_coef(C), mean_silhouette=float(sil.mean()), ARI_vs_integrative=ari(cl, integ),
              misassigned=list(kres.sample[kres.kinome_cluster != kres.groupby("integrative").kinome_cluster.transform(lambda x: x.mode()[0])]))
    json.dump(ks, open(os.path.join(TAB, "06_FigS8BC_kinome_NMF_summary.json"), "w"), indent=1, default=float); print(ks)
    # 40 kinases (Fig. 8B)
    e = is_emt(km.columns).astype(float)
    r, _ = M.moderated_t(km.values, np.column_stack([np.ones(len(e)), e]), [0, 1]); r.index = km.index
    r = r.sort_values("P")
    p40 = ["BAZ1B", "STK24", "CASK", "KALRN", "STK39", "COQ8A", "RPS6KA3", "TAOK3", "BRD2", "RPS6KA1", "MST1R", "MET", "ABR", "MAP2K3",
           "MAPK13", "SRC", "EPHB2", "GAK", "PAK1", "CAMK2D", "LYN", "CDK2", "CDK1", "CSNK2A1", "CDK7", "TGM2", "STK3", "BMP2K", "STK4",
           "NRBP1", "SCYL2", "MAP3K7", "CDK9", "YES1", "MAP3K20", "ILK", "BRD4", "PTK7", "CAMK1", "IGF1R"]
    r["in_paper40"] = r.index.isin(p40); r.to_csv(os.path.join(TAB, "06_Fig8B_kinase_modt.csv"))
    rank = {k: i + 1 for i, k in enumerate(r.index)}
    chk = pd.DataFrame(dict(kinase=p40, in_kinase_list=[k in km.index for k in p40], rank_this_run=[rank.get(k) for k in p40]))
    chk.to_csv(os.path.join(TAB, "06_Fig8B_paper40_check.csv"), index=False)
    x40 = km.loc[[k for k in p40 if k in km.index]]
    from scipy.cluster import hierarchy
    hc = hierarchy.fcluster(hierarchy.linkage((x40.sub(x40.mean(1), axis=0)).T.values, "complete"), 2, "maxclust")
    s8 = dict(top40_overlap=int(r.index[:40].isin(p40).sum()), testable=int(chk.in_kinase_list.sum()),
              all_testable_in_top40=bool((chk.dropna().rank_this_run <= 40).all()), ARI_hclust40_vs_integrative=ari(hc, integ))
    json.dump(s8, open(os.path.join(TAB, "06_Fig8B_summary.json"), "w"), indent=1); print(s8)
    o = np.argsort([intlab[line_of(c)] + c for c in x40.columns])
    zz = x40.sub(x40.mean(1), axis=0).div(x40.std(1), axis=0)
    fig, ax = plt.subplots(figsize=(10, 7)); ax.imshow(zz.values[:, o], cmap="cividis", vmin=-2, vmax=2, aspect="auto")
    ax.set_yticks(range(len(zz))); ax.set_yticklabels(zz.index, fontsize=6); ax.set_xticks(range(len(o)))
    ax.set_xticklabels([f"{x40.columns[i]} ({intlab[line_of(x40.columns[i])][:3]})" for i in o], rotation=90, fontsize=5)
    ax.set_title("Paper's 40 kinases (WP), z-score"); savefig(fig, "06_Fig8B_kinase_heatmap")

# ====================================================================== PP kinase activity (6 matrices)
if step in ("all", "activity"):
    ks = pd.read_csv(os.path.join(DB, "human_ksrr.csv"))
    P = PP.copy(); P.index = pa["GeneSymbol"] + "_" + pa["Residue"]
    wpg = pd.Series(range(len(WP)), index=wa["Gene"].values); wpg = wpg[~wpg.index.duplicated()]
    has = P.index.str.split("_").str[0].isin(wpg.index)
    Pn = P[has] - WP.iloc[wpg[P.index[has].str.split("_").str[0]].values].values
    mats = {"raw": P.sub(P.mean(1), axis=0), "norm": Pn.sub(Pn.mean(1), axis=0)}
    print("sites raw", len(P), "normalised", len(Pn))
    act = {}
    for nm, m in mats.items():
        reg = NK.build_regulons(list(m.index), ks)
        act[f"{nm}_mean"] = NK.activity_mean(m.values, reg).set_axis(m.columns, axis=1)
        act[f"{nm}_ksea"] = NK.activity_ksea(m.values, reg).set_axis(m.columns, axis=1)
        act[f"{nm}_mlr"] = NK.activity_mlr(m.values, reg).set_axis(m.columns, axis=1)
        print(nm, "kinases with substrates:", len(reg), flush=True)
    os.makedirs(os.path.join(TAB, "kinase_activity"), exist_ok=True)
    rows = []
    for n, a in act.items():
        a.to_csv(os.path.join(TAB, "kinase_activity", f"06_PP_activity_{n}.csv"))
        e = is_emt(a.columns)
        for k, v in a.iterrows():
            v = v.values.astype(float)
            if np.nanstd(v) == 0: continue
            p = stats.ttest_ind(v[e], v[~e], equal_var=False).pvalue
            lfc = np.log2(np.mean(2 ** v[e]) / np.mean(2 ** v[~e])) if n.endswith("mean") else v[e].mean() - v[~e].mean()
            rows.append(dict(matrix=n, kinase=k, log2FC_EMT_vs_Met=lfc, P=p))
    D = pd.DataFrame(rows)
    D["assigned"] = np.where((D.P < 0.05) & (D.log2FC_EMT_vs_Met.abs() > 1), np.where(D.log2FC_EMT_vs_Met > 0, "EMT", "Metabolism"), "")
    D.to_csv(os.path.join(TAB, "06_PP_kinase_differential_all_matrices.csv"), index=False)
    cnt = D[D.assigned != ""].groupby(["kinase", "assigned"]).size().unstack(fill_value=0)
    for c in ["EMT", "Metabolism"]:
        if c not in cnt: cnt[c] = 0
    pe = ["MTOR", "CDK7", "LATS1", "PIK3CB", "PKN2", "SRPK1", "ULK1", "AKT3", "MARK2", "PDGFRB", "PIK3CA", "PRKACB", "RIPK3", "AKT1",
          "CIT", "CLK1", "CLK2", "IRAK1", "MAPK1", "MAPK3", "PAK1", "ROCK1", "ROCK2", "RPS6KB1", "STK17A", "STK4"]
    pm = ["STK3", "MKNK2", "STK4", "MAP2K2", "EIF2AK1", "MAP3K15"]
    cnt["paper_EMT"] = cnt.index.isin(pe); cnt["paper_Met"] = cnt.index.isin(pm)
    cnt = cnt.sort_values(by=["EMT", "Metabolism"], ascending=False)
    cnt.to_csv(os.path.join(TAB, "06_Fig8A_assignment_counts.csv"))
    chk = pd.DataFrame(dict(kinase=sorted(set(pe) | set(pm))))
    chk["paper_group"] = ["both" if k in pe and k in pm else ("EMT" if k in pe else "Metabolism") for k in chk.kinase]
    chk["n_EMT"] = [int(cnt.EMT.get(k, 0)) for k in chk.kinase]; chk["n_Met"] = [int(cnt.Metabolism.get(k, 0)) for k in chk.kinase]
    chk["in_n_matrices"] = [sum(k in a.index for a in act.values()) for k in chk.kinase]
    chk["consistent"] = [(r.paper_group in ("EMT", "both") and r.n_EMT > 0) or (r.paper_group in ("Metabolism", "both") and r.n_Met > 0) for r in chk.itertuples()]
    chk.to_csv(os.path.join(TAB, "06_Fig8A_paper_kinases_check.csv"), index=False); print(chk.to_string())
    print(cnt.head(30).to_string())
    top = cnt[(cnt.EMT + cnt.Metabolism) > 0].head(35)
    fig, ax = plt.subplots(figsize=(7, 8)); y = np.arange(len(top))[::-1]
    ax.barh(y, -top.EMT, color="#e8312f"); ax.barh(y, top.Metabolism, color="#1f3c5a")
    ax.set_yticks(y); ax.set_yticklabels([k + (" *" if (r.paper_EMT or r.paper_Met) else "") for k, r in top.iterrows()], fontsize=7)
    ax.axvline(0, color="k"); ax.set_xlabel("times assigned (of 6)  EMT <- | -> Metabolism"); ax.set_title("PP kinase activity (* also in paper Fig. 8A)")
    savefig(fig, "06_Fig8A_PP_kinase_assignments")
    # Fig. S8A-type heatmap: KSEA on normalised PP, assigned kinases
    a = act["norm_ksea"]; ksel = [k for k in cnt.index if k in a.index][:30]
    o = np.argsort([intlab[line_of(c)] + c for c in a.columns])
    zz = a.loc[ksel]; zz = zz.sub(zz.mean(1), axis=0).div(zz.std(1) + 1e-9, axis=0)
    fig, ax = plt.subplots(figsize=(10, 7)); ax.imshow(zz.values[:, o], cmap="cividis", vmin=-2, vmax=2, aspect="auto")
    ax.set_yticks(range(len(ksel))); ax.set_yticklabels(ksel, fontsize=6); ax.set_xticks(range(len(o)))
    ax.set_xticklabels([f"{a.columns[i]} ({intlab[line_of(a.columns[i])][:3]})" for i in o], rotation=90, fontsize=5)
    ax.set_title("KSEA activity, protein-normalised PP"); savefig(fig, "06_FigS8A_KSEA_norm_heatmap")

# ====================================================================== INKA
if step in ("all", "inka"):
    kd = pd.read_csv(os.path.join(DB, "PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv"), keep_default_na=False)
    kd["site_key"] = kd["SUB_GENE"] + "-" + kd["SUB_MOD_RSD"]
    psp = kd[kd.Source == "PhosphoSitePlus"][["GENE", "site_key"]].rename(columns={"GENE": "kinase"}).drop_duplicates()
    nw = kd[kd.Source == "NetworKIN"].copy(); nw["networkin_score"] = pd.to_numeric(nw.networkin_score, errors="coerce")
    nw = nw.sort_values("networkin_score", ascending=False).drop_duplicates("site_key")   # top-scoring kinase per site
    nwk = nw[["GENE", "site_key"]].rename(columns={"GENE": "kinase"})
    al = pd.read_csv(os.path.join(DB, "phomics_global_output.txt"), sep="\t", keep_default_na=False)
    loop = set()
    for g, s in zip(al["Gene Name"], al["Loop Sites"]):
        for p in str(s).split(";"):
            if p.strip().isdigit():
                loop.add((g.split()[0], "Y" + p.strip())); loop.add((g.split()[0], "S" + p.strip())); loop.add((g.split()[0], "T" + p.strip()))
    keys = (ya["GeneSymbol"] + "-" + ya["Residue"]).values
    inten = 2 ** YP.values
    comp = NK.inka_scores(inten, ya["GeneSymbol"].values, keys, KIN, loop, psp, nwk)
    for n, d in comp.items():
        d.columns = YP.columns; d.to_csv(os.path.join(DC, f"inka_{n}.csv"))
    S = comp["inka"]; S = S[(S > 0).any(axis=1)]
    tk = [k for k in S.index if group_of(k) == "TK"]
    St = S.loc[tk]
    St.to_csv(os.path.join(TAB, "06_INKA_scores_TK.csv"))
    top10 = {c: list(St[c].sort_values(ascending=False).index[:10]) for c in St.columns}
    freq = pd.Series([k for v in top10.values() for k in v]).value_counts()
    act_ = np.log2(St + 1).mean(axis=1)
    F = pd.DataFrame(dict(frequency_top10=freq, mean_log2_INKA=act_)).fillna(0).sort_values("frequency_top10", ascending=False)
    F.to_csv(os.path.join(TAB, "06_Fig7A_INKA_frequency_activity.csv")); print(F.head(15))
    paper_top = ["SRC", "ABL1", "EGFR", "EPHA2", "PTK2", "MET"]
    json.dump(dict(top6_this_run=list(F.index[:6]), paper_tier=paper_top, overlap=len(set(F.index[:6]) & set(paper_top)),
                   SRC_rank_by_frequency=int(list(F.index).index("SRC") + 1) if "SRC" in F.index else None,
                   n_TK_scored=len(St)), open(os.path.join(TAB, "06_Fig7A_summary.json"), "w"), indent=1)
    fig, ax = plt.subplots(figsize=(7, 6)); ax.scatter(F.mean_log2_INKA, F.frequency_top10, s=12, color="grey")
    for k in F.index[:15]:
        ax.annotate(k, (F.mean_log2_INKA[k], F.frequency_top10[k]), fontsize=7, color="red" if k in paper_top else "k")
    ax.set_xlabel("mean log2(INKA score + 1)"); ax.set_ylabel("# samples with kinase in top 10 (of 28)")
    ax.set_title("Tyrosine kinases, INKA (self-implemented); red = paper's first tier"); savefig(fig, "06_Fig7A_INKA_frequency_activity")
    src = pd.DataFrame({"mean": [np.log2(S.loc["SRC", [c for c in S.columns if line_of(c) == l]] + 1).mean() for l in LINES],
                        "min": [np.log2(S.loc["SRC", [c for c in S.columns if line_of(c) == l]] + 1).min() for l in LINES],
                        "max": [np.log2(S.loc["SRC", [c for c in S.columns if line_of(c) == l]] + 1).max() for l in LINES]}, index=LINES)
    paper_src_order = ["KatoIII", "N87", "Ocum1", "TMK1", "SNU719", "AZ521", "NCC24", "YCC3", "SNU1750", "YCCEL1", "MKN28", "AGS", "SNU601", "eAGS"]
    src["paper_rank"] = [paper_src_order.index(l) + 1 for l in src.index]; src["this_rank"] = src["mean"].rank().astype(int)
    src.to_csv(os.path.join(TAB, "06_Fig7D_SRC_INKA.csv"))
    rho = stats.spearmanr(src.paper_rank, src.this_rank)
    json.dump(dict(spearman_rho=rho[0], P=rho[1]), open(os.path.join(TAB, "06_Fig7D_SRC_rank_agreement.json"), "w"), indent=1)
    print("SRC rank agreement with Fig. 7D:", rho)
    o = src.sort_values("mean")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(range(14), o["mean"], yerr=[o["mean"] - o["min"], o["max"] - o["mean"]], color=["#e8312f" if intlab[l] == "EMT" else "#1f3c5a" for l in o.index])
    ax.set_xticks(range(14)); ax.set_xticklabels(o.index, rotation=60); ax.set_ylabel("log2(SRC INKA + 1)")
    ax.set_title(f"SRC INKA per line (range of 2 replicates); Spearman vs paper order rho={rho[0]:.2f}"); savefig(fig, "06_Fig7D_SRC_INKA")
    rad = {c: list(St[c].sort_values(ascending=False).index[:10]) for c in St.columns}
    pd.DataFrame(rad).to_csv(os.path.join(TAB, "06_FigS7A_top10_TK_per_sample.csv"), index=False)

# ====================================================================== PTM-SEA
if step in ("all", "ptmsea"):
    seqs, cur = {}, None
    for line in open(os.path.join(DC, "uniprot", "human_sp_isoforms.fasta")):
        if line.startswith(">"): cur = line.split("|")[1]; seqs[cur] = []
        else: seqs[cur].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}

    def win(acc, res):
        m = re.match(r"([STY])(\d+)$", str(res))
        if not m: return None
        aa, pos = m.group(1), int(m.group(2))
        for a in [acc, str(acc).split("-")[0]]:
            s = seqs.get(a)
            if s and pos <= len(s) and s[pos - 1] == aa:
                return ("_" * 7 + s + "_" * 7)[pos - 1:pos + 14]
        return None
    w = [win(a, r) for a, r in zip(pa["UniprotID"], pa["Residue"])]
    ok = np.array([x is not None for x in w])
    X = PP[ok].sub(PP[ok].mean(1), axis=0)
    ids = [x + "-p" for x in np.array(w, dtype=object)[ok]]
    X = X.groupby(ids).mean()
    print("PP sites with windows:", ok.sum(), "unique windows:", len(X))
    sig = M.read_gmt(os.path.join(DB, "ptm.sig.db.all.flanking.human.v1.9.1.gmt"))
    names, ES, NES, Pp = M.ssgsea_scores(X.values, list(X.index), sig, weight=0.75, nperm=1000, min_overlap=10, seed=1)
    N = pd.DataFrame(NES, index=names, columns=X.columns); N.to_csv(os.path.join(TAB, "06_PTMSEA_NES.csv"))
    pd.DataFrame(Pp, index=names, columns=X.columns).to_csv(os.path.join(TAB, "06_PTMSEA_P.csv"))
    e = is_emt(N.columns).astype(float)
    r, _ = M.moderated_t(N.values, np.column_stack([np.ones(len(e)), e]), [0, 1]); r.index = names
    r = r.sort_values("P"); r.to_csv(os.path.join(TAB, "06_PTMSEA_EMT_vs_Met.csv"))
    s5d = ["KINASE-PSP_HIPK2", "KINASE-PSP_PKCA/PRKCA", "PATH-NP_PROLACTIN_PATHWAY", "PATH-NP_EGFR1_PATHWAY", "KINASE-PSP_Src/SRC",
           "PATH-NP_CCR7_PATHWAY", "KINASE-PSP_JNK1/MAPK8", "PATH-NP_TSLP_PATHWAY", "KINASE-PSP_p70S6K/RPS6KB1", "KINASE-PSP_CDK2",
           "KINASE-PSP_CDK1", "PATH-NP_GLP1_PATHWAY", "PATH-NP_B_CELL_RECEPTOR_PATHWAY", "KINASE-PSP_CDK4", "KINASE-PSP_CDK6",
           "KINASE-PSP_CDK5", "KINASE-PSP_ERK1/MAPK3", "KINASE-PSP_P38B/MAPK11"]
    chk = r.reindex(s5d); chk["tested"] = chk.t.notna(); chk.to_csv(os.path.join(TAB, "06_FigS5D_paper_signatures_check.csv"))
    print(chk.to_string())
    sel = [s for s in s5d if s in N.index]
    o = np.argsort([intlab[line_of(c)] + c for c in N.columns])
    zz = N.loc[sel]; zz = zz.sub(zz.mean(1), axis=0).div(zz.std(1) + 1e-9, axis=0)
    fig, ax = plt.subplots(figsize=(10, 6)); ax.imshow(zz.values[:, o], cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel, fontsize=6); ax.set_xticks(range(len(o)))
    ax.set_xticklabels([f"{N.columns[i]} ({intlab[line_of(N.columns[i])][:3]})" for i in o], rotation=90, fontsize=5)
    ax.set_title("PTM-SEA NES (z-scored) for the paper's Fig. S5D signatures"); savefig(fig, "06_FigS5D_PTMSEA")
