"""Targeted analyses:
 (1) protein-level 'amplicon' scan by cytoband (positive controls ERBB2/N87, FGFR2/KatoIII)
 (2) lineage transcription-factor programme (HNF1A/HNF4A/FOXA2) vs the published subtypes and
     its relation to the oxidative / fatty-acid metabolism protein score
 (3) the EBV-like proteome module (ER / N-glycosylation up, transcriptional co-activators down)
 (4) marker panel heatmap
"""
import json
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from common import (CAT, CELL_META, DATA, INK2, MUTED, SEED, SUBTYPE_COLORS, TAB, bh,
                    diverging_cmap, read_gmt, savefig, style)

style()
rng = np.random.default_rng(SEED)
out = {}
L = pd.read_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
gene = L.pop("gene").fillna("")
sym = gene.str.split(";").str[0]
lines = list(L.columns)
sub = CELL_META.paper_subtype.reindex(lines)
EBV_NAT = ["NCC24", "SNU719", "YCCEL1"]
ORDER = ["AGS", "eAGS", "AZ521", "SNU1750", "MKN28", "N87", "TMK1",
         "NCC24", "SNU719", "YCCEL1", "KatoIII", "Ocum1", "SNU601", "YCC3"]

# --------------------------------------------------------------------------
# (1) amplicon scan: per-line z-score vs other lines, aggregated by cytoband
# --------------------------------------------------------------------------
hg = pd.read_csv(DATA / "genesets/hgnc_complete_set.txt", sep="\t", low_memory=False,
                 usecols=["symbol", "location", "locus_group"])
hg = hg[hg.locus_group == "protein-coding gene"].drop_duplicates("symbol").set_index("symbol")
band = sym.map(hg.location)
Lp = L[L.notna().sum(axis=1) >= 10]
z = {}
for ln in lines:
    others = Lp.drop(columns=ln)
    z[ln] = (Lp[ln] - others.median(axis=1)) / (others.std(axis=1) + 0.1)
Z = pd.DataFrame(z)
Z["band"] = band.reindex(Z.index)
Z["arm_band"] = Z.band.str.extract(r"^(\d+|X|Y)([pq]\d+)")[0] + \
    Z.band.str.extract(r"^(?:\d+|X|Y)([pq]\d+)")[0]
amp = []
for ln in lines:
    g = Z.dropna(subset=["arm_band"]).groupby("arm_band")[ln]
    stat = g.agg(["mean", "count"])
    stat = stat[stat["count"] >= 3]
    for b, row in stat.iterrows():
        members = Z[(Z.arm_band == b)][ln].sort_values(ascending=False)
        hi = members[members > 3]
        if len(hi) >= 3:
            amp.append((ln, b, int(row["count"]), len(hi), float(row["mean"]),
                        ", ".join(sym[hi.index][:12])))
amp = pd.DataFrame(amp, columns=["line", "cytoband", "proteins_in_band", "n_proteins_z>3",
                                 "mean_z_band", "proteins_z>3"]).sort_values("n_proteins_z>3", ascending=False)
amp.to_csv(TAB / "protein_outlier_cytobands.csv", index=False)
print(amp.head(25).to_string())
out["cytoband_hotspots_top"] = amp.head(15).values.tolist()

fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"wspace": 0.3})
for ax, ln, gs, title in [(axes[0], "N87", ["ERBB2", "GRB7", "CDK12", "PGAP3", "MIEN1", "STARD3", "PNMT", "ORMDL3"], "A  N87: 17q12 (ERBB2 amplicon)"),
                          (axes[1], "KatoIII", ["FGFR2", "WDR11", "ATE1", "NSMCE4A", "TACC2", "BAG3", "PLEKHA1", "HTRA1"], "B  KatoIII: 10q26 (FGFR2 amplicon)")]:
    rows = [(g, i) for g in gs for i in L.index[sym == g]]
    data = [(g, L.loc[i, ln], L.loc[i].drop(ln).dropna().values) for g, i in rows]
    for k, (g, v, oth) in enumerate(data):
        ax.scatter(np.full(len(oth), k) + rng.uniform(-0.15, 0.15, len(oth)), oth, s=10,
                   color="#c3c2b7", linewidths=0)
        if pd.notna(v):
            ax.scatter([k], [v], s=36, color=CAT[1], zorder=3, edgecolor="white")
    ax.set_xticks(range(len(data)), [d[0] for d in data], rotation=45)
    ax.set_ylabel("log2 intensity")
    ax.set_title(title)
axes[0].scatter([], [], color=CAT[1], label="index line"); axes[0].scatter([], [], color="#c3c2b7", label="other 13 lines")
axes[0].legend(fontsize=7, loc="lower left")
savefig(fig, "11_amplicon_positive_controls")
plt.close(fig)

# --------------------------------------------------------------------------
# (2) lineage TF programme
# --------------------------------------------------------------------------
floor = np.nanmin(L.values) - 1
Lf = L.fillna(floor)  # undetected -> below detection floor (presence information kept)


def get(g):
    idx = L.index[sym == g]
    return Lf.loc[idx].mean() if len(idx) else None


lineage_genes = ["HNF1A", "HNF4A", "FOXA2", "ANXA10", "LGALS4", "AGR2", "MUC13", "CES2", "CDH17", "VIL1"]
lin = pd.DataFrame({g: get(g) for g in lineage_genes}).T
linz = lin.sub(lin.mean(axis=1), axis=0).div(lin.std(axis=1), axis=0)
TFscore = linz.loc[["HNF1A", "HNF4A", "FOXA2"]].mean()
# metabolic score from hallmarks (computed in 02)
hs = pd.read_csv(TAB / "hallmark_scores_per_line.csv", index_col=0)
met_score = hs.loc[["FATTY_ACID_METABOLISM", "OXIDATIVE_PHOSPHORYLATION", "PEROXISOME", "ADIPOGENESIS"]].mean()
r_p, p_p = pearsonr(TFscore[lines], met_score[lines])
r_s, p_s = spearmanr(TFscore[lines], met_score[lines])
# within-group partial check: correlation inside each published subtype
within = {s: spearmanr(TFscore[sub[sub == s].index], met_score[sub[sub == s].index]) for s in ("EMT", "Metabolism")}
out["lineageTF_vs_metabolic_score"] = {"pearson_r": r_p, "p": p_p, "spearman_rho": r_s, "p_spearman": p_s,
                                       "within_EMT_rho": within["EMT"][0], "within_EMT_p": within["EMT"][1],
                                       "within_Metab_rho": within["Metabolism"][0], "within_Metab_p": within["Metabolism"][1]}
# proteome-wide correlation with HNF1A presence/level and HNF4A
Lc = L[L.notna().sum(axis=1) >= 12]
rows = []
for p in Lc.index:
    v = Lc.loc[p]
    ok = v.notna()
    rows.append((p, sym[p], spearmanr(TFscore[ok.index[ok]], v[ok])[0]))
cor = pd.DataFrame(rows, columns=["protein", "gene", "rho_with_lineageTF"]).set_index("protein")
cor.to_csv(TAB / "correlation_with_lineageTF_score.csv.gz", float_format="%.4f")
from stats import rank_set_test
gsets = {}
for f in ["h.all.v2023.1.Hs.symbols.gmt", "c2.cp.kegg.v7.5.1.symbols.gmt", "c2.cp.reactome.v7.5.1.symbols.gmt",
          "c5.go.bp.v7.5.1.symbols.gmt"]:
    gsets.update(read_gmt(DATA / "genesets" / f, set(sym), 10, 300))
gs = rank_set_test(cor.groupby("gene").rho_with_lineageTF.mean(), gsets)
gs.to_csv(TAB / "geneset_correlation_with_lineageTF.csv", index=False, float_format="%.4g")
out["lineageTF_correlated_sets_pos"] = gs[(gs.q < 0.05) & (gs.z > 0)].head(15).set.tolist()
out["lineageTF_correlated_sets_neg"] = gs[(gs.q < 0.05) & (gs.z < 0)].head(15).set.tolist()
print(json.dumps({k: out[k] for k in ["lineageTF_vs_metabolic_score", "lineageTF_correlated_sets_pos",
                                      "lineageTF_correlated_sets_neg"]}, indent=1, default=float))

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"wspace": 0.35, "width_ratios": [1.4, 1]})
ax = axes[0]
cols = [l for l in ORDER]
M = lin[cols].where(lin[cols] > floor)
v = linz[cols]
im = ax.imshow(v.values, cmap=diverging_cmap(), vmin=-2, vmax=2, aspect="auto")
for i in range(M.shape[0]):
    for jx in range(M.shape[1]):
        if pd.isna(M.iloc[i, jx]):
            ax.text(jx, i, "n.d.", ha="center", va="center", fontsize=5.5, color=INK2)
ax.set_xticks(range(len(cols)), cols, rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_yticks(range(len(lineage_genes)), lineage_genes)
fig.colorbar(im, ax=ax, shrink=0.8, label="z-score (n.d. = not detected,\nset below detection floor)")
ax.set_title("A  GI-lineage transcription factors and targets")
ax = axes[1]
for s, c in SUBTYPE_COLORS.items():
    m = sub == s
    ax.scatter(TFscore[lines][m.values], met_score[lines][m.values], s=40, color=c, label=s, edgecolor="white")
for ln in lines:
    ax.annotate(ln, (TFscore[ln], met_score[ln]), fontsize=7, color=INK2, xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("HNF1A/HNF4A/FOXA2 protein score (z)")
ax.set_ylabel("FAO + OXPHOS + peroxisome +\nadipogenesis hallmark score")
ax.set_title(f"B  Lineage TFs vs metabolic programme\nSpearman ρ = {r_s:.2f} (p = {p_s:.2g}), n = 14 lines")
ax.legend(title="Hu et al.", fontsize=7)
savefig(fig, "12_lineage_TF_programme")
plt.close(fig)

# --------------------------------------------------------------------------
# (3) EBV-like module
# --------------------------------------------------------------------------
d = pd.read_csv(TAB / "diff_clusterEBVlike_vs_rest.csv.gz", index_col=0)
up = d[(d.q < 0.05) & (d.logFC > 0.5)].index
dn = d[(d.q < 0.05) & (d.logFC < -0.5)].index
Lz = L.sub(L.mean(axis=1), axis=0).div(L.std(axis=1), axis=0)
mod_score = Lz.loc[up].mean() - Lz.loc[dn].mean()
# leave-TMK1-out: derive module only from natural EBV lines vs EBV-negative lines (excluding TMK1, eAGS)
dn2 = pd.read_csv(TAB / "diff_EBVnatural_vs_EBVneg.csv.gz", index_col=0)
from stats import moderated_t
neg_noTMK = [l for l in lines if l not in EBV_NAT + ["TMK1", "eAGS"]]
r2 = moderated_t(L, EBV_NAT, neg_noTMK, min_per_group=3)
up2 = r2[(r2.q < 0.05) & (r2.logFC > 0.5)].index
dn2_ = r2[(r2.q < 0.05) & (r2.logFC < -0.5)].index
mod2 = Lz.loc[up2].mean() - Lz.loc[dn2_].mean()
# permutation: how extreme is TMK1 vs. the 10 EBV-negative lines on a module it did not help define?
null = mod2[neg_noTMK]
out["EBV_module_TMKout"] = {"n_up": len(up2), "n_down": len(dn2_),
                            "score": mod2.round(3).to_dict(),
                            "TMK1_rank_among_EBVneg_plus_TMK1": int((mod2[neg_noTMK + ["TMK1"]].rank(ascending=False))["TMK1"]),
                            "TMK1_z_vs_EBVneg": float((mod2["TMK1"] - null.mean()) / null.std())}
# fair null: leave-one-out. For each EBV-negative line (and TMK1), the module is re-derived
# WITHOUT that line (and without TMK1/eAGS), then the held-out line is scored.
loo = {}
for held in neg_noTMK + ["TMK1"]:
    train_neg = [l for l in neg_noTMK if l != held]
    rr = moderated_t(L, EBV_NAT, train_neg, min_per_group=3)
    u = rr[(rr.q < 0.05) & (rr.logFC > 0.5)].index
    dd = rr[(rr.q < 0.05) & (rr.logFC < -0.5)].index
    loo[held] = {"score": float(Lz.loc[u, held].mean() - Lz.loc[dd, held].mean()),
                 "n_up": len(u), "n_down": len(dd),
                 "EBVnat_mean": float((Lz.loc[u, EBV_NAT].mean() - Lz.loc[dd, EBV_NAT].mean()).mean())}
loo = pd.DataFrame(loo).T
loo.to_csv(TAB / "EBV_module_leave_one_out.csv")
negs = loo.loc[neg_noTMK, "score"]
out["EBV_module_LOO"] = {"TMK1_score": float(loo.loc["TMK1", "score"]),
                         "heldout_EBVneg_scores": negs.round(3).to_dict(),
                         "TMK1_exceeds_all_heldout_EBVneg": bool((loo.loc["TMK1", "score"] > negs).all()),
                         "TMK1_z_vs_heldout_EBVneg": float((loo.loc["TMK1", "score"] - negs.mean()) / negs.std()),
                         "EBVnat_mean_score_range": [float(loo.EBVnat_mean.min()), float(loo.EBVnat_mean.max())]}
print(json.dumps(out["EBV_module_TMKout"], indent=1))
print(json.dumps(out["EBV_module_LOO"], indent=1))
r2.insert(0, "gene", gene.reindex(r2.index))
r2.sort_values("p").to_csv(TAB / "diff_EBVnatural_vs_EBVneg_excl_TMK1.csv.gz", float_format="%.5g")

# GO-CC/ Reactome composition of the module
er_sets = {k: v for k, v in gsets.items() if k in (
    "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION", "GOBP_ERAD_PATHWAY",
    "GOBP_ENDOPLASMIC_RETICULUM_TO_GOLGI_VESICLE_MEDIATED_TRANSPORT",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_MRNA_SPLICING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")}
sc = pd.DataFrame({k: Lz.loc[L.index[sym.isin(v)]].mean() for k, v in er_sets.items()}).T

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"wspace": 0.45, "width_ratios": [1, 1.5]})
ax = axes[0]
vals = mod2[ORDER]
colors = [CAT[6] if l in EBV_NAT else (CAT[4] if l == "TMK1" else (CAT[3] if l == "eAGS" else "#c3c2b7")) for l in ORDER]
ax.bar(range(len(ORDER)), vals.values, color=colors, width=0.75)
ax.axhline(0, color="#c3c2b7", lw=0.6)
ax.set_xticks(range(len(ORDER)), ORDER, rotation=90)
ax.set_ylabel("module score (z up − z down)")
ax.set_title("A  Natural-EBV module (defined without TMK1)\n%d up / %d down proteins" % (len(up2), len(dn2_)))
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=CAT[6], label="natural EBV+ (Table S1)"), Patch(color=CAT[4], label="TMK1 (EBV− per Table S1)"),
                   Patch(color=CAT[3], label="eAGS (in-vitro EBV)"), Patch(color="#c3c2b7", label="EBV−")],
          fontsize=7, loc="upper left")
ax = axes[1]
im = ax.imshow(sc[ORDER].values, cmap=diverging_cmap(), vmin=-1.2, vmax=1.2, aspect="auto")
ax.set_xticks(range(len(ORDER)), ORDER, rotation=90)
ax.set_yticks(range(len(sc)), [s.replace("REACTOME_", "R_").replace("HALLMARK_", "H_").replace("GOBP_", "GO_")[:45] for s in sc.index], fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.8, label="mean protein z-score")
ax.set_title("B  Pathway scores across lines")
savefig(fig, "13_EBV_module")
plt.close(fig)

# --------------------------------------------------------------------------
# (4) marker heatmap
# --------------------------------------------------------------------------
panel = {"epithelial": ["CDH1", "EPCAM", "KRT8", "KRT18", "CLDN4"],
         "mesenchymal": ["VIM", "CDH2", "FN1", "FSCN1", "ZEB1"],
         "gastric/foveolar": ["MUC5AC", "TFF1", "ANXA10", "CLDN18", "LGALS4", "AGR2"],
         "intestinal": ["CDX2", "CDH17", "VIL1", "MUC13"],
         "lineage TF": ["HNF1A", "HNF4A", "FOXA2", "GATA6", "KLF5", "ELF3", "SOX2"],
         "EBV-GC related": ["CDKN2A", "CD274", "HLA-A", "HLA-B", "B2M", "STAT1", "ARID1A", "EP300"],
         "drivers/targets": ["ERBB2", "GRB7", "CDK12", "FGFR2", "MET", "EGFR", "SRC", "MTOR"]}
rows, labels, grp = [], [], []
for gname, gl in panel.items():
    for g in gl:
        idx = L.index[sym == g]
        if len(idx):
            i = L.loc[idx].notna().sum(axis=1).idxmax()
            rows.append(L.loc[i, ORDER]); labels.append(g); grp.append(gname)
H = pd.DataFrame(rows, index=labels)
Hz = H.sub(H.mean(axis=1), axis=0).div(H.std(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(8.5, 12))
im = ax.imshow(Hz.values, cmap=diverging_cmap(), vmin=-2.2, vmax=2.2, aspect="auto")
for i in range(H.shape[0]):
    for jx in range(H.shape[1]):
        if pd.isna(H.iloc[i, jx]):
            ax.add_patch(plt.Rectangle((jx - 0.5, i - 0.5), 1, 1, facecolor="white",
                                       edgecolor="#c3c2b7", hatch="////", lw=0))
ax.set_xticks(range(len(ORDER)), ORDER, rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_yticks(range(len(labels)), labels, fontsize=7.5)
ax.tick_params(axis="y", pad=2)
prev = None
for i, g in enumerate(grp):
    if g != prev:
        ax.axhline(i - 0.5, color="white", lw=2)
        ax.text(-4.2, i - 0.45, g, fontsize=7, color=INK2, va="top", fontweight="bold")
        prev = g
fig.colorbar(im, ax=ax, shrink=0.5, orientation="horizontal", pad=0.08,
             label="z-score across lines (hatched = not detected)")
ax.set_title("Marker proteins across the 14 lines\n(label colour: Hu et al. subtype)")
savefig(fig, "14_marker_heatmap")
plt.close(fig)
H.to_csv(TAB / "marker_panel_log2.csv", float_format="%.3f")
json.dump(out, open(TAB / "targeted_summary.json", "w"), indent=2, default=float)
