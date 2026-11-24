#!/usr/bin/env python3
"""Step 06b - sensitivity of the PP kinase-activity assignment (Fig. 8A) to two implementation details that
the paper does not report but that the original PhosMap code fixes:
  (1) regulon construction: exact kinase-symbol matching (our main analysis) vs PhosMap's substring matching
      (R grepl(kinase, "kinase|site|substrate") -> e.g. 'CDK1' also collects CDK12/CDK13 substrates and any entry
      whose substrate symbol contains 'CDK1'),
  (2) input scaling: per-site centred log2 (main) vs uncentred log2 intensities.
Methods: mean and MLR (self-implemented, networks_kinases); KSEA for the substring/centred variant.
Same assignment rule: Welch t-test P < 0.05 and |log2FC| > 1, counted over raw and normalised PP."""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from common import *
import networks_kinases as NK

do_ksea = "--ksea" in sys.argv
INT = pd.read_csv(os.path.join(TAB, "03_integrative_clusters.csv")); intlab = dict(zip(INT.line, INT.integrative_subtype))
WP, wa = load_WP(); PP, pa = load_PP()
ks = pd.read_csv(os.path.join(DB, "human_ksrr.csv"))
P = PP.copy(); P.index = pa["GeneSymbol"] + "_" + pa["Residue"]
wpg = pd.Series(range(len(WP)), index=wa["Gene"].values); wpg = wpg[~wpg.index.duplicated()]
has = P.index.str.split("_").str[0].isin(wpg.index)
Pn = P[has] - WP.iloc[wpg[P.index[has].str.split("_").str[0]].values].values


def regulons_substring(site_ids, ks):
    """PhosMap get_ka_by_mean_or_mlr regulon logic: collect 'kinase|site|substrate' strings of all detected sites,
    then for each kinase K take every string matching grepl(K, .) (regular-expression substring)."""
    key = pd.Series(range(len(site_ids)), index=site_ids)
    k2 = ks.assign(sid=ks["substrate"].astype(str) + "_" + ks["site"].astype(str))
    k2 = k2[k2.sid.isin(key.index)]
    strings = (k2.kinase.astype(str) + "|" + k2.site.astype(str) + "|" + k2.substrate.astype(str)).values
    sids = k2.sid.values
    reg = {}
    for kin in ks.kinase.astype(str).unique():
        hit = pd.Series(strings).str.contains(kin, regex=True).values
        if hit.any():
            reg[kin] = np.unique(key[sids[hit]].values)
    return reg


def assign(act, name):
    e = np.array([intlab[line_of(c)] == "EMT" for c in PP.columns]); rows = []
    for k, v in act.iterrows():
        v = v.values.astype(float)
        if np.nanstd(v) == 0: continue
        p = stats.ttest_ind(v[e], v[~e], equal_var=False).pvalue
        lfc = np.log2(np.mean(2 ** v[e]) / np.mean(2 ** v[~e])) if "mean" in name else v[e].mean() - v[~e].mean()
        rows.append(dict(variant=name, kinase=k, log2FC=lfc, P=p,
                         assigned="EMT" if (p < 0.05 and lfc > 1) else ("Metabolism" if (p < 0.05 and lfc < -1) else "")))
    return pd.DataFrame(rows)


out = []
for match in ["exact", "substring"]:
    for scaling in ["centred", "uncentred"]:
        for lvl, m in [("raw", P), ("norm", Pn)]:
            X = m.sub(m.mean(axis=1), axis=0) if scaling == "centred" else m
            reg = NK.build_regulons(list(X.index), ks) if match == "exact" else regulons_substring(list(X.index), ks)
            tag = f"{match}|{scaling}|{lvl}"
            out.append(assign(NK.activity_mean(X.values, reg).set_axis(X.columns, axis=1), tag + "|mean"))
            out.append(assign(NK.activity_mlr(X.values, reg).set_axis(X.columns, axis=1), tag + "|mlr"))
            if do_ksea and match == "substring" and scaling == "centred":
                out.append(assign(NK.activity_ksea(X.values, reg).set_axis(X.columns, axis=1), tag + "|ksea"))
            print(tag, "kinases:", len(reg), flush=True)
D = pd.concat(out); D.to_csv(os.path.join(TAB, "06b_activity_sensitivity_all.csv"), index=False)
pe = {"MTOR", "CDK7", "LATS1", "PIK3CB", "PKN2", "SRPK1", "ULK1", "AKT3", "MARK2", "PDGFRB", "PIK3CA", "PRKACB", "RIPK3", "AKT1",
      "CIT", "CLK1", "CLK2", "IRAK1", "MAPK1", "MAPK3", "PAK1", "ROCK1", "ROCK2", "RPS6KB1", "STK17A", "STK4"}
pm = {"STK3", "MKNK2", "STK4", "MAP2K2", "EIF2AK1", "MAP3K15"}
D["setting"] = D.variant.str.rsplit("|", n=2).str[0]
summ = []
for st, g in D.groupby("setting"):
    a = g[g.assigned != ""]
    emt, met = set(a[a.assigned == "EMT"].kinase), set(a[a.assigned == "Metabolism"].kinase)
    summ.append(dict(setting=st, n_assigned_EMT=len(emt), n_assigned_Met=len(met),
                     paper_EMT_recovered=len(emt & pe), paper_Met_recovered=len(met & pm),
                     paper_kinases_opposite=len((emt & (pm - pe)) | (met & (pe - pm))),
                     recovered=";".join(sorted((emt & pe) | (met & pm)))))
S = pd.DataFrame(summ); S.to_csv(os.path.join(TAB, "06b_activity_sensitivity_summary.csv"), index=False); print(S.to_string())
