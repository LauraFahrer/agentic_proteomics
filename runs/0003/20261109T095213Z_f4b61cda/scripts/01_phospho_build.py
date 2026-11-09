"""Build the (phospho)proteome matrices and protein-normalised phosphosite values.

Source: publication Supplementary Table S3 (author-processed PP/YP/WP matrices); the quantms
e60b9e9aa4c8 release itself has no PTM layer (stated in the report).
Outputs: tables/pp_linemean.csv.gz, yp_linemean.csv.gz, pp_protnorm_linemean.csv.gz,
         phospho_site_meta.csv.gz, phospho_coverage.json + coverage figure.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CAT, CELL_META, SUBTYPE_COLORS, TAB, savefig, style
from phospho_common import load_PP, load_WP_linemean, load_YP

style()
order = list(CELL_META.index)
ppm, ppx = load_PP(order)
ypm, ypx = load_YP(order)
wp_acc, wp_gene = load_WP_linemean(order)
cov = {}

# keep sites whose 10-mer central residue matches the annotated residue type
ppm["win_ok"] = ppm.win10.notna() & (ppm.win10.str[5] == ppm.restype)
cov["PP_sites"] = int(len(ppm)); cov["PP_sites_window_ok"] = int(ppm.win_ok.sum())
cov["PP_restype"] = ppm.restype.value_counts().to_dict()
cov["YP_sites"] = int(len(ypm)); cov["YP_window_ok"] = int(ypm.win10.notna().sum())

# protein-normalise PP: phosphosite(line) - parent protein(line), matched by UniProt accession
prot = wp_acc.reindex(ppm.UniprotID.values)
prot.index = ppm.index
has_prot = prot.notna().all(axis=1)
ppnorm = ppx.copy()
ppnorm[has_prot.values] = (ppx[has_prot.values].values - prot[has_prot.values].values)
ppnorm[~has_prot.values] = np.nan
cov["PP_sites_protein_normalisable"] = int(has_prot.sum())

# save
meta_cols = ["site", "GeneSymbol", "UniprotID", "pos", "restype", "win10", "win_ok"]
ppm_out = ppm[meta_cols].copy()
ppm_out.to_csv(TAB / "phospho_site_meta_PP.csv.gz", index=False)
ypm[["site", "GeneSymbol", "UniprotID", "pos", "restype", "win10"]].to_csv(TAB / "phospho_site_meta_YP.csv.gz", index=False)
ppx.assign(site=ppm.site.values).to_csv(TAB / "pp_linemean.csv.gz", index=False, float_format="%.4f")
ypx.assign(site=ypm.site.values).to_csv(TAB / "yp_linemean.csv.gz", index=False, float_format="%.4f")
ppnorm.assign(site=ppm.site.values, prot_normalised=has_prot.values).to_csv(
    TAB / "pp_protnorm_linemean.csv.gz", index=False, float_format="%.4f")
json.dump(cov, open(TAB / "phospho_coverage.json", "w"), indent=2)
print(json.dumps(cov, indent=2))

# coverage figure: sites quantified per line + S/T/Y composition
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={"wspace": 0.3})
ax = axes[0]
ppn = ppx.notna().sum(); ypn = ypx.notna().sum()
x = np.arange(len(order)); w = 0.4
ax.bar(x - w / 2, [ppn.get(l, 0) for l in order], width=w, color=CAT[0], label="phosphoproteome (S/T/Y)")
ax.bar(x + w / 2, [ypn.get(l, 0) for l in order], width=w, color=CAT[1], label="tyrosine (pY) profiling")
ax.set_xticks(x, order, rotation=90)
for t, l in zip(ax.get_xticklabels(), order):
    t.set_color(SUBTYPE_COLORS[CELL_META.paper_subtype[l]])
ax.set_ylabel("phosphosites quantified"); ax.legend(fontsize=8)
ax.set_title("A  Phosphosite coverage per cell line")
ax = axes[1]
comp = ppm.restype.value_counts()
ax.bar(["Ser", "Thr", "Tyr\n(global)"], [comp.get("S", 0), comp.get("T", 0), comp.get("Y", 0)], color=CAT[2])
ax.bar(["pY\n(IP)"], [len(ypm)], color=CAT[1])
ax.set_ylabel("phosphosites"); ax.set_title("B  Modified-residue composition")
for i, v in enumerate([comp.get("S", 0), comp.get("T", 0), comp.get("Y", 0), len(ypm)]):
    ax.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
savefig(fig, "01_phospho_coverage")
plt.close(fig)
