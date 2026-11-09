"""Shared paths, metadata and helpers for the PXD048538 (quantms e60b9e9aa4c8) analysis."""
import os
from pathlib import Path

import numpy as np
import pandas as pd

RUN_DIR = Path(__file__).resolve().parent.parent
RES = RUN_DIR / "results"
FIG = RES / "figures"
TAB = RES / "tables"
for p in (FIG, TAB):
    p.mkdir(parents=True, exist_ok=True)

# Data directory holding the downloaded files (see 00_download.sh); not committed.
DATA = Path(os.environ.get("PXD048538_DATA", "data_cache"))
QDIR = DATA / "quantms"

SEED = 20261105

# --------------------------------------------------------------------------
# Cell-line annotation.
# * histology_sdrf  : quantms SDRF (characteristics[disease])
# * histology_s1, ebv, msi, alterations : Hu et al. 2024 (Mol Oncol) Table S1
# * paper_subtype   : integrative subtype, Hu et al. 2024 section 3.5
# * paper_wp_cluster: not given per line in the text; omitted
# --------------------------------------------------------------------------
CELL_META = pd.DataFrame(
    [
        # line,   sdrf name,        histology (S1),             EBV,        MSI,       key alterations (S1),                                   paper subtype
        ("AGS",     "AGS",            "Poorly differentiated",   "neg",      "MSS",     "PIK3CA E453K/E545A; KRAS G12D; CTNNB1 G34E",            "EMT"),
        ("eAGS",    "AGS-EBV (eAGS)", "Poorly differentiated",   "acquired", "NA",      "AGS + in-vitro EBV infection",                          "EMT"),
        ("AZ521",   "AZ-521",         "Duodenal adenocarcinoma", "neg",      "NA",      "CTNNB1 S37F",                                           "EMT"),
        ("KatoIII", "KATO III",       "Signet ring cell",        "neg",      "MSS",     "FGFR2 amp; NOTCH1 W2532*",                              "Metabolism"),
        ("MKN28",   "MKN28",          "Well differentiated",     "neg",      "MSS",     "RASA1 L761fs; APC R1450*; TP53 I251L",                  "EMT"),
        ("N87",     "NCI-N87",        "Well differentiated",     "neg",      "MSS",     "ERBB2 amp; TP53 R248Q; ELF3 S308fs",                    "EMT"),
        ("NCC24",   "SNU-NCC-24",     "Signet ring cell",        "pos",      "MSS",     "TP53 G266R/R273H",                                      "Metabolism"),
        ("Ocum1",   "OCUM-1",         "Signet ring cell",        "neg",      "MSI-low", "TGFBR2 R537C; CTNNB1 S45C; MAP2K1 Q56P",                "Metabolism"),
        ("SNU601",  "SNU-601",        "Signet ring cell",        "neg",      "MSS",     "PIK3CA E542K; KRAS G12D; TP53 R273H; NF1 E2174*",       "Metabolism"),
        ("SNU719",  "SNU-719",        "Well differentiated",     "pos",      "NA",      "NRAS Q61L; CTNNB1 G34V; PIK3CA P104R; RB1 R621S",       "Metabolism"),
        ("SNU1750", "SNU-1750",       "Poorly differentiated",   "neg",      "MSS",     "TP53 R273P",                                            "EMT"),
        ("TMK1",    "TMK-1",          "Poorly differentiated",   "neg",      "MSS",     "TP53 V173M",                                            "EMT"),
        ("YCC3",    "YCC-3",          "NA",                      "neg",      "NA",      "TP53 R175H",                                            "Metabolism"),
        ("YCCEL1",  "YCCEL1",         "NA",                      "pos",      "NA",      "-",                                                     "Metabolism"),
    ],
    columns=["line", "sdrf_name", "histology_s1", "ebv", "msi", "alterations_s1", "paper_subtype"],
).set_index("line")

SUBTYPE_COLORS = {"EMT": "#eb6834", "Metabolism": "#2a78d6"}
# categorical slots (validated reference palette, light mode)
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"


def diverging_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "bgr", ["#104281", "#3987e5", "#cde2fb", "#f0efec", "#f9d0c4", "#e34948", "#8f1d1d"])


def style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white", "axes.titlesize": 10,
        "axes.titleweight": "bold", "axes.titlelocation": "left"})


def short_line(sdrf_name):
    m = {v: k for k, v in CELL_META["sdrf_name"].items()}
    return m[sdrf_name]


def read_gmt(path, universe=None, min_size=10, max_size=500):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes = set(parts[2:])
            if universe is not None:
                genes &= universe
            if min_size <= len(genes) <= max_size:
                sets[parts[0]] = genes
    return sets


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.minimum(ranked, 1)
    out[ok] = res
    return out


def savefig(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
