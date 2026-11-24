"""Shared paths, loaders and published reference labels for run 0007."""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DC = os.path.join(RUN, "data_cache")
PROC = os.path.join(DC, "processed")
DB = os.path.join(DC, "db")
FIG = os.path.join(RUN, "results", "figures")
TAB = os.path.join(RUN, "results", "tables")
os.makedirs(FIG, exist_ok=True); os.makedirs(TAB, exist_ok=True)

LINES = ["AGS", "AZ521", "eAGS", "KatoIII", "MKN28", "N87", "NCC24", "Ocum1",
         "SNU1750", "SNU601", "SNU719", "TMK1", "YCC3", "YCCEL1"]
# published assignments (read from the paper's figures; identical to run 0006 common.R)
PAPER_WP = dict(AGS=1, eAGS=1, MKN28=1, N87=1, AZ521=2, SNU1750=2, TMK1=2,
                KatoIII=3, NCC24=3, Ocum1=3, SNU601=3, SNU719=3, YCC3=3, YCCEL1=3)
PAPER_PP = dict(AGS=1, AZ521=1, eAGS=1, MKN28=1, N87=1, SNU1750=1, TMK1=1,
                KatoIII=2, NCC24=2, Ocum1=2, SNU601=2, SNU719=2, YCC3=2, YCCEL1=2)
PAPER_INT = {k: ("EMT" if v == 1 else "Metabolism") for k, v in PAPER_PP.items()}
PAPER_YP = dict(AGS=1, AZ521=1, MKN28=1, SNU1750=1, TMK1=1, YCC3=1, KatoIII=2, N87=2,
                eAGS=3, NCC24=3, Ocum1=3, SNU601=3, SNU719=3, YCCEL1=3)
WP_NAMES = {1: "Cell prolif/EMT", 2: "EMT", 3: "Metabolism"}


def line_of(s):
    return re.sub(r"_r\d+$", "", s)


def paper_label(d, samples):
    ref = {"WP": PAPER_WP, "WPq": PAPER_WP, "PP": PAPER_PP, "YP": PAPER_YP}[d]
    lab = [ref[line_of(s)] for s in samples]
    if d == "YP":  # Fig. 4C: Ocum1_r1 in YP cluster 1
        lab = [1 if s == "Ocum1_r1" else l for s, l in zip(samples, lab)]
    return np.array(lab)


def load(name):
    d = pd.read_csv(os.path.join(PROC, name), sep="\t", keep_default_na=False, na_values=[""])
    s = [c for c in d.columns if re.match(r".+_r\d$", c)]
    m = d[s].astype(float); m.index = d["id"].astype(str)
    ann = d[[c for c in d.columns if c not in s]].set_index(m.index)
    return m, ann


def load_WP(): return load("WP_S3.tsv")
def load_PP(): return load("PP_S3.tsv")
def load_YP(): return load("YP_S3.tsv")
def load_WPq(): return load("WP_quantms_processed.tsv")


def ari(a, b):
    a = pd.Series(a).astype(str).values; b = pd.Series(b).astype(str).values
    tab = pd.crosstab(a, b).values
    from scipy.special import comb
    s1 = comb(tab, 2).sum(); s2 = comb(tab.sum(1), 2).sum(); s3 = comb(tab.sum(0), 2).sum()
    e = s2 * s3 / comb(tab.sum(), 2)
    return (s1 - e) / (0.5 * (s2 + s3) - e)


def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)
