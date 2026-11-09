"""Minimal shared helpers for the synthesis run (self-contained; no dependency on prior-run code)."""
from pathlib import Path

import numpy as np
import pandas as pd

RUN_DIR = Path(__file__).resolve().parent.parent
RES = RUN_DIR / "results"
FIG = RES / "figures"
TAB = RES / "tables"
for p in (FIG, TAB):
    p.mkdir(parents=True, exist_ok=True)

REPO = RUN_DIR.parent.parent.parent  # repo root (runs/0005/<id>/ -> repo)
RUN0001 = REPO / "runs/0001/20261105T121500Z_22794b7b"
RUN0002 = REPO / "runs/0002/20261107T120709Z_1d6e8ba6"
RUN0003 = REPO / "runs/0003/20261109T095213Z_f4b61cda"
RUN0004 = REPO / "runs/0004/20261109T132833Z_5b600c70"

# cell line -> paper integrative subtype (Hu et al. 2024)
SUBTYPE = {"AGS": "EMT", "eAGS": "EMT", "AZ521": "EMT", "SNU1750": "EMT", "MKN28": "EMT",
           "N87": "EMT", "TMK1": "EMT", "NCC24": "Metabolism", "SNU719": "Metabolism",
           "YCCEL1": "Metabolism", "KatoIII": "Metabolism", "Ocum1": "Metabolism",
           "SNU601": "Metabolism", "YCC3": "Metabolism"}
ORDER = ["AGS", "eAGS", "AZ521", "SNU1750", "MKN28", "N87", "TMK1",
         "NCC24", "SNU719", "YCCEL1", "KatoIII", "Ocum1", "SNU601", "YCC3"]
SUBTYPE_COLORS = {"EMT": "#eb6834", "Metabolism": "#2a78d6"}
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"


def diverging_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "bgr", ["#104281", "#3987e5", "#cde2fb", "#f0efec", "#f9d0c4", "#e34948", "#8f1d1d"])


def style():
    import matplotlib as mpl
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "figure.facecolor": "white",
        "axes.facecolor": "white", "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.titlelocation": "left"})


def savefig(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")


def read_gmt(path, universe=None, min_size=3, max_size=500):
    sets = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        g = set(p[2:])
        if universe is not None:
            g &= universe
        if min_size <= len(g) <= max_size:
            sets[p[0]] = g
    return sets
