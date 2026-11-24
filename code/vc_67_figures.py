import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_67_figures.py - manuscript figures F1/F3 (Physiological Genomics version).

红线：
- 所有数字从落盘产物读出并断言，不手填结果性数值；
- vc_60 标定数字先对 _w60_calib_report.md 做存在断言再画；
- 产物落 data/derived/GSE319266/T65/figs/。
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

BASE = _T65
W1C = os.path.join(_T65, "reference")
REP = os.path.join(_T65, "reference", "calibration_report.md")
OUT = os.path.join(BASE, "figs")
os.makedirs(OUT, exist_ok=True)

C_AD = "#C62828"   # red  = adipogenic pole
C_OS = "#1565C0"   # blue = osteogenic pole
GREY = "#455A64"

def say(m):
    print("[vc_67]", m, flush=True)

def load_w1c():
    p = os.path.join(W1C, "w1c_verdict.csv")
    rows = list(csv.DictReader(open(p, encoding="utf-8")))
    assert len(rows) == 12, len(rows)
    assert all(r["verdict"] == "PASS" for r in rows), "W1c 12/12 须全 PASS"
    return rows

def fig_f1(rows):
    """F1: (a) 坐标轴示意 (b) D1 臂 LODO/CV (c) 全 12 项判定。"""
    fig = plt.figure(figsize=(7.2, 4.6), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], hspace=0.55, wspace=0.28)

    # ---- (a) schematic ----
    ax = fig.add_subplot(gs[0, :])
    ax.set_xlim(0, 10); ax.set_ylim(0, 1.75); ax.axis("off")
    rng = np.random.default_rng(7)
    xs = 0.5 * rng.standard_normal(220) + 2.6
    xs2 = 0.5 * rng.standard_normal(220) + 7.4
    ax.plot(xs, 1.30 + 0.10 * rng.standard_normal(220), ".", ms=1.0, color=C_AD, alpha=.5)
    ax.plot(xs2, 1.30 + 0.10 * rng.standard_normal(220), ".", ms=1.0, color=C_OS, alpha=.5)
    ax.text(2.6, 1.68, "AD program (11 genes)", ha="center", fontsize=6.8, color=C_AD)
    ax.text(7.4, 1.68, "OS program (12 genes)", ha="center", fontsize=6.8, color=C_OS)
    ax.add_patch(FancyArrowPatch((0.8, 0.95), (9.2, 0.95), arrowstyle="-|>",
                                 mutation_scale=12, lw=1.4, color="k"))
    ax.text(0.8, 0.62, "adipogenic pole", ha="left", fontsize=7.5, color=C_AD)
    ax.text(9.2, 0.62, "osteogenic pole", ha="right", fontsize=7.5, color=C_OS)
    ax.text(5.0, 0.66, "$L = z(\\mathrm{AD}) - z(\\mathrm{OS})$", ha="center",
            va="top", fontsize=8)
    ax.text(5.0, 0.12, "CRISPRi perturbation $\\Rightarrow$ shift of the cell-state distribution $d$ along $L$",
            ha="center", fontsize=6.8, color="k")
    ax.set_title("(a) Lineage axis and readout", fontsize=8.5, loc="left")

    # ---- (b) D1 arm LODO vs CV ----
    ax = fig.add_subplot(gs[1, 0])
    arm = [r for r in rows if r["arm"] == "D1"]
    labels, lodo, cv, cols = [], [], [], []
    for r in arm:
        labels.append("%s\n%s" % (r["resp"].upper(), r["ytrack"]))
        lodo.append(float(r["lodo_med"])); cv.append(float(r["cv_med"]))
        cols.append(C_AD if r["resp"] == "AD" else C_OS)
    x = np.arange(len(labels))
    ax.bar(x - .19, lodo, width=.36, color=cols, alpha=.9, label="LODO median")
    ax.bar(x + .19, cv, width=.36, color=cols, alpha=.38, label="cross-validated")
    for xi, v in zip(x - .19, lodo):
        ax.text(xi, v + .025, "%.2f" % v, ha="center", fontsize=5.6, zorder=4)
    ax.axhline(.5, color="grey", lw=.7, ls=":", zorder=1)
    ax.text(3.62, .415, "chance", fontsize=5.6, color="grey", va="top")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=5.6)
    ax.set_ylim(0, 1.08); ax.set_ylabel("classification accuracy")
    ax.tick_params(axis="y", labelsize=6.5)
    ax.legend(fontsize=5.8, loc="upper left", frameon=False)
    ax.set_xlabel("response (marker: UCell / centroid)", fontsize=7)
    ax.set_title("(b) D1 external validation (W1c, 12/12 PASS)", fontsize=8.5, loc="left")

    # ---- (c) all arms ----
    ax = fig.add_subplot(gs[1, 1])
    for i, r in enumerate(rows):
        col = C_AD if r["resp"] == "AD" else C_OS
        ax.scatter(float(r["lodo_med"]), i, s=16, color=col, zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(["%s|%s|%s" % (r["mver"], r["arm"], r["resp"].upper()) for r in rows],
                       fontsize=5.4)
    ax.invert_yaxis()
    ax.axvline(.5, color="grey", lw=.7, ls=":")
    ax.set_xlim(0, 1.0)
    ax.tick_params(axis="x", labelsize=6.5)
    ax.set_xlabel("LODO median accuracy")
    ax.set_title("(c) All 12 donor-validation cells", fontsize=8.5, loc="left")

    p = os.path.join(OUT, "F1_coordinate_validation.png")
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    say("F1 -> %s" % p)

def fig_f3():
    """F3: 检出限 vs 设计，DEU 锚点带。数字先对报告断言。"""
    rep = open(REP, encoding="utf-8").read()
    for s in ["0.049 / 0.079 / 0.107", "**0.029**", "0.021", "0.034",
              "0.016–0.025", "σ̂=0.197", "-0.0228", "-0.0251", "-0.0160"]:
        assert s in rep, "报告缺期望串: %s" % s
    configs = [
        ("single donor\n150/group", [0.049, 0.079, 0.107], "#90A4AE"),
        ("3 donors\n150x3/group", [0.034], "#78909C"),
        ("4 donors\n150x4/group", [0.029], "#1565C0"),
        ("4 donors\n300x4/group", [0.021], "#0D47A1"),
    ]
    fig, ax = plt.subplots(figsize=(5.4, 3.2), dpi=300)
    xs, ls = [], []
    x = 0
    for name, vals, col in configs:
        for v in vals:
            ax.bar(x, v, width=.62, color=col, alpha=.85)
            ax.text(x, v + .002, "%.3f" % v, ha="center", fontsize=6.2)
            xs.append(x); ls.append(name); x += 1
        x += .55
    # natural-drift band 1 DEU = 0.016–0.025
    ax.axhspan(.016, .025, color="#C62828", alpha=.14, zorder=0)
    ax.text(x - .4, .0162, "natural disease drift: 1 DEU = |$\\Delta$c| 0.016–0.025",
            fontsize=6.4, color=C_AD, ha="right", va="bottom")
    ax.axhline(.0083, color="grey", lw=.6, ls=":")
    ax.set_xticks(xs)
    ax.set_xticklabels([l.replace("\n", " ") for l in ls], fontsize=5.6, rotation=20,
                       ha="right")
    ax.set_ylabel("detection limit  |$\\Delta$c|\n($\\alpha$=.05/30, power=.8, $\\hat\\sigma$=0.197)",
                  fontsize=7.5)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylim(0, .13)
    ax.set_title("Calibration of perturbation effects in natural-drift units (DEU)",
                 fontsize=8.5, loc="left")
    p = os.path.join(OUT, "F3_deu_calibration.png")
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    say("F3 -> %s" % p)

if __name__ == "__main__":
    rows = load_w1c()
    fig_f1(rows)
    fig_f3()
    say("done")
