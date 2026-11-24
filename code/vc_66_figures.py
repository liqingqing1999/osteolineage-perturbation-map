import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_66_figures.py - manuscript figures F2/F4/S2. Reads deposited results only; idempotent.

F2  扰动响应地图：(a) 全 robust 基因级火山（D_g vs 置换 p）＋ (b) top-15 森林（带 null 95% CI）
F4  阴性校准：(a) 真实集 vs NTC 伪单元违反率 × n_cells 四分位（Wilson 95% CI）＋ (b) 总体率对照
S2  G2 先验基因森林（升级版：CI＋一致/反例着色）
产出：T65/figs/{F2_map,F4_negative_calibration,S2_g2_forest}.png ＋ T65/gene_level_all.csv
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = _T65
NUL = os.path.join(BASE, "null")
FIG = os.path.join(BASE, "figs")
os.makedirs(FIG, exist_ok=True)

PRIOR = {  # 基因 -> (类别, 预期 D_g 符号)
    "DLX5": ("pro-osteogenic", +1), "SMAD5": ("pro-osteogenic", +1),
    "RSPO3": ("pro-osteogenic", +1), "HOXC8": ("pro-osteogenic", +1),
    "HOXC9": ("pro-osteogenic", +1), "TCF7L1": ("Wnt-repressor (inverted)", -1),
    "RARB": ("RA / anti-adipogenic", +1), "NR3C1": ("pro-adipogenic", -1),
    "HMGA2": ("pro-adipogenic", -1),
}
C_OS = "#2166AC"   # 负位移 = 向成骨
C_AD = "#B2182B"   # 正位移 = 向成脂
C_GY = "#9E9E9E"


def load_null(r):
    fp = os.path.join(NUL, "%s__%s.csv" % (r["cre_unit"], r["day"]))
    return pd.read_csv(fp)["d_perm"].to_numpy(float)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(c - h, 0), min(c + h, 1))


def gene_table(obs):
    """全 robust 基因级表：与 vc_65d G2 完全同口径（robust ∧ pass_floor，n_cells 加权）。"""
    sub_all = obs[(obs["robust_kd"] == True) & (obs["pass_floor"] == True)]
    rows = []
    for g, sub in sub_all.groupby("top_gene"):
        if not g or pd.isna(g):
            continue
        w = sub["n_cells"].to_numpy(float)
        d = sub["d_obs"].to_numpy(float)
        W = w.sum()
        dg = float(np.dot(w, d) / W)
        dperm_agg = np.zeros(200)
        for _, r in sub.iterrows():
            dperm_agg += load_null(r) * r["n_cells"] / W
        p_two = (1 + int((np.abs(dperm_agg) >= abs(dg)).sum())) / 201.0
        lo, hi = np.percentile(dperm_agg, [2.5, 97.5])
        cat, exp = PRIOR.get(g, ("none", 0))
        rows.append(dict(gene=g, prior=cat, expected_sign=exp, D_g=dg,
                         ci_lo=float(lo), ci_hi=float(hi),
                         n_units=len(sub), n_cells=int(W), p_perm_two_sided=p_two,
                         concordant=(int(np.sign(dg) == exp) if exp != 0 else np.nan)))
    return pd.DataFrame(rows)


def fig_F2(gt):
    gt = gt.sort_values("p_perm_two_sided").reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), dpi=300)

    # (a) 火山
    ax = axes[0]
    is_prior = gt["prior"] != "none"
    x = gt["D_g"].to_numpy(float)
    y = -np.log10(gt["p_perm_two_sided"].to_numpy(float))
    sizes = 8 + 30 * np.sqrt(gt["n_cells"] / gt["n_cells"].max())
    ax.scatter(x[~is_prior], y[~is_prior], s=sizes[~is_prior], c=C_GY, alpha=.55,
               lw=.3, edgecolor="white", label="eBMD / other targets (n=%d)" % (~is_prior).sum())
    for cat, col, mk in [("pro-osteogenic", C_OS, "o"), ("pro-adipogenic", C_AD, "s"),
                         ("RA / anti-adipogenic", "#7B52AB", "^"), ("Wnt-repressor (inverted)", "#00897B", "D")]:
        m = gt["prior"] == cat
        if m.any():
            ax.scatter(x[m], y[m], s=sizes[m] * 1.5, c=col, marker=mk, lw=.4,
                       edgecolor="white", label=cat, zorder=3)
    for g in ["RSPO3", "NR3C1", "RARB", "HOXC8"]:
        r = gt[gt["gene"] == g]
        if len(r):
            ax.annotate(g, (r["D_g"].iloc[0], -np.log10(r["p_perm_two_sided"].iloc[0])),
                        textcoords="offset points", xytext=(6, 4), fontsize=6.5)
    ax.axvline(0, color="k", lw=.7, ls="--")
    ax.set_xlabel("$D_g$  (weighted displacement, $L=z_{AD}-z_{OS}$)")
    ax.set_ylabel("$-$log$_{10}$ permutation $p$ (two-sided)")
    ax.set_title("Perturbation-response map: %d genes" % len(gt), fontsize=9)
    ax.legend(fontsize=5.8, frameon=False, loc="upper left", handletextpad=.2)

    # (b) top-15 森林
    ax = axes[1]
    top = gt.reindex(gt["D_g"].abs().sort_values(ascending=False).index).head(15)
    top = top.iloc[::-1]
    ypos = np.arange(len(top))
    cols = [C_AD if v > 0 else C_OS for v in top["D_g"]]
    err = np.vstack([top["D_g"] - top["ci_lo"], top["ci_hi"] - top["D_g"]]).clip(min=0)
    ax.errorbar(top["D_g"], ypos, xerr=err, fmt="none", ecolor="#616161",
                elinewidth=.9, capsize=1.6, zorder=2)
    ax.scatter(top["D_g"], ypos, c=cols, s=26, zorder=3, lw=.4, edgecolor="white")
    ylab = ["%s \u2605" % r["gene"] if r["prior"] != "none" else r["gene"]
            for _, r in top.iterrows()]
    ax.set_yticks(ypos, ylab, fontsize=6.8)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.axvline(0, color="k", lw=.7, ls="--")
    ax.set_xlabel("$D_g$ with null 95% CI   (\u2605 pre-specified prior)")
    ax.set_xlim(-0.42, 0.75)
    ax.set_title("Top 15 by $|D_g|$", fontsize=9)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    fp = os.path.join(FIG, "F2_map.png")
    fig.savefig(fp, bbox_inches="tight")
    plt.close(fig)
    return fp


def fig_F4():
    sc = pd.read_csv(os.path.join(BASE, "scores.csv.gz"))
    obs = pd.read_csv(os.path.join(BASE, "obs_displacement.csv"))
    g1 = obs[(obs["robust_kd"] == False) & (obs["pass_floor"] == True)
             & (obs["cre_unit"].str.startswith("SNP-"))].copy()
    g1["viol"] = [abs(r["d_obs"]) > np.percentile(np.abs(load_null(r)), 95)
                  for _, r in g1.iterrows()]
    ps = pd.read_csv(os.path.join(BASE, "g1v2_pseudo.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)

    # (a) 四分位对比
    ax = axes[0]
    qr = pd.qcut(g1["n_cells"], 4, duplicates="drop")
    qp = pd.qcut(ps["n_cells"], 4, duplicates="drop")
    rq = g1.groupby(qr, observed=True)["viol"].agg(["mean", "sum", "size"])
    pq = ps.groupby(qp, observed=True)["viol"].agg(["mean", "sum", "size"])
    labels = ["%s\n%d\u2013%d cells" % (q, int(iv.left), int(iv.right))
              for q, iv in zip(["Q1", "Q2", "Q3", "Q4"], rq.index)]
    xpos = np.arange(len(rq))
    wr = [wilson(int(r["sum"]), int(r["size"])) for _, r in rq.iterrows()]
    wp = [wilson(int(r["sum"]), int(r["size"])) for _, r in pq.iterrows()]
    rq_lo = np.array([m - lo for m, (lo, hi) in zip(rq["mean"], wr)])
    rq_hi = np.array([hi - m for m, (lo, hi) in zip(rq["mean"], wr)])
    pq_lo = np.array([m - lo for m, (lo, hi) in zip(pq["mean"], wp)])
    pq_hi = np.array([hi - m for m, (lo, hi) in zip(pq["mean"], wp)])
    ax.errorbar(xpos - .12, rq["mean"].to_numpy(), yerr=np.vstack([rq_lo, rq_hi]),
                fmt="s", color=C_AD, ms=5, capsize=2.4, lw=1.1,
                label="non-robust GWAS-CRE units (real)")
    ax.errorbar(xpos + .12, pq["mean"].to_numpy(), yerr=np.vstack([pq_lo, pq_hi]),
                fmt="o", color="#455A64", ms=5, capsize=2.4, lw=1.1,
                label="NTC pseudo-units (noise floor)")
    ax.axhline(.05, color="#9E9E9E", lw=.8, ls=":", label="nominal 5%")
    ax.axhline(1.5 * float(ps["viol"].mean()), color=C_AD, lw=.8, ls="--")
    ax.text(len(rq) - .55, 1.5 * float(ps["viol"].mean()) + .004,
            "G1-v2 threshold (1.5\u00d7 floor)", fontsize=5.8, color=C_AD, ha="right")
    ax.set_xticks(xpos, labels, fontsize=6.5)
    ax.set_ylabel("violating own null p95")
    ax.set_ylim(0, .22)
    ax.set_title("Violation rate vs. unit power (quartiles)", fontsize=9)
    ax.legend(fontsize=6, frameon=False, loc="upper left")

    # (b) 总体率
    ax = axes[1]
    rr, pp = float(g1["viol"].mean()), float(ps["viol"].mean())
    vals, cs = [rr, pp], [C_AD, "#455A64"]
    ypos = [0, 1]
    for yi, v, c in zip(ypos, vals, cs):
        ax.barh(yi, v, height=.52, color=c)
        ax.text(v + .004, yi, "%.1f%% (%d/%d)" % (100 * v, round(v * len(g1 if yi == 0 else ps)),
                                                  len(g1 if yi == 0 else ps)),
                va="center", fontsize=7)
    ax.axvline(1.5 * pp, color=C_AD, lw=.9, ls="--")
    ax.axvline(.05, color="#9E9E9E", lw=.8, ls=":")
    ax.text(.051, -0.42, "nominal 5%", fontsize=6, color="#757575")
    ax.set_yticks(ypos, ["real\nnegative set", "NTC\npseudo floor"], fontsize=7.5)
    ax.set_xlim(0, .16)
    ax.set_ylim(-0.55, 1.6)
    ax.set_xlabel("overall violation rate")
    ax.set_title("G1-v2: FAIL (%.1f%% > %.1f%% = 1.5\u00d7floor);\nflat floor \u21d2 readout unbiased" % (100 * rr, 100 * 1.5 * pp), fontsize=8.5)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    fp = os.path.join(FIG, "F4_negative_calibration.png")
    fig.savefig(fp, bbox_inches="tight")
    plt.close(fig)
    return fp


def fig_S2(gl):
    gl = gl.sort_values("D_g")
    fig, ax = plt.subplots(figsize=(4.6, 3.2), dpi=300)
    ypos = np.arange(len(gl))
    cols = ["#2E7D32" if c == 1 else "#C62828" for c in gl["concordant"]]
    err = np.vstack([gl["D_g"] - gl["ci_lo"], gl["ci_hi"] - gl["D_g"]]).clip(min=0)
    ax.errorbar(gl["D_g"], ypos, xerr=err, fmt="none", ecolor="#616161",
                elinewidth=.9, capsize=1.8, zorder=2)
    ax.scatter(gl["D_g"], ypos, c=cols, s=30, zorder=3, lw=.4, edgecolor="white")
    ax.set_yticks(ypos, ["%s\n(%s)" % (r["gene"], r["prior"]) for _, r in gl.iterrows()],
                  fontsize=6.2)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.axvline(0, color="k", lw=.7, ls="--")
    ax.set_xlim(-.5, .8)
    ax.set_xlabel("$D_g$ with null 95% CI   (green = concordant, red = discordant)")
    ax.set_title("G2: prior-anchored direction concordance\n8/9, one-sided sign-test p = .0195", fontsize=9)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fp = os.path.join(FIG, "S2_g2_forest.png")
    fig.savefig(fp, bbox_inches="tight")
    plt.close(fig)
    return fp


def main():
    obs = pd.read_csv(os.path.join(BASE, "obs_displacement.csv"))
    gt = gene_table(obs)
    gt.to_csv(os.path.join(BASE, "gene_level_all.csv"), index=False)
    prior_gl = gt[gt["prior"] != "none"].copy()
    # 与 gene_level.csv（vc_65d 产物）对账
    ref = pd.read_csv(os.path.join(BASE, "gene_level.csv"))
    mg = ref.merge(prior_gl, on="gene", suffixes=("_ref", ""))
    assert len(mg) == 9
    assert np.allclose(mg["D_g_ref"], mg["D_g"], atol=1e-12), "D_g 与 vc_65d 不一致"
    assert np.allclose(mg["p_perm_two_sided_ref"], mg["p_perm_two_sided"], atol=1e-12)
    k = int((mg["expected_sign_ref"] == np.sign(mg["D_g"])).sum())
    print("对账 OK：9 先验基因 D_g/p 逐位一致，concordant=%d/9" % k)

    f2 = fig_F2(gt)
    f4 = fig_F4()
    s2 = fig_S2(prior_gl.sort_values("D_g"))
    print("robust genes mapped: %d" % len(gt))
    for f in (f2, f4, s2):
        print("written:", f)


if __name__ == "__main__":
    main()
