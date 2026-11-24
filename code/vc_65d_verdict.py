import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_65d (S6): 判定 Q1/G1/G2/G3（预注册 30_ §三，含修订 A1/A2）。

产出：T65/verdict.tsv、T65/gene_level.csv、T65/vc65_forest.png、T65/vc65_g1floor.png
幂等：全量重算（输入只读，秒级）。
"""
import os
import time
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = _T65
NUL = os.path.join(BASE, "null")
LOG = os.path.join(_T65, "vc65d_run.log")

PRIOR = {  # 基因 -> (先验类别, 预期 D_g 符号)
    "DLX5": ("pro-osteogenic", +1), "SMAD5": ("pro-osteogenic", +1),
    "RSPO3": ("pro-osteogenic", +1), "HOXC8": ("pro-osteogenic", +1),
    "HOXC9": ("pro-osteogenic", +1), "TCF7L1": ("Wnt-repressor(inverted)", -1),
    "RARB": ("RA/anti-adipogenic", +1), "NR3C1": ("pro-adipogenic", -1),
    "HMGA2": ("pro-adipogenic", -1),
}


def say(msg):
    m = "[vc_65d %s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def load_null(r):
    fp = os.path.join(NUL, "%s__%s.csv" % (r["cre_unit"], r["day"]))
    return pd.read_csv(fp)["d_perm"].to_numpy(float)


def main():
    sc = pd.read_csv(os.path.join(BASE, "scores.csv.gz"))
    obs = pd.read_csv(os.path.join(BASE, "obs_displacement.csv"))
    rows = []

    # ---- Q1：CD55-MOI 敏感性 QC（逐 day 2D 占据位移，标签置换 200 次）----
    moi = sc[sc["cre_unit"] == "CD55-MOI"].reset_index(drop=True)
    ntc = sc[sc["cre_unit"] == "NTC"].reset_index(drop=True)
    rng = np.random.default_rng(20261106)
    dists_obs, dists_perm = [], []
    for d in sorted(moi["day"].unique()):
        A = moi[moi["day"] == d][["AD_UCell", "OS_UCell"]].to_numpy(float)
        B = ntc[ntc["day"] == d][["AD_UCell", "OS_UCell"]].to_numpy(float)
        P = np.vstack([A, B])

        def occ_dist(idx_a, data):
            oa = np.histogram2d(data[idx_a, 0], data[idx_a, 1], bins=10,
                                range=[[P[:, 0].min(), P[:, 0].max()], [P[:, 1].min(), P[:, 1].max()]])[0]
            ob = np.histogram2d(data[~idx_a, 0], data[~idx_a, 1], bins=10,
                                range=[[P[:, 0].min(), P[:, 0].max()], [P[:, 1].min(), P[:, 1].max()]])[0]
            oa = oa / max(oa.sum(), 1); ob = ob / max(ob.sum(), 1)
            oa, ob = oa.ravel() + 1e-9, ob.ravel() + 1e-9
            return 1 - float(np.dot(oa, ob) / (np.linalg.norm(oa) * np.linalg.norm(ob)))

        ia = np.concatenate([np.ones(len(A), bool), np.zeros(len(B), bool)])
        d0 = occ_dist(ia, P)
        dp = np.array([occ_dist(rng.permutation(ia), P) for _ in range(200)])
        w = len(A)
        dists_obs.append(d0 * w); dists_perm.append(dp * w)
    d_obs_q1 = float(np.sum(dists_obs) / len(moi))
    dp_q1 = np.concatenate(dists_perm) / len(moi)
    p_q1 = (1 + int((dp_q1 >= d_obs_q1).sum())) / 201.0
    rows.append(("Q1", "CD55-MOI 2D occupancy shift (within-day, weighted)", "dist=%.4f null95=%.4f p=%.3f" % (d_obs_q1, np.percentile(dp_q1, 95), p_q1), "PASS" if p_q1 < 0.05 else "FAIL"))

    # ---- G1：阴性楼层（无 robust 标签的 SNP-CRE 单元，过功效线）----
    g1 = obs[(obs["robust_kd"] == False) & (obs["pass_floor"] == True) & (obs["cre_unit"].str.startswith("SNP-"))]
    viol = 0
    for _, r in g1.iterrows():
        dperm = load_null(r)
        if abs(r["d_obs"]) > np.percentile(np.abs(dperm), 95):
            viol += 1
    frac_ok = 1 - viol / len(g1)
    rows.append(("G1", "negative floor: non-robust SNP units within own null p95", "within=%.3f (%d/%d)" % (frac_ok, len(g1) - viol, len(g1)), "PASS" if frac_ok >= 0.90 else "FAIL"))

    # ---- G2：先验锚定方向一致性 ----
    gene_rows = []
    for g, (cat, exp) in PRIOR.items():
        sub = obs[(obs["top_gene"] == g) & (obs["robust_kd"] == True) & (obs["pass_floor"] == True)]
        w = sub["n_cells"].to_numpy(float); d = sub["d_obs"].to_numpy(float)
        dg = float(np.dot(w, d) / w.sum())
        dperm_agg = np.zeros(200)
        for _, r in sub.iterrows():
            dperm_agg += load_null(r) * r["n_cells"] / w.sum()
        p_two = (1 + int((np.abs(dperm_agg) >= abs(dg)).sum())) / 201.0
        conc = int(np.sign(dg) == exp)
        gene_rows.append(dict(gene=g, prior=cat, expected_sign=exp, D_g=dg,
                              n_units=len(sub), n_cells=int(w.sum()), concordant=conc,
                              p_perm_two_sided=p_two))
    gl = pd.DataFrame(gene_rows)
    gl.to_csv(os.path.join(BASE, "gene_level.csv"), index=False)
    k = int(gl["concordant"].sum()); n = len(gl)
    p_g2 = float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue)
    rows.append(("G2", "prior-anchored direction concordance (sign test)", "%d/%d concordant, one-sided p=%.4f" % (k, n, p_g2), "PASS" if (n >= 8 and p_g2 < 0.05) else "FAIL"))

    # ---- G3：留一 batch（Full-1..6）＋ day2/4 敏感性 ----
    ntc_mean_all = sc[sc["cre_unit"] == "NTC"].groupby("day")["L"].mean()
    folds_ok = 0
    for run in ["Full-1", "Full-2", "Full-3", "Full-4", "Full-5", "Full-6"]:
        kk = 0; nn = 0
        for g, (cat, exp) in PRIOR.items():
            sub = sc[(sc["gene_unit"] == g) & (sc["run"] != run) & (sc["cre_unit"] != "NTC")]
            sub = sub.merge(ntc_mean_all.rename("nm"), left_on="day", right_index=True).reset_index(drop=True)
            agg = (sub.assign(Ld=sub["L"] - sub["nm"]).groupby(["cre_unit", "day"])
                      .agg(n=("Ld", "size"), m=("Ld", "mean")).reset_index())
            agg = agg[agg["n"] >= 10]
            if len(agg) == 0:
                continue
            dg = np.average(agg["m"], weights=agg["n"])
            nn += 1; kk += int(np.sign(dg) == exp)
        if nn >= 8 and stats.binomtest(kk, nn, 0.5, alternative="greater").pvalue < 0.05:
            folds_ok += 1
    # day2+4 敏感性
    kk4 = 0; nn4 = 0
    for g, (cat, exp) in PRIOR.items():
        sub = obs[(obs["top_gene"] == g) & (obs["robust_kd"] == True) & (obs["pass_floor"] == True) & (obs["day"].isin(["day2", "day4"]))]
        if len(sub) == 0:
            continue
        dg = np.average(sub["d_obs"], weights=sub["n_cells"])
        nn4 += 1; kk4 += int(np.sign(dg) == exp)
    p_d24 = float(stats.binomtest(kk4, nn4, 0.5, alternative="greater").pvalue) if nn4 >= 8 else np.nan
    rows.append(("G3", "LOBO stability (6 Full runs) + day2/4 sensitivity", "LOBO significant folds=%d/6; day2/4: %d/%d p=%.4f" % (folds_ok, kk4, nn4, p_d24), "REPORT"))

    vd = pd.DataFrame(rows, columns=["gate", "definition", "result", "verdict"])
    vd.to_csv(os.path.join(BASE, "verdict.tsv"), sep="\t", index=False)
    say("\n" + vd.to_string(index=False))
    ok = dict((r["gate"], r["verdict"]) for _, r in vd.iterrows())
    if ok["Q1"] == "PASS" and ok["G1"] == "PASS":
        a3 = "A3 PASS" if ok["G2"] == "PASS" else "A3 INCONCLUSIVE"
    else:
        a3 = "A3 FAIL"
    say("== %s ==" % a3)
    with open(os.path.join(BASE, "A3_verdict.txt"), "w") as f:
        f.write(a3 + "\n")

    # ---- 图：先验基因 forest ----------------------------------------
    gl2 = gl.sort_values("D_g")
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)
    ypos = np.arange(len(gl2))
    colors = ["#B2182B" if c else "#4393C3" for c in gl2["concordant"]]
    ax.errorbar(gl2["D_g"], ypos, xerr=1.96 * gl2["p_perm_two_sided"].clip(0.001) * 0 + gl2["D_g"] * 0,
                fmt="none", ecolor="gray")
    ax.scatter(gl2["D_g"], ypos, c=colors, s=42, zorder=3)
    for yi, (_, r) in zip(ypos, gl2.iterrows()):
        ax.text(r["D_g"], yi + 0.28, "p=%.3f" % r["p_perm_two_sided"], ha="center", fontsize=6.5, color="gray")
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_yticks(ypos, ["%s (%s)" % (r["gene"], r["prior"]) for _, r in gl2.iterrows()], fontsize=8)
    ax.set_xlabel("Displacement along lineage axis L = z(AD) - z(OS)  |  red = concordant with prior")
    ax.set_title("vc_65 G2: prior-anchored direction concordance (%d/%d, one-sided p=%.4f)" % (k, n, p_g2), fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(BASE, "vc65_forest.png"))
    # G1 楼层直方图
    fig2, ax2 = plt.subplots(figsize=(5.4, 3.8), dpi=200)
    z = []
    for _, r in g1.iterrows():
        dperm = load_null(r)
        z.append(abs(r["d_obs"]) / max(np.percentile(np.abs(dperm), 95), 1e-9))
    ax2.hist(np.clip(z, 0, 4), bins=40, color="#4d4d4d")
    ax2.axvline(1, color="#B2182B", lw=1.2, ls="--", label="null p95 (=1)")
    ax2.set_xlabel("|d_obs| / own null p95  (non-robust SNP units)"); ax2.set_ylabel("# units")
    ax2.legend(fontsize=7); ax2.set_title("G1 negative floor: %.1f%% within p95" % (100 * frac_ok), fontsize=9)
    plt.tight_layout(); fig2.savefig(os.path.join(BASE, "vc65_g1floor.png"))
    say("图与判定表落盘完成")


if __name__ == "__main__":
    main()
