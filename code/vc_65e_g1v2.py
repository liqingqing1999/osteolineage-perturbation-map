import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_65e: G1-v2 —— NTC 自校准楼层（修订附注 A3，2026-11-06 锁定后执行）。

对 796 个真实阴性单元各生成 1 个伪单元（同 day NTC 池内无放回抽 n_cells；
d_pseudo = mean(sel) − mean(rest)），零分布 = 同池 200 次同构抽取。
判定：r_real ≤ 1.5 × r_pseudo ⇒ G1-v2 PASS。
产出：T65/g1v2_pseudo.csv、T65/g1v2_verdict.tsv；追加 _vc65e_log.txt。
"""
import os
import time
import numpy as np
import pandas as pd

BASE = _T65
N_PERM = 200
SEED = 20261107
LOG = os.path.join(_T65, "vc65e_run.log")


def say(msg):
    m = "[vc_65e %s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def main():
    rng = np.random.default_rng(SEED)
    sc = pd.read_csv(os.path.join(BASE, "scores.csv.gz"))
    obs = pd.read_csv(os.path.join(BASE, "obs_displacement.csv"))
    g1 = obs[(obs["robust_kd"] == False) & (obs["pass_floor"] == True) & (obs["cre_unit"].str.startswith("SNP-"))]
    # 真实集违反率（vc_65d 已算 10.3%，此处重算保证同源）
    def real_viol(r):
        dp = pd.read_csv(os.path.join(BASE, "null", "%s__%s.csv" % (r["cre_unit"], r["day"])))["d_perm"].to_numpy(float)
        return abs(r["d_obs"]) > np.percentile(np.abs(dp), 95)
    g1 = g1.copy()
    g1["viol"] = [real_viol(r) for _, r in g1.iterrows()]
    r_real = float(g1["viol"].mean())
    say("真实阴性集 r_real=%.4f（%d/%d）" % (r_real, int(g1["viol"].sum()), len(g1)))

    ntc_L = {d: v.to_numpy(float) for d, v in sc[sc["cre_unit"] == "NTC"].groupby("day")["L"]}
    out = []
    for _, r in g1.iterrows():
        pool = ntc_L[r["day"]]
        k = min(int(r["n_cells"]), len(pool) // 2)   # 保证 rest 非空且有足够对照性
        sel = rng.choice(len(pool), size=k, replace=False)
        mask = np.zeros(len(pool), bool); mask[sel] = True
        d_ps = pool[mask].mean() - pool[~mask].mean()
        dp = np.empty(N_PERM)
        for i in range(N_PERM):
            s2 = rng.choice(len(pool), size=k, replace=False)
            m2 = np.zeros(len(pool), bool); m2[s2] = True
            dp[i] = pool[m2].mean() - pool[~m2].mean()
        out.append(dict(cre_unit=r["cre_unit"], day=r["day"], n_cells=k,
                        d_pseudo=d_ps, viol=abs(d_ps) > np.percentile(np.abs(dp), 95)))
    ps = pd.DataFrame(out)
    ps.to_csv(os.path.join(BASE, "g1v2_pseudo.csv"), index=False)
    r_ps = float(ps["viol"].mean())
    q = pd.qcut(ps["n_cells"], 4, duplicates="drop")
    tab_q = ps.groupby(q, observed=True)["viol"].agg(["mean", "size"]).round(4)
    qq = pd.qcut(g1["n_cells"], 4, duplicates="drop")
    tab_rq = g1.groupby(qq, observed=True)["viol"].agg(["mean", "size"]).round(4)
    say("伪单元楼层 r_pseudo=%.4f（%d/%d）" % (r_ps, int(ps["viol"].sum()), len(ps)))
    say("伪单元×四分位：\n%s" % tab_q)
    say("真实集×四分位：\n%s" % tab_rq)
    thr = 1.5 * r_ps
    verdict = "PASS" if r_real <= thr else "FAIL"
    say("G1-v2 判定：r_real=%.4f ≤ 1.5×r_pseudo=%.4f ? ⇒ %s" % (r_real, thr, verdict))

    lines = ["G1-v2 (NTC self-calibrated floor)",
             "r_real=%.4f  r_pseudo=%.4f  threshold=1.5x=%.4f" % (r_real, r_ps, thr),
             "verdict=%s" % verdict]
    with open(os.path.join(BASE, "g1v2_verdict.tsv"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    # A3 重映射：A3 = PASS ⇔ Q1 PASS ∧ G1-v2 PASS ∧ G2 PASS
    q1_pass, g2_pass = True, True   # 31_ 判定表：Q1 p=.015 PASS；G2 8/9 p=.0195 PASS
    a3 = "A3 PASS" if (q1_pass and verdict == "PASS" and g2_pass) else "A3 FAIL"
    say("== %s ==" % a3)
    with open(os.path.join(BASE, "A3_verdict.txt"), "w") as f:
        f.write(a3 + "\n")


if __name__ == "__main__":
    main()
