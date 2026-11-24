import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_65c (S5): 分层置换零分布（预注册 30_ §二，修订 A2）。

每 unit×day 行：同日 NTC 池内无放回抽 n_cells 个伪 KD 集 × 200 次，
d_perm = mean(伪 KD) − mean(NTC 池)；checkpoint = T65/null/<cre_unit>__<day>.csv（单列 d_perm）。
幂等：已有 checkpoint 的行跳过（支持任意时刻断电续跑）。
"""
import os
import time
import numpy as np
import pandas as pd

BASE = _T65
NUL = os.path.join(BASE, "null")
N_PERM = 200
SEED = 20261106
LOG = os.path.join(_T65, "vc65c_run.log")


def say(msg):
    m = "[vc_65c %s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def main():
    os.makedirs(NUL, exist_ok=True)
    rng = np.random.default_rng(SEED)
    sc = pd.read_csv(os.path.join(BASE, "scores.csv.gz"))
    obs = pd.read_csv(os.path.join(BASE, "obs_displacement.csv"))
    ntc_L = {d: v.to_numpy(float) for d, v in sc[sc["cre_unit"] == "NTC"].groupby("day")["L"]}
    say("NTC 池规模：%s" % {k: len(v) for k, v in ntc_L.items()})

    todo = 0
    t0 = time.time()
    for _, r in obs.iterrows():
        tag = "%s__%s.csv" % (r["cre_unit"], r["day"])
        fp = os.path.join(NUL, tag)
        if os.path.exists(fp):
            continue
        pool = ntc_L[r["day"]]
        k = int(r["n_cells"])
        k = min(k, len(pool) - 1)
        draws = np.empty(N_PERM)
        for i in range(N_PERM):
            sel = rng.choice(len(pool), size=k, replace=False)
            draws[i] = pool[sel].mean() - pool.mean()
        pd.DataFrame({"d_perm": draws}).to_csv(fp, index=False)
        todo += 1
    say("本轮新算 %d 行，累计 %d/%d（%.1f min）" % (todo, len(os.listdir(NUL)), len(obs), (time.time() - t0) / 60))
    assert len(os.listdir(NUL)) == len(obs), "checkpoint 数与行数不等"


if __name__ == "__main__":
    main()
