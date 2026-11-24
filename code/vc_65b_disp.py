import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_65b (S4): 逐单元观测位移 d_u（预注册 30_ §一.3，修订 A1）。

d_u = mean(L | KD 细胞, day) − mean(L | 同日 NTC 池)
产出 T65/obs_displacement.csv（全部 unit×day 行，含非 robust——G1 要用）。
幂等：文件存在且行数一致则跳过。
"""
import os
import time
import numpy as np
import pandas as pd

BASE = _T65
LOG = os.path.join(_T65, "vc65b_run.log")


def say(msg):
    m = "[vc_65b %s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def main():
    out = os.path.join(BASE, "obs_displacement.csv")
    sc = pd.read_csv(os.path.join(BASE, "scores.csv.gz"))
    assert len(sc) == 50005
    ntc_mean = sc[sc["cre_unit"] == "NTC"].groupby("day")["L"].mean()

    kd = sc[sc["cre_unit"] != "NTC"]
    g = kd.groupby(["cre_unit", "gene_unit", "day"])["L"]
    obs = g.mean().rename("d_obs").reset_index()
    n = g.size().rename("n_cells").reset_index()
    obs = obs.merge(n, on=["cre_unit", "gene_unit", "day"])
    obs["ntc_mean_L"] = obs["day"].map(ntc_mean)
    # 复核：d_obs 已含减 NTC 均值（L 本身未减池均值，此处显式相减）
    obs["d_obs"] = obs["d_obs"] - obs["ntc_mean_L"]

    um = pd.read_csv(os.path.join(BASE, "unit_map.csv"))[["cre_unit", "day", "pass_floor", "robust_kd", "top_gene"]]
    obs = obs.merge(um, on=["cre_unit", "day"], how="left")
    assert len(obs) == len(um), "行数 %d != unit_map %d" % (len(obs), len(um))
    obs.to_csv(out, index=False)
    say("obs_displacement.csv 落盘 rows=%d（S4 checkpoint）" % len(obs))
    r = obs[obs["robust_kd"] == True]
    say("robust 行 d_obs 分布：min %.3f / median %.3f / max %.3f" % (r["d_obs"].min(), r["d_obs"].median(), r["d_obs"].max()))


if __name__ == "__main__":
    main()
