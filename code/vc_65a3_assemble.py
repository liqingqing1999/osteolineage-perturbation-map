import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_65a3 (S2/S3 part 2): 组装打分表 + 单元功效盘点。

预注册 30_ §一/§二（修订 A1/A2）：
  读出 L = z(AD_UCell) − z(OS_UCell)，z 化用同日 NTC 池
产出：
  T65/scores.parquet     S2 checkpoint：cellkey, run, day, cre_unit, gene_unit, AD, OS, L
  T65/unit_map.csv       S3 checkpoint：(cre_unit, day) 功效盘点 + 纳入判定
幂等：两产物存在且行数断言过则跳过。
"""
import os
import time
import numpy as np
import pandas as pd

BASE = _T65
EXP = os.path.join(BASE, "export")
LOG = os.path.join(_T65, "vc65a3_run.log")

# KD 臂 ≥10 细胞（预注册 §一.6）；NTC 池按 day 全池 ≥30 断言式满足
KD_CELL_FLOOR = 10


def say(msg):
    m = "[vc_65a3 %s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def main():
    sp_path = os.path.join(BASE, "scores.csv.gz")
    um_path = os.path.join(BASE, "unit_map.csv")
    if os.path.exists(sp_path) and os.path.exists(um_path):
        n = len(pd.read_csv(sp_path))
        assert n == 50005
        say("两产物已存在（scores n=%d），幂等跳过" % n)
        return

    meta = pd.read_csv(os.path.join(EXP, "meta.csv"), keep_default_na=False)
    sc = pd.read_csv(os.path.join(EXP, "ucell_scores.csv"))
    assert len(meta) == 50005 and len(sc) == 50005
    sc = meta.merge(sc, on="cell", validate="1:1")
    assert len(sc) == 50005

    # L = z(AD) − z(OS)，z 用同日 NTC 池（修订 A1）
    ntc = sc[sc["cre_unit"] == "NTC"]
    stats = ntc.groupby("day")[["AD_UCell", "OS_UCell"]].agg(["mean", "std"])
    say("NTC 池统计：\n%s" % stats.round(4))

    def zrow(r):
        mAD, sAD = stats.loc[r["day"], ("AD_UCell", "mean")], stats.loc[r["day"], ("AD_UCell", "std")]
        mOS, sOS = stats.loc[r["day"], ("OS_UCell", "mean")], stats.loc[r["day"], ("OS_UCell", "std")]
        return (r["AD_UCell"] - mAD) / sAD - (r["OS_UCell"] - mOS) / sOS

    sc["L"] = sc.apply(zrow, axis=1)
    sc = sc[["cellkey", "run", "day", "cre_unit", "gene_unit", "guide", "AD_UCell", "OS_UCell", "L"]]
    sc.to_csv(sp_path, index=False)
    say("scores.csv.gz 落盘 n=%d（S2 checkpoint）" % len(sc))

    # ---- S3：单元功效盘点 -------------------------------------------
    kd = sc[sc["cre_unit"] != "NTC"].copy()
    um = (kd.groupby(["cre_unit", "gene_unit", "day"])
            .agg(n_cells=("L", "size"), n_run=("run", "nunique"),
                 L_mean=("L", "mean"), L_sd=("L", "std"))
            .reset_index())
    um["pass_floor"] = um["n_cells"] >= KD_CELL_FLOOR
    # robust 标签并回（s2h_unit_labels：cre_unit × day）
    lab = pd.read_csv(os.path.join(os.environ["GSE319266_S0V2"], "s2h_unit_labels.csv"))
    lab = lab[["cre_unit", "day", "robust_kd", "top_gene", "n_run"]]
    um = um.merge(lab, on=["cre_unit", "day"], how="left")
    um["robust_kd"] = um["robust_kd"].fillna(False)
    um.to_csv(um_path, index=False)

    n_units = len(um)
    n_rob = int((um["robust_kd"]).sum())
    n_pass = int((um["pass_floor"] & um["robust_kd"]).sum())
    genes_pass = um.loc[um["pass_floor"] & um["robust_kd"], "top_gene"].nunique()
    prior_genes = ["DLX5", "SMAD5", "RSPO3", "HOXC8", "HOXC9", "TCF7L1", "RARB", "NR3C1", "HMGA2"]
    gp = um[um["pass_floor"] & um["robust_kd"] & um["top_gene"].isin(prior_genes)]
    say("单元×day 行=%d｜robust 行=%d｜过功效线=%d（基因 %d 个）" % (n_units, n_rob, n_pass, genes_pass))
    say("先验 9 基因过功效线：%s" % sorted(gp["top_gene"].unique().tolist()))
    miss = [g for g in prior_genes if g not in set(gp["top_gene"])]
    say("先验中因功效出局：%s（G2 下限 n≥8，容忍 1）" % miss)


if __name__ == "__main__":
    main()
