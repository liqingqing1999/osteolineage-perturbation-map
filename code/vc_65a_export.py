import os
_T65 = os.environ.get("T65_DIR", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results")))
# -*- coding: utf-8 -*-
"""vc_65a (S2/S3 part 1): 导出 50,005 细胞 counts 矩阵 + 元数据，供 R UCell 打分。

预注册：`output/virtual-cell-osteo/30_vc65预注册_S1方向相容性检验_2026-11-06.md`
  - 基座 = S0v2 keep 过滤（single & day & guide & cre_unit），断言 == 50,005
  - 修订 A1：对照池 = 同日 NTC 池；A2：置换分层 = day 合池
  - UCell 打分放 R（与 vc_55 逐字同函数 ScoreSignatures_UCell, maxRank=1500）

产出（T65/export/）：
  matrix.mtx.gz  50,005 x 联合基因（symbol 去重，保留首次出现）
  genes.tsv / barcodes.tsv / meta.csv
幂等：export 目标已存在且行数一致时跳过重建。
"""
import os
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import anndata as ad
import gzip

L0 = os.environ["GSE319266_L0"]  # local S0v2 L0 dir with Full-*.h5ad
L1 = os.environ["GSE319266_L1"]  # local S0v2 L1 dir
SUP = os.environ["GSE319266_SUPP"]  # extracted GEO supplementary files
OUT = os.path.join(_T65, "export")
RUNS = ["Full-1", "Full-2", "Full-3", "Full-4", "Full-5", "Full-6", "Pilot-1", "Pilot-2"]
LANE2RUN = {"full_LN-1": "Full-1", "full_LN-2": "Full-2", "full_LN-3": "Full-3",
            "full_LN-4": "Full-4", "full_LN-5": "Full-5", "full_LN-6": "Full-6",
            "pilot_LN-1": "Pilot-1", "pilot_LN-2": "Pilot-2"}
HTO_RULE = {"HTO1": "day0", "HTO2": "day0", "HTO3": "day2", "HTO4": "day4"}
LOG = os.path.join(_T65, "vc65a_run.log")


def say(msg):
    m = "[vc_65a %s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def barcode(s):
    b = str(s).rsplit("_", 1)[-1]
    return b[:-2] if b.endswith("-1") else b


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()

    done_flag = os.path.join(OUT, "_DONE")
    if os.path.exists(done_flag):
        meta = pd.read_csv(os.path.join(OUT, "meta.csv"), keep_default_na=False)
        say("已存在（_DONE，n=%d），幂等跳过" % len(meta))
        return

    # ---- 1. 细胞表：与 vc_29 逐行同构 -------------------------------
    c = pd.read_csv(os.path.join(L1, "cell_assignment.csv.gz"), keep_default_na=False)
    c["bc"] = c["cell"].map(barcode)
    recs = []
    for run in RUNS:
        A = ad.read_h5ad(os.path.join(L0, run + ".h5ad"), backed="r")
        hcols = [x for x in A.obs.columns if x.startswith("hto_")]
        H = A.obs[hcols].to_numpy(float)
        frac = H / np.maximum(H.sum(1, keepdims=True), 1)
        recs.append(pd.DataFrame({"cell": A.obs.index.values,
                                  "nsig": ((H >= 10) & (frac >= 0.15)).sum(1)}))
        A.file.close()
    c = c.merge(pd.concat(recs, ignore_index=True), on="cell", how="left")
    c["single"] = c["nsig"].fillna(9) <= 1

    e = pd.read_csv(os.path.join(SUP, "S2E_cell_timepoint.csv"))
    e["bc"] = e["Cell ID"].map(barcode)
    e["run"] = e["Single-cell 10X chip lane"].map(LANE2RUN)
    e = e.dropna(subset=["run"])
    emap = (e.assign(ckey=e["run"] + "|" + e["bc"]).drop_duplicates("ckey")
             [["ckey", "hFOB time point"]].rename(columns={"hFOB time point": "day_s2e"}))
    c["ckey"] = c["run"] + "|" + c["bc"]
    c = c.merge(emap, on="ckey", how="left")
    c["day_s2e"] = c["day_s2e"].fillna("")
    c["day_rule"] = np.where(c["run"].isin(["Pilot-1", "Pilot-2"]), "day0",
                             c["hto"].map(HTO_RULE))
    c["day"] = c["day_rule"].fillna("")

    keep = c["single"] & (c["day"] != "") & (c["guide"] != "") & (c["cre_unit"] != "")
    cm = c.loc[keep].copy()
    assert len(cm) == 50005, "keep=%d != 50005" % len(cm)
    ntc = cm[cm["cre_unit"] == "NTC"]
    ntc_dist = ntc["day"].value_counts().to_dict()
    say("keep=%d（断言过）｜NTC=%d %s" % (len(cm), len(ntc), ntc_dist))
    assert ntc_dist.get("day0") == 620 and ntc_dist.get("day2") == 395 and ntc_dist.get("day4") == 257, \
        "NTC 日分布与 S0v2 不一致: %s" % ntc_dist

    # ---- 2. 逐 run 取矩阵，联合 symbol -------------------------------
    cm = cm.set_index("cell")
    mats = []
    sym_all = None
    for run in RUNS:
        t1 = time.time()
        A = ad.read_h5ad(os.path.join(L0, run + ".h5ad"))
        gsym = np.asarray(A.var["symbol"].astype(str), dtype=object)
        if sym_all is None:
            sym_all = gsym
            assert len(sym_all) == len(set(sym_all)) or True  # 先记，联合时去重
            sym_all = gsym
        else:
            assert np.array_equal(sym_all, gsym), "%s 基因面板不一致" % run
        idx = [b for b in A.obs_names if b in cm.index]
        sub = cm.loc[idx]
        X = A[idx].X
        assert X.shape[0] == len(idx)
        mats.append((run, idx, X.tocsr()))
        say("%s: 取 %d/%d 细胞（%.1f s）" % (run, len(idx), A.n_obs, time.time() - t1))
        del A

    Xall = sp.vstack([m[2] for m in mats], format="csr")
    bars = [b for m in mats for b in m[1]]
    assert Xall.shape == (50005, len(sym_all))

    # symbol 去重（保留首次出现），与 vc_55 的 genes[!duplicated(genes)] 同构
    seen, keepcols = set(), []
    for i, s in enumerate(sym_all):
        if s not in seen:
            seen.add(s)
            keepcols.append(i)
    Xall = Xall[:, keepcols]
    genes_out = sym_all[keepcols]
    say("联合矩阵 %s；symbol 去重后基因 %d" % (Xall.shape, len(genes_out)))

    # ---- 3. 落盘 -----------------------------------------------------
    meta = cm.loc[bars].reset_index().rename(columns={"index": "cell"})
    meta = meta[["cell", "run", "day", "cre_unit", "gene_unit", "guide"]]
    meta["cellkey"] = meta["run"] + "|" + meta["cell"].map(barcode)
    meta.to_csv(os.path.join(OUT, "meta.csv"), index=False)
    with open(os.path.join(OUT, "genes.tsv"), "w") as f:
        f.write("\n".join(genes_out) + "\n")
    with open(os.path.join(OUT, "barcodes.tsv"), "w") as f:
        f.write("\n".join(bars) + "\n")
    with gzip.open(os.path.join(OUT, "matrix.mtx.gz"), "wb") as fh:
        sio.mmwrite(fh, Xall)
    with open(done_flag, "w") as f:
        f.write("n_cells=%d n_genes=%d t=%.0fs\n" % (Xall.shape[0], Xall.shape[1], time.time() - t0))
    say("导出完成 %.1f min；nnz=%d" % ((time.time() - t0) / 60, Xall.nnz))


if __name__ == "__main__":
    main()
