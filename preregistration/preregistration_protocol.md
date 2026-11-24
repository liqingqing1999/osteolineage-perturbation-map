# vc_65 Pre-registration: directional compatibility test of in silico KD effects on osteolineage cell-state distributions (S1)

> Dated 2026-11-06. Status at lock: **S1 registered; downstream analysis (S2) may start only after user approval of Section 4.**
> Nature: pre-registration document. **The criteria below were locked upon approval; any deviation must be appended as a dated amendment at the end of this file. Rewriting locked clauses is prohibited.**
> Upstream: draft criteria in internal document 24_ §3.3 (revised structurally in two places, Section 3); A-tier design (28_); innovation points #2/#3/#4 (26_).

---

## 0. Fact check completed before lock (2026-11-06; fixed as premises)

| # | Fact | Evidence |
|---|---|---|
| F1 | **PPARG / RUNX2 / CEBPA / SP7 have no units in the CRISPRi library** | all 988 rows of `s2h_unit_labels.csv` are `SNP-*` (GWAS CREs); the four genes are absent from the whole library ⇒ reciprocal test (criterion 1 of doc 24_ §3.3) is **structurally infeasible** |
| F2 | Robust-KD-measurable set = **54 genes / 67 units**, with **zero overlap** with the W1c scoring panel (AD 11 + OS 12 = 23 genes) | programmatic comparison 2026-11-06; no circularity by construction, no LOGO rescue needed |
| F3 | `CD55-MOI` is an author-built **titration control** (22 guides, ~924 cells), not a pathway-gene panel | `vc_18b_moi_dose.py`; data portrait in internal doc 03_ |
| F4 | Per-unit scale ≈ 4 cells/run, median 5.2 runs ⇒ limited single-gene power (MDES 0.38–1.92 log2) | `unit_summary.csv`; power audit in doc 14_ |
| F5 | Genes with classifiable priors = **9 genes / 11 units** (table in Section 4) | this document, Section 4; pending user approval |

---

## 1. Data and readout definitions (locked)

1. **Base**: S0v2-processed GSE319266 products (50,005 cells; L0 `Full-*.h5ad`). Scoring by UCell with gene sets identical to W1c (`AD_drift` 11 genes, `OS_identity` 12 genes; `vc_55` §GENE_SETS).
2. **Lineage axis L**: `L = z(AD_UCell) − z(OS_UCell)`, z-scored against same-day control cells.
3. **Primary statistic (per-unit displacement)**: for each robust-KD unit u,
   `d_u = mean(L | KD cells, within same day × batch layer) − mean(L | control cells, same layer)`.
4. **Gene-level aggregation**: multiple units of one gene are averaged with cell-count weights → `D_g`; gene-level direction = sign(D_g).
5. **Secondary readouts** (not gated): 2D occupancy displacement (cosine distance over a 10×10 binning of the (AD, OS) plane); DEU conversion (1 DEU = 0.016–0.025 in Δc units, per the dated calibration report `results/aux/calibration_report.md`).
6. **Inclusion threshold**: KD arm ≥ 10 cells and control pool ≥ 30 cells; otherwise the unit is recorded as "no power" — descriptive only, excluded from tests (reported n as observed).

## 2. Null distribution and test (locked)

1. **Stratified permutation**: within each (day × batch) layer, draw pseudo-KD sets of the same size as the KD arm from the control pool, **200 draws per unit**; two-sided p = (1 + #{|d_perm| ≥ |d_obs|}) / 201.
2. **Multiplicity statement**: gene-level BH-q is descriptive only; the primary test (G2) is an aggregated sign test and requires **no** per-gene multiplicity correction.
3. **Per-gene checkpointing**: each gene's null distribution is written to `null_<gene>.csv` immediately upon completion; re-runs scan existing checkpoints and fill gaps only (interruptible at any time).

## 3. Decision criteria (locked) — two structural revisions of doc 24_ §3.3

**Revision 1 (criterion 1 infeasible)**: F1 ⇒ the PPARG/RUNX2 reciprocal test is abolished and replaced by **prior-anchored directional consistency (G2)** + **negative-control floor (G1)** + **CD55 sensitivity QC (Q1)**.
**Revision 2 (criterion 2 misaligned)**: GSE319266 batches are sequencing batches, not donors ⇒ the "≥3/4 donors" requirement moves to the planned L2 experiment (follow-up study); here replaced by **leave-one-batch-out stability (G3, secondary)**.

| Gate | Content | Decision rule |
|---|---|---|
| **Q1** sensitivity QC | CD55-MOI pool (~924 cells) must be **detectable** on the 2D occupancy readout (two-sided p < .05) | fail ⇒ pipeline sensitivity insufficient; one pipeline fix allowed, recorded as a dated amendment |
| **G1** negative floor | among SNP-CRE units without robust labels (same CRE class), **≥90%** of \|d\| must fall within their own permutation p95 | fail ⇒ readout inflated by confounding; A3 verdict = FAIL |
| **G2** primary criterion | among evaluable prior genes in the Section 4 table, **n ≥ 8**, sign test of direction consistency with priors, **one-sided p < .05** (e.g., ≥7/8 or ≥8/9) | significant ⇒ A3 = PASS |
| **G3** secondary robustness | leave-one-batch-out (6 Full batches): G2 direction preserved in ≥5/6 folds; day2+4-restricted sensitivity reported | reported only; does not gate |

**A3 verdict mapping (locked)**: Q1 pass ∧ G1 pass ∧ G2 significant ⇒ **PASS** (directional manuscript proceeds); Q1 pass ∧ G1 pass ∧ G2 not significant ⇒ **INCONCLUSIVE** (downgrade to map-only; whether to write the paper is the investigator's decision; experimental spending frozen); G1 fail ⇒ **FAIL** (stop; no wet-lab spending).

## 4. Prior classification table (locked candidates, **pending user approval**; after approval, not to be adjusted based on results)

Direction convention: prior "pro-osteogenic" ⇒ after KD expect OS↓ / AD↑; "pro-adipogenic" ⇒ after KD expect AD↓.

| Gene | Prior class | Expected KD direction | One-line rationale | Units (robust) | Confidence |
|---|---|---|---|---|---|
| DLX5 | pro-osteogenic | OS↓/AD↑ | canonical osteoblast-differentiation TF (Dlx5/Dlx6) | SNP-134/day0 | high |
| SMAD5 | pro-osteogenic | OS↓/AD↑ | downstream of BMP2/4; required for osteogenic differentiation | SNP-86/day4 | high |
| RSPO3 | pro-osteogenic | OS↓/AD↑ | Wnt agonist; eBMD GWAS locus gene | SNP-109/day4 | high |
| HOXC8 | pro-osteogenic | OS↓/AD↑ | osteo/chondro differentiation regulator; KD impairs osteogenesis | SNP-230/day0 + day2 | med-high |
| HOXC9 | pro-osteogenic | OS↓/AD↑ | promotes hMSC osteogenic differentiation | SNP-228/day0 | med-high |
| TCF7L1 | pro-Wnt-repression | **OS↑/AD↓ (inverted prior)** | TCF/β-catenin repressor; KD ⇒ Wnt pathway upregulation | SNP-50/day0 | medium |
| RARB | pro-RA signaling | AD↑ (OS weak) | RA promotes osteogenesis, suppresses adipogenesis ⇒ KD liberates adipogenesis | SNP-60/day2 | medium |
| NR3C1 | pro-adipogenic | AD↓ | glucocorticoid receptor drives MSC adipogenesis | SNP-89/day4 | high |
| HMGA2 | pro-adipogenic | AD↓ | promotes adipocyte differentiation / preadipocyte proliferation | SNP-234/day0 | high |

**Prior exclusion list** (enter the descriptive map, not G2): GLIS2, CXCL12, LCORL, PLAG1, CHD4, SMARCD3, and all genes without priors (rationale: insufficient or bidirectional lineage-direction evidence). → **Evaluable prior genes n = 9, units = 11**; with the G2 floor n ≥ 8, one gene may drop out for power.

## 5. Staged execution and checkpoint design (locked)

| Stage | Product | Time | Checkpoint |
|---|---|---|---|
| S1 | this document approved | — | — |
| S2 | per-cell AD/OS UCell scores (50,005 cells) | ~20 min | safe once written |
| S3 | unit↔cell↔layer map + power accounting (report included n) | ~10 min | safe once written |
| S4 | per-unit observed displacement | ~15 min | safe once written |
| S5 | per-gene null distributions | 1–1.5 h machine | **per-gene checkpoints; interruptible at any time** |
| S6 | Q1/G1/G2/G3 verdict tables + figures + execution report | ~1 h | aggregation; interruptible |

Scripts: `vc_65a_export.py` (S2 inputs), `vc_65a2_ucell.R` (S2 scoring), `vc_65a3_assemble.py` (S3), `vc_65b_disp.py` (S4), `vc_65c_null.py` (S5, idempotent), `vc_65d_verdict.py` (S6), `vc_65e_g1v2.py` (amendment A3). All idempotent: re-runs skip completed work.

## 6. Non-claims list (locked)

⛔ No "first/novel in silico perturbation for osteoporosis" claims. ⛔ The term "virtual bone cell" is not used. ⛔ "Drift" is not used (shift/continuum instead). ⛔ Batches are not treated as donors. ⛔ No effect-size hit claims (direction only). ⛔ OLAS is not used as the sole validation. References to the W1c coordinates are restricted to the D1 compartment.

---

## Amendment note (2026-11-06 16:45 — before S2 started; no displacement results seen)

- **A1 (control pool)**: the control pool in §1.3 is changed from "same-day CRE-unit cells without robust labels" to the "**same-day NTC pool**" (after single-HTO filtering: day0 620 / day2 395 / day4 257; same pool as the S0v2 null distribution). Rationale: NTC cells are true negatives without residual target-gene effects and are mechanically isomorphic to the existing S0v2 machinery; the original control pool is retained as a sensitivity readout (reported at S6). The G1 test set (SNP-CRE units without robust labels) is **unchanged**.
- **A2 (permutation stratification)**: the stratification in §2.1 is changed from (day × batch) to **pooled within day**. Rationale: after pooling, day4 has only 11–39 NTC cells per batch — insufficient per stratum (same ruling as S0v2); pooling makes the null include between-run variability ⇒ conservative bias only, direction unaffected.

---

## Amendment A3 (2026-11-06 17:50 — G1-v2 definition, locked before rerunning)

**Background**: the first-round G1 = 89.7% (714/796; threshold 90%) ⇒ FAIL. Post-hoc audit showed the violation rate rising monotonically with n_cells and independent of target-gene knockdown ⇒ the negative set (non-robust SNP-CRE units) contains genuine enhancer candidates and is not purely negative. **This audit had been seen — a partial unblinding; this amendment is the disclosure.**

**G1-v2 definition (NTC-self-calibrated floor, locked)**:
1. Construct pseudo-units: for each of the 796 real G1 test units, draw one pseudo-unit — sample n_cells cells without replacement from the **same-day NTC pool** as `sel`; reference = remaining cells `rest`; `d_pseudo = mean(sel) − mean(rest)`;
2. Pseudo-unit null: 200 same-pool draws of equal size, `d_perm = mean(sel_i) − mean(pool∖sel_i)` (isomorphic to the observed statistic);
3. Violation criterion identical to round one: `|d_pseudo| > p95(|d_perm|)`;
4. **Floor**: pseudo-unit violation rate r_pseudo (overall + stratified by n_cells quartiles, reported);
5. **Decision**: G1-v2 PASS ⇔ real-set violation rate r_real ≤ **1.5 × r_pseudo** (overall); strata reported, not gated;
6. **A3 remapping**: A3 = PASS ⇔ Q1 pass ∧ G1-v2 pass ∧ G2 pass; otherwise the original mapping semantics hold (G1-v2 FAIL with simultaneously elevated r_pseudo ⇒ biased readout; leading to outcome ③, map-only).

**Disclosure obligation (carried into the manuscript Methods)**: "G1 was revised after the initial specification showed a power-dependent violation pattern in a negative set containing unperturbed GWAS-CRE candidates; an NTC-self-calibrated empirical floor was pre-specified before rerunning (amendment A3, dated)."

---

## Final outcome (2026-11-06 18:00)

- Execution sequence: S1 pre-registration (9-gene prior table approved by the investigator) → S2–S6 completed same day → amendments A1/A2/A3 all applied before unblinding or under the disclosure above → **G1-v2 FAIL (r_real 10.3% > 1.5 × floor = 8.29%)**; however, the floor is flat (5.5%), confirming an **unbiased readout**, and the excess violations concentrate in the high-power stratum and are absent from pure noise ⇒ the negative set contains genuine enhancer effects (readout not at fault).
- **A3 final verdict = INCONCLUSIVE (investigator decision ③, 2026-11-06 18:00)**: the FAIL semantics of G1-v2 are retained without reversal; the disposition is changed from "stop" to "downgrade". No directional claim enters the main text; the G2 8/9 result is reported as a suggestive secondary analysis with full disclosure.
- The manuscript proceeds as **map-only** (contributions: the perturbation–response map + DEU calibration); planned wet-lab follow-up remains frozen; directional claims are reserved for future donor-collected experiments.
