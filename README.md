# A pre-registered perturbation–response map of the human osteolineage

Code, dated pre-registration protocol, and result files for the manuscript
*"A pre-registered perturbation–response map of the human osteolineage"*
(submitted to *Physiological Genomics*).

The study asks a single question: **when a regulatory element near a bone-mineral-density
GWAS candidate gene is perturbed (CRISPRi knockdown), does the cell-state distribution of
human osteolineage cells shift along an adipogenic–osteogenic lineage axis, and in which
direction?** The manuscript reports the resulting perturbation–response map (1,245
perturbation units; 54 robust target genes) and a pre-registered directional compatibility
test whose final verdict — including a failed negative-floor gate and a dated amendment —
is reported transparently as INCONCLUSIVE.

## Data sources

| Dataset | Source | Use |
|---|---|---|
| GSE319266 | GEO (Kahlani et al. 2026, bioRxiv 2026.02.11.705362) | CRISPRi perturbation-seq: 50,005 cells, day 0/2/4, 1,245 GWAS-CRE perturbation units |
| GSE169396 | GEO (Fu et al. 2026, SSRN 10.2139/ssrn.7477279) | Donor-matched osteolineage reference: derivation of the 11-gene adipogenic / 12-gene osteogenic scoring panel (W1c coordinates, D1 compartment) |

## Repository layout

```
code/             analysis pipeline (run order: 65a → 65a2 → 65a3 → 65b → 65c → 65d → 65e → 66 → 67)
  vc_65a_export.py      export per-cell matrices from the S0v2-processed h5ad objects
  vc_65a2_ucell.R       UCell scoring (AD 11-gene / OS 12-gene panels, maxRank=1500)
  vc_65a3_assemble.py   assemble scores, unit map, power accounting
  vc_65b_disp.py        per-unit lineage-axis displacement d_u (layered by day; batch-descriptive)
  vc_65c_null.py        stratified permutation null, 200 draws/unit, per-gene checkpoint files
  vc_65d_verdict.py     gates Q1 / G1 / G2 / G3 verdict tables
  vc_65e_g1v2.py        amendment A3: NTC-self-calibrated negative floor (G1-v2)
  vc_66_figures.py      figures F2 / F4 / S2 (+ Supplemental Table S1 regeneration)
  vc_67_figures.py      figures F1 / F3
results/          deposited outputs (see below)
preregistration/  dated pre-registration protocol with amendments A1–A3 and final verdict
figures/          manuscript figures F1–F4 and S2 (byte-identical to `code/` output)
```

## Deposited results

- `results/obs_displacement.csv` — 1,245 units × observed displacement, cell counts, gate flags
- `results/gene_level.csv` — the 9 pre-specified prior genes (gate G2 analysis table)
- `results/gene_level_all.csv` — **Supplemental Table S1**: all 54 robust target genes
  (D_g, permutation 95% CI, cell counts, permutation p)
- `results/g1v2_pseudo.csv`, `results/g1v2_verdict.tsv`, `results/verdict.tsv`,
  `results/A3_verdict.txt` — amendment-A3 floor analysis and gate verdicts
- `results/null/` — 1,245 per-unit permutation checkpoints (200 draws each)
- `results/scores.csv.gz` — per-cell AD/OS UCell scores and lineage axis L (50,005 cells)
- `results/reference/w1c_verdict.csv` — 12/12 validation verdicts of the scoring panel (W1c)
- `results/reference/calibration_report.md` — dated DEU calibration report (asserted by `vc_67`)

## Reproducing the figures

All five manuscript figures regenerate **byte-identically** from the deposited results:

```bash
pip install numpy pandas scipy matplotlib          # Python 3.13; figures need these only
export T65_DIR="$(pwd)/results"
python code/vc_66_figures.py    # F2, F4, S2  -> results/figs/
python code/vc_67_figures.py    # F1, F3      -> results/figs/
```

`T65_DIR` defaults to `<repo>/results` if unset.

Re-running the full pipeline from raw data additionally requires the GEO download and the
S0v2 preprocessing products; set `GSE319266_L0`, `GSE319266_L1`, `GSE319266_S0V2`,
`GSE319266_SUPP`, and `T65_EXPORT` as documented at the top of each script.
Software: Python 3.13 (numpy, pandas, scipy, matplotlib), R 4.6.1 (UCell) — see
`environment.txt`.

## Pre-registration and amendments

`preregistration/preregistration_protocol.md` is the dated protocol locked before any
displacement result was computed, with three dated amendments (A1/A2 before any
unblinding; A3 after partial unblinding, disclosed in the document and in the manuscript
Methods). The final gate verdicts, including the failed negative-floor gate and the
resulting INCONCLUSIVE status with the map-only analysis plan, are part of the protocol
record and are reported in the manuscript.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use this code or the deposited result tables, please cite the manuscript and this
repository (DOI: see the Zenodo badge above once minted).
