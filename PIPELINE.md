# Track 1 build — TCR–epitope binding, epitope-held-out

Galileo Hackathon, Track 1 (individual). Due Nov 18, 2026, 8:00am PST.
Task: predict whether a TCR binds a peptide epitope, from sequence alone, generalising to **unseen epitopes**.

## What is in this folder

| File | What it does |
|---|---|
| `download_iedb.py` | Tries the documented IEDB export routes; validates what it gets. See "Getting the real data" below — the full export is browser-gated. |
| `normalise_iedb.py` | Maps an IEDB export (`tcr_full_v3.xlsx` or `airr_full.tsv`) to the canonical schema in `SCHEMA.md`. `--inspect` lists the actual columns first — it never silently guesses a mapping. |
| `make_sample.py` | Generates `sample_data.csv`: a small synthetic dataset with the same schema, designed so the leakage demo shows a real gap. |
| `build_splits.py` | Filtering + dedup + epitope-held-out train/val/test splits + CDR3 decontamination. Writes `train.csv`, `validation.csv`, `test_inputs.csv` (no labels), `test_labels.csv` (evaluator-only), and `DATA_CARD.md`. |
| `baseline.py` | Two baselines: k-NN on CDR3+epitope 3-mers (TITAN-style) and logistic regression on the same features. Saves run records. |
| `evaluate.py` | Macro-averaged AUPRC/AUROC across held-out epitopes, pooled metrics, per-epitope table, shuffle check, constant-baseline check. |
| `leakage_demo.py` | Same baseline, random split vs epitope-held-out split. The gap is the evidence the task measures generalisation, not memorisation. |
| `run_all.sh` | Sample end-to-end run: make sample → splits → baselines → evaluate → leakage demo. |
| `SCHEMA.md` | Canonical column schema + IEDB mapping notes. |
| `RUN_LOG.md` | Actual counts from the sample run in this sandbox (Nov 8, 2026). |

## Getting the real data — read this first

Verified source: IEDB Database Export, https://www.iedb.org/database_export_v3.php
Weekly builds named there: `tcr_full_v3.xlsx`, `bcr_full_v3.xlsx`, `airr_full.zip`.
Also: homepage query → Receptors tab → Export Results (XLSX/TSV/CSV/JSON, or Export Type = AIRR).

**What I found in this sandbox (Nov 8, 2026): the full export could NOT be fetched by direct `curl`.** The export page is JS-driven; `curl` returns the page HTML, and guessed direct-download URLs (`downloader.php?file_name=tcr_full_v3.xlsx`, etc.) also return HTML, not a file. So do not assume a scripted download works — download `tcr_full_v3.xlsx` (or `airr_full.zip`) **in a browser** from the Database Export page, save it as `data/tcr_full_v3.xlsx`, then run:

```
python3 normalise_iedb.py --inspect data/tcr_full_v3.xlsx
python3 normalise_iedb.py data/tcr_full_v3.xlsx --out data/normalised.csv
python3 build_splits.py data/normalised.csv --out data/splits
```

`normalise_iedb.py --inspect` prints the real column names. If they differ from `SCHEMA.md`, update the mapping table in that file — do not edit the split logic to fit an uninspected file.

Verified IEDB scope (Scheffer et al., bioRxiv 2026-06-06; data as of 2026-05-28): 195,808 unique αβ TCR–epitope pairs, 181,168 TCRs, 2,849 peptidic epitopes. Licence: CC BY 4.0 since May 2017, cite IEDB. Do **not** substitute VDJdb data here — its redistribution licence was not verified.

**Correction (Nov 8, after inspecting the real files): `tcr_full_v3.xlsx` is positive-only** — 302,807 rows, 74 columns in a 2-row header, no assay-outcome column. The binding/nonbinding data live in `airr_full.zip` (`tcr_full_airr.tsv`: 339,842 chain rows, 302,511 cell IDs, 3,349 epitopes; positive 338,432 vs negative only 1,410). Use `normalise_airr.py` on that file, not `normalise_iedb.py` on the xlsx, for the binary task. Synthetic negatives (1:1, cross-paired, labelled synthetic) are required because real negatives are 0.4% of the data. Direct download that works: `https://www.iedb.org/downloader.php?file_name=doc/airr_full.zip` — the `doc/` prefix is essential.

## Quick start (sample data, no download needed)

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash run_all.sh
```

That generates the sample, builds the splits, runs both baselines, evaluates them, and runs the leakage demo. See `RUN_LOG.md` for the numbers this produced here.

## Filtering decisions (encoded in `build_splits.py`)

1. Human αβ TCRs with peptidic epitopes only. Mouse, TRA-only / TRB-only-only rows are dropped and counted.
2. Deduplicate by receptor group + CDR3α + CDR3β + epitope + MHC. One row per unique triple after dedup.
3. Labels come from IEDB assay outcomes only (Positive / Negative). No invented negatives in the test set. Synthetic random-pair negatives, if ever added, must be labelled synthetic, kept out of test, and reported separately.
4. 10X Genomics single-cell subsets flagged in the source are excluded from the main splits (they carry a known false-positive warning in the IEDB preprint) — counted, not silently dropped.
5. Dominant epitopes / HLA alleles are capped (default cap 200 pairs per epitope in train; see `--cap-per-epitope`) so one epitope does not dominate the score. The cap and the per-epitope counts go in `DATA_CARD.md`.
6. Splits are by epitope: train / validation / test epitopes are disjoint. Then CDR3 decontamination removes any test TCR whose CDR3β is within 1 amino-acid edit (or ≥90% 3-mer overlap) of a train CDR3β for a *different* epitope. The threshold is a parameter, reported in the data card — it is not tuned to make the baseline look good.

## Metrics

Primary: **AUPRC, macro-averaged across held-out epitopes** (each epitope weighted equally), plus pooled AUPRC. Secondary: AUROC (macro + pooled), balanced accuracy at the validation-chosen threshold, per-epitope table, calibration (mean predicted probability for positives vs negatives).

Validity checks in `evaluate.py`: output format (one probability per test row, in [0,1], no missing/extra rows), shuffle check (labels shuffled → metric collapses to chance), constant / epitope-prior check (predicting only train prevalence must score near chance on held-out epitopes — if it does not, the split has leaked).

## Known limits (do not hide these in the submission)

- The sample data is synthetic and motif-based; its numbers prove the pipeline runs, not that the real task is solvable at that level.
- The real held-out epitope set will be small/skewed (only 2,849 peptidic epitopes in IEDB, long tail). Fix the held-out list after the real export and shrink it rather than padding with near-duplicates.
- `task.yaml` / SIDERIUS packaging is a separate step (see the outline in the parent folder). No canonical Track 1 template was verified; re-check manifest fields against the canonical `yuema137/SIDERIUS` repo before submitting.
- A real SIDERIUS run needs an NVIDIA GPU with CUDA (no CPU mode) and bubblewrap.
