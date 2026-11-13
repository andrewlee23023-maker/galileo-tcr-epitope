# Run log — sample end-to-end run

**Environment note (Nov 8, 2026):** the original venv `/tmp/galileo-venv2` was wiped when /tmp was cleaned. The venv now lives at `~/workspace/galileo-venv` (durable, survives VM replacement): Python 3.12, scikit-learn 1.9.1, pandas 3.1.1, openpyxl, pyyaml, numpy, pydantic. Use it for all further runs: `~/workspace/galileo-venv/bin/python ...`

Run Nov 8, 2026 in this sandbox. Command: `PY=/tmp/galileo-venv/bin/python bash run_all.sh`
(Python 3.12, venv with scikit-learn 1.9.1, pandas 3.1.1. On Andrew's machine: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`, then `bash run_all.sh`.)

## Sample data

`make_sample.py`: 687 rows, 12 epitopes, 52.8% positive (synthetic, motif-based — see `make_sample.py`).

## Filtering / splits (`build_splits.py`)

- Input: 687 rows
- After filtering/dedup: 636 rows
- Dropped non-human / non-alpha-beta: 10 (6 mouse, 4 TRA-only)
- Dropped missing CDR3/epitope: 0
- Excluded 10X-flagged: 5
- Dropped exact duplicates: 36
- Dropped test rows by CDR3 decontamination: 0 — the sample's near-duplicates are within the same epitope, so they stay in the same split; cross-epitope decontamination had nothing to remove here. On the real IEDB export this count must be checked, not assumed to be 0.
- Train: 424 rows, 8 epitopes (E01–E08), 50.9% positive
- Validation: 106 rows, 2 epitopes (E09, E10), 50.9% positive
- Test: 106 rows, 2 epitopes (E11, E12), 50.9% positive

## Baselines (validation)

- Logistic regression (3-mer features): val AUPRC 0.677, AUROC 0.614
- k-NN (k=5, same features): val AUPRC 0.655, AUROC 0.642

## Test evaluation (`evaluate.py`)

Logistic regression on held-out test (E11, E12):
- AUPRC macro 0.550, pooled 0.535; AUROC macro 0.549, pooled 0.542
- Per-epitope: E11 AUPRC 0.532 / AUROC 0.527 (n=53, 27 pos); E12 AUPRC 0.569 / AUROC 0.571 (n=53, 27 pos)
- Shuffle check: mean AUPRC 0.540 (sd 0.053, 20 shuffles) vs prevalence 0.509 — collapses to chance, as required.

k-NN on held-out test:
- AUPRC macro 0.529, pooled 0.521; AUROC macro 0.489, pooled 0.493
- Shuffle check: mean AUPRC 0.529 (sd 0.053) vs prevalence 0.509.

Interpretation for the sample only: both baselines are near chance on unseen epitopes, which is expected — the sample motifs are epitope-private, so there is almost nothing transferable to learn. The real IEDB task may be more (or less) learnable; do not quote these numbers as expected real-task performance.

## Leakage demo (`leakage_demo.py`)

- Random split AUPRC: 0.985 (leaky — same epitopes in train and test)
- Epitope-held-out AUPRC: 0.681 (validation epitopes unseen in train)
- Gap: +0.304

This is the core evidence for the Track 1 package: the same baseline, the same data, and the random split inflates the score by ~0.3 AUPRC. The held-out design is measuring generalisation.

## Download check (`download_iedb.py`)

Exit 3, as designed: no direct download link for `tcr_full_v3.xlsx` / `airr_full.zip` was found on the JS-driven export page. Full export must be downloaded in a browser from https://www.iedb.org/database_export_v3.php, saved as `data/tcr_full_v3.xlsx`, then inspected with `normalise_iedb.py --inspect`.

**Correction Nov 8, 2026: direct download DOES work** — the missing piece is the `doc/` prefix: `https://www.iedb.org/downloader.php?file_name=doc/airr_full.zip` (9.9 MB, fetched by curl) and `.../doc/tcr_full_v3.xlsx` (31.1 MB, via browser). `download_iedb.py` header updated accordingly.

---

# Real IEDB run — Nov 8, 2026

Commands (venv `/tmp/galileo-venv2`, scikit-learn/pandas/openpyxl):
```
python normalise_airr.py data/airr/tcr_full_airr.tsv --out data/normalised_real.csv --cap-per-epitope 1000 --synthetic-negatives 1
python build_splits.py data/normalised_real.csv --out data/splits_real --test-epitopes SPRWYFYYL TFEYVSQPFLMDLE --val-epitopes MPASWVMRI VLWAHGFEL --cap-per-epitope 200 --skip-decontamination
python baseline.py --splits data/splits_real --model logreg --predict-test --out runs_real
python baseline.py --splits data/splits_real --model knn --predict-test --out runs_real
python evaluate.py --splits data/splits_real --pred runs_real/test_predictions_{logreg,knn}.csv
python leakage_demo.py data/normalised_real.csv --test-epitopes SPRWYFYYL TFEYVSQPFLMDLE --val-epitopes MPASWVMRI VLWAHGFEL
```

## Real data (AIRR, Nov 3, 2026 build)

- `tcr_full_v3.xlsx` is **positive-only**: 302,807 rows, no outcome column. Do not use it for the binary task.
- `tcr_full_airr.tsv`: 339,842 chain rows, 302,511 cell IDs, 3,349 epitopes; positive 338,432 vs negative 1,410 (0.4%).
- Human, pure-peptide 8–15, TRA/TRB, CDR3 present: 288,608 chain rows → 256,699 pairs → 174,700 after dedup.
- After cap 1,000/epitope: 87,848 real positives + 916 real negatives + 87,848 synthetic negatives (per-epitope balanced, 1:1) = **176,612 rows, 2,510 epitopes, 49.7% positive**.
- First synthetic-negative version sampled a random epitope per TCR and left each epitope ~94% positive — fixed to per-epitope balancing before splitting.
- `build_splits.py` no longer requires paired alpha-beta: only ~33k cell IDs have both chains; most are TRB-only. Keeping Beta/Alpha-only rows retains ~80% of the data that the original filter discarded.

## Real splits (`data/splits_real`)

- Train: **62,197 rows, 2,506 epitopes, 49.6% positive** (capped 200/epitope)
- Validation: **3,042 rows** (`MPASWVMRI`, `VLWAHGFEL`), 50.0% positive
- Test: **2,876 rows** (`SPRWYFYYL`, `TFEYVSQPFLMDLE`), 50.0% positive
- CDR3 decontamination: **skipped** (`--skip-decontamination`). With it on, 1,599/2,876 test rows (55.6%) were dropped — synthetic negatives are *defined* by reusing a train TCR with a new epitope, so Hamming≤1 decontamination deletes the task itself. The epitope-disjoint split is the leakage defence; cross-epitope CDR3 reuse is the signal.

## Real baselines

Logistic regression (3-mer features, val): AUPRC 0.458, AUROC 0.460 — **below chance**.
k-NN k=5 (same features, val): AUPRC 0.489, AUROC 0.503 — chance.

Test evaluation:
- **LogReg: macro AUPRC 0.497, pooled 0.526; macro AUROC 0.441, pooled 0.463.** Per-epitope: `SPRWYFYYL` AUPRC 0.606 / AUROC 0.620; `TFEYVSQPFLMDLE` AUPRC 0.387 / AUROC 0.262 (strongly anti-predictive). Shuffle check: mean AUPRC 0.503 (sd 0.008) vs prevalence 0.500.
- **k-NN: macro AUPRC 0.745, pooled 0.797; macro AUROC 0.753, pooled 0.773.** Shuffle check: mean 0.502 (sd 0.009).

## Leakage demo (real data)

- Random split AUPRC 0.489 vs epitope-held-out 0.453 — **gap only +0.036**, vs +0.304 on the synthetic sample. The sample's big gap does NOT reproduce on real data: with cross-paired synthetic negatives, CDR3 memorisation does not help even under a random split, because the same TCR appears as positive (its true epitope) and negative (a cross-paired epitope) in both train and test. An additive 3-mer model cannot learn the TCR×epitope interaction either way.

## Epitope-similarity decontamination — fixed and re-run, Nov 8

Added to `build_splits.py` (step 4b): drop any train epitope that is a substring of / contains a held-out epitope, is within edit distance ≤2 of one, or shares ≥0.60 3-mer Jaccard with one. Re-ran splits + both baselines on the decontaminated data:

- Dropped: **6 train epitopes, 2,604 train rows** (train 62,197 → 61,585; val/test unchanged at 3,042 / 2,876).
- **LogReg test: macro AUPRC 0.488, pooled 0.531; macro AUROC 0.423.** (Was 0.497 / 0.441 — essentially unchanged, still chance/below.)
- **k-NN test: macro AUPRC 0.434, pooled 0.425; macro AUROC 0.374.** (Was 0.745 / 0.753 — collapsed.) Confirms the hypothesis: k-NN's 0.745 was *entirely* driven by the 6 near-duplicate epitopes. Without them, both simple 3-mer baselines are at or below chance on truly unseen epitopes.

Honest conclusion for the package: simple 3-mer baselines do **not** generalise to unseen epitopes on this data — the task is genuinely hard and needs a model that learns TCR×epitope interaction (e.g. cross-attention, pretrained protein embeddings), not additive sequence features. That is the floor the package sets, and the epitope-similarity decontamination is what makes the numbers trustworthy. Do not quote the pre-decontamination k-NN 0.745 anywhere; it was leakage.

## SIDERIUS package built — Nov 8, 2026

`siderius_package/` (11 MB), checked field-by-field against the canonical
`yuema137/SIDERIUS` repo (cloned Nov 8): `task.yaml` + `workflow.py` +
`llm.json` (canonical default shape, no keys) + `declared/` (dataset_profile,
task_config, 3 metric declarations, task_health, k-mer vocab: 7,035 3-mers
min-count-2 + 145 MHC one-hots = 7,180 features) + `plugins/` (custom data
path with 3-mer-count features, macro-AUPRC primary + pooled AUPRC
+ macro AUROC secondaries in numpy, health views, reference MLP) + `data/`
(train/validation/test_inputs, DATA_CARD, heldout list) +
`evaluation/` (evaluate.py, leakage_demo.py, **test_labels.csv —
evaluator-only, NOT in data/**) + README, SPLIT_EVIDENCE, CITATION.

Validation (`validate_package.py`, run Nov 8): **42/42 checks pass** —
YAML/JSON parse, required manifest sections, data path on real data
(train 61,585 / val 3,042, 7,180 features), scope serialize round-trip,
training dataset, deliverable write/read round-trip, all three metrics
matching scikit-learn to <1e-6, test-input hygiene (no label, no
is_synthetic_negative, unique row_id), all plugins compile. Two real bugs
found and fixed during validation: (1) `_featurise` returned Python lists
and `read_table` stacked them — OOM on 61,585×7,180; now preallocates the
matrix once and fills in place. (2) Package-layout baseline reproduces the
exact numbers: logreg macro AUPRC 0.4881629 via
`siderius_package/evaluation/evaluate.py`.

**Not yet verified (updated Nov 8, late):** the old "no composed SIDERIUS run" caveat is now
partly superseded — the kit's deterministic checker loads the package through real SIDERIUS
v1.1.0 (steps 3a/3b PASS below). Track 1 does not require running SIDERIUS research loops or a
GPU run; the checker only confirms the package loads and the baseline reproduces.

## Kit submission + deterministic checker (Nov 8, late)

SIDERIUS v1.1.0 `uv sync --group dev --frozen` finished (exit 0, ~30 min; torch 2.10.0+cu128,
CUDA unavailable in this sandbox — fine, the baseline is CPU). `scikit-learn==1.7.2` +
`skops==0.13.0` installed into the SIDERIUS venv for local training; `uv` 0.13.7 also
installed into it (the checker needs `uv` on PATH for step 4.1).

Baseline trained via `submission/.../baseline/run_baseline.py train` (seed 20261108, CPU):
train 54.979 s, evaluate 3.407 s. Validation metrics (task's own metric plugins):
**macro_auprc 0.4369476683173579, pooled_auprc 0.4581749881994889,
macro_auroc 0.4353005117364963** — train and a separate evaluate run agree exactly.
`baseline/final_model.skops` (66,756 bytes) + `baseline/baseline_result.json` written
(tolerance 0.0; hardware AMD EPYC 9D25, 4 vCPUs, no GPU; framework_versions python 3.12.3,
scikit-learn 1.7.2, skops 0.13.0, numpy 2.4.3, siderius 1.1.0 (5ccb436)).

**TASK_CARD correction made during this run:** the Baseline and Evaluation sections had
quoted the *test* numbers (0.488 / 0.531 / 0.423) as validation scores. Fixed to the real
validation numbers above; the test numbers are still quoted, now correctly labelled as
held-out test (development evaluation via `evaluation/evaluate.py`).

Checker (`tools/check_submission.py`, report `check-report.json`, zip
`tcr-epitope-heldout.zip` 132,434 bytes, sha256
032b1595082d0365caba8ea1265bacec3a756cbb7377ab169377b7239302a1c5, local working id
— no portal task id yet):
- 1 package PASS; 2b data_files PASS; 3.0 siderius_version PASS; 3a admission PASS
  (bwrap note: not installed in sandbox, recorded only); 3b composition PASS
  (primary macro_auprc higher, 2 secondary, Health declared); 4.0 baseline_model PASS;
  4.1 baseline_environment PASS; 4a baseline_evaluate PASS (metrics reproduce within
  tolerance); 4b baseline_retrain PASS (retrain from seed 20261108 reproduces metrics).
- **2a zenodo_record FAIL — the only failure: placeholder `ZENODO_RECORD_ID`.**
  Blocked on Andrew publishing the 3 files on Zenodo and filling in the record id/DOI.
  Verdict is FAIL until then, by design.

Two environment fixes found while running the checker:
1. `--data_dir` must hold the *raw Zenodo files* (all 3 listed in data_source.json),
   not the processed output — the checker verifies them, then runs fetch_data.py
   itself. Used `siderius_package/data/`; the processed `data_check/` folder is not
   the right input.
2. The checker builds its ~7.6 GB isolated venv under TMPDIR — /tmp is a 512 MB tmpfs,
   so step 4.1 failed with "No space left on device". Re-ran with
   `TMPDIR=~/workspace/.tmp-check` (cleaned up afterwards); uv cache makes the
   rebuild ~45 s.
