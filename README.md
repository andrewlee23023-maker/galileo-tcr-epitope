# Galileo Hackathon — Track 1: TCR–epitope binding prediction

A scientific ML challenge for the [Galileo Hackathon](https://siderius.scitra.org/hackathon):
predict whether a T-cell receptor binds a peptide epitope, evaluated on
**epitopes the model has never seen** (epitope-held-out split with
epitope-similarity decontamination). Built on IEDB data, packaged for
[SIDERIUS](https://github.com/yuema137/SIDERIUS) v1.1.0.

**Status (Nov 12, 2026):** the submission is built and the deterministic
checker passes every step except 2a, which needs the Zenodo record (see
`HANDOFF.md`). Track 1 is due **Nov 18, 2026, 8:00am PST**. Decision:
**PUBLIC** — the task may be released and become a Track 2 challenge.

## Repo layout

| Path | What it is |
|---|---|
| `HANDOFF.md` | **Start here if you are continuing this work.** Current state, open items, exact commands, gotchas. |
| `submission/tcr-epitope-heldout/TASK_CARD.md` | The nine-section challenge brief (question, design, data, evaluation, baseline, novelty, difficulty, leakage risks, contact). |
| `submission/tcr-epitope-heldout/` | **The submission.** Kit-spec layout: `TASK_CARD.md`, `data_source.json`, `fetch_data.py`, `task/` (SIDERIUS package), `baseline/` (logreg + `baseline_result.json`). |
| `siderius_package/` | The build draft the submission was derived from (task.yaml, workflow.py, llm.json, declared/, plugins/, data/, evaluation/, baseline/). **Note:** `siderius_package/data/train.csv` is stored as 10 shards (`train.csv.00.part` … `train.csv.09.part`, ~1 MB each) because of a per-push size limit — reassemble with `cat train.csv.??.part > train.csv` (see "Data" below) and verify the md5. |
| `*.py`, `SCHEMA.md`, `PIPELINE.md`, `sample_data.csv`, `data/splits/` | The data pipeline: IEDB download → normalisation → splits → baselines → evaluation → leakage demo. `RUN_LOG.md` is the full run log. |
| `SUBMISSION_CHECKLIST.md` | Pre-submission checklist. |
| `check-report.json` | Latest deterministic checker report. |

## The task in one paragraph

2,510 peptidic epitopes (8–15 aa) from IEDB AIRR data; 176,612 TCR–epitope rows,
~50% positive. Real IEDB negatives are only 0.4% of rows, so negatives are
mostly **synthetic and labelled as such**: each positive gets one cross-paired
negative (a real TCR from another epitope, never reported to bind this one),
balanced 1:1 per epitope. Train 61,585 rows (~2,500 epitopes); validation
3,042 rows (2 epitopes); test 2,876 rows (2 epitopes). Primary metric:
macro-averaged AUPRC across held-out epitopes. Train epitopes within edit
distance ≤ 2 or 3-mer Jaccard ≥ 0.60 of a held-out epitope were removed (6
epitopes, 2,604 rows) — without that, k-NN scored 0.745 macro AUPRC entirely
from leakage; after decontamination it collapses to 0.434. The honest result:
additive 3-mer baselines are at chance on unseen epitopes (logreg validation
macro AUPRC 0.437, test 0.488). Beating them requires learning TCR × epitope
interaction. Full numbers in `RUN_LOG.md` and `submission/.../BASELINE_RESULTS.md`.

## Quickstart

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt        # pandas, numpy, scikit-learn, openpyxl
python download_iedb.py                # IEDB export (needs network; files NOT in this repo)
python normalise_airr.py data/airr/tcr_full_airr.tsv --out data/normalised_real.csv --cap-per-epitope 1000 --synthetic-negatives 1
python build_splits.py data/normalised_real.csv --out data/splits_real \
  --test-epitopes SPRWYFYYL TFEYVSQPFLMDLE --val-epitopes MPASWVMRI VLWAHGFEL \
  --cap-per-epitope 200 --skip-decontamination
python baseline.py --splits data/splits_real --model logreg
```

Check the submission (needs a SIDERIUS v1.1.0 checkout and `uv`; see `HANDOFF.md`
for the exact command, including the `TMPDIR` workaround):

```bash
"$SIDERIUS/.venv/bin/python" tools/check_submission.py submission/tcr-epitope-heldout \
  --siderius "$SIDERIUS" --data_dir siderius_package/data --retrain
```

## What is NOT in this repo (deliberately)

- **`test_labels.csv`** — the held-out test answer key. It is evaluator-only:
  not in the submission zip, not on Zenodo, and not here. Publishing it would
  make the task trivially gameable. It lives only in the local build folder.
- **`final_model.skops`** — the trained baseline model (binary). Regenerate it
  with `submission/.../baseline/run_baseline.py train` (see `HANDOFF.md`).
- **Raw IEDB exports** (`data/airr/`, `tcr_full_v3.xlsx`, `normalised_real.csv`,
  ~200 MB) — regenerate with `download_iedb.py` + `normalise_airr.py`, or fetch
  the three processed CSVs from Zenodo once published (see `data_source.json`).
- `data_check/`, `runs/`, `runs_real/` — regenerable intermediate outputs.

## Data

Source: Immune Epitope Database (IEDB), AIRR export, CC BY 4.0 (cite IEDB).
Processed files in `siderius_package/data/`: `train.csv` (61,585 rows),
`validation.csv` (3,042), `test_inputs.csv` (2,876, no labels). The same three
files will be published on Zenodo for the submission (see
`submission/.../data_source.json` for sizes and md5s).

`train.csv` is stored in this repo as 10 shards (`train.csv.00.part` …
`train.csv.09.part`) because a single 9.7 MB file exceeds the per-push size
limit of the GitHub upload. Reassemble it after cloning:

```bash
cd siderius_package/data
cat train.csv.00.part train.csv.01.part train.csv.02.part train.csv.03.part \
  train.csv.04.part train.csv.05.part train.csv.06.part train.csv.07.part \
  train.csv.08.part train.csv.09.part > train.csv
md5sum train.csv   # must print cd53b10080cda95dd762b0c77e0bcaea
```

## License

Task text and code: MIT (the task is intended for public release under MIT if
selected for Track 2). Data: IEDB, CC BY 4.0.
