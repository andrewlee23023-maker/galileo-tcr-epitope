# TCR–epitope binding prediction, epitope-held-out

**Galileo Hackathon — Track 1 (Create a scientific challenge).**
Author: Andrew Lee, UPenn Bioengineering '29 (Goodman lab).
Submission due: **Nov 18, 2026, 8:00am PST (11:00am ET)**.

## The task in one paragraph

Predict whether a given T-cell receptor (TCR) binds a given peptide epitope presented by MHC, from sequence alone. Input: a TCR (CDR3α and/or CDR3β, V/J genes where available) + an epitope peptide (8–15 amino acids) + MHC allele. Output: a binding probability. The scientific question is whether a model generalises to **unseen epitopes** — the case that matters for a new pathogen, a new tumour neoantigen, or a new CAR-T target — rather than memorising TCRs it has already seen for a familiar epitope.

## Data — IEDB

Source: Immune Epitope Database (IEDB), https://www.iedb.org — AIRR export `airr_full.zip` from the Database Export page (direct link: `https://www.iedb.org/downloader.php?file_name=doc/airr_full.zip`, Nov 3, 2026 build). Licence: IEDB data are freely available under CC BY 4.0 since May 2017; cite IEDB (see `CITATION.md`).

**Correction to the original plan:** the Excel export `tcr_full_v3.xlsx` is **positive-only** (302,807 rows, no outcome column) and cannot support a binary task. The AIRR export is the source: 339,842 chain rows, 302,511 cell IDs, 3,349 epitopes, with `iedb_assay_outcome` — but real negatives are only **1,410 of 339,842 rows (0.4%)**.

### Filtering and negatives (decisions, not defaults)

- Human only (NCBITAXON:9606), TRA/TRB only, CDR3 present, epitope pure amino acids, length 8–15.
- Paired alpha-beta is the minority (only ~33k of 302k cell IDs have both chains); Beta-only and Alpha-only rows are kept, using whichever CDR3 is present. Restricting to paired receptors would discard ~80% of the data.
- Dedup by (receptor group, CDR3α, CDR3β, epitope, MHC): 256,699 pairs → 174,700.
- Cap positives at 1,000 per epitope (the top epitope alone has 60,677 rows).
- **Synthetic negatives are required and are labelled as such.** Each positive gets one cross-paired negative: a TCR from another epitope, assigned to this epitope, never reported to bind it. They are balanced **per epitope** (1:1), not globally — an earlier version sampled a random epitope per TCR and left each epitope ~94% positive, which would make AUPRC meaningless. The `is_synthetic_negative` flag is kept in `train.csv`/`validation.csv` so solvers can stratify, but it is **removed from `test_inputs.csv`** (it would reveal the label: a synthetic negative is label 0 by construction). Real negatives (916) are kept as-is.
- Final: **176,612 rows, 2,510 epitopes, 49.7% positive** (87,848 real positives + 916 real negatives + 87,848 synthetic negatives).

## Splits — the core design

- **Train:** 61,585 rows, 2,500 epitopes, 49.6% positive (capped 200/epitope).
- **Validation:** 3,042 rows, epitopes `MPASWVMRI` and `VLWAHGFEL`, 50.0% positive. Used for model selection. Never tune on test.
- **Test:** 2,876 rows, epitopes `SPRWYFYYL` and `TFEYVSQPFLMDLE`, 50.0% positive. Inputs in `data/test_inputs.csv` (no labels); labels evaluator-only in `evaluation/test_labels.csv`.
- **Epitope-disjoint by exact sequence, plus epitope-similarity decontamination:** any train epitope that is a substring of, within edit distance ≤ 2 of, or sharing ≥ 0.60 3-mer Jaccard with a held-out epitope was removed — 6 epitopes, 2,604 rows. Exact-sequence splitting alone was **not** enough (see below).
- **CDR3 decontamination was considered and deliberately NOT applied:** synthetic negatives are defined by reusing a train TCR with a new epitope, so Hamming ≤ 1 CDR3 decontamination drops 1,599/2,876 test rows (55.6%) — mostly the synthetic negatives the task is built on. The epitope-disjoint split is the leakage defence; cross-epitope CDR3 reuse is the signal, not a leak. Full reasoning and the failing-counterexample check are in `SPLIT_EVIDENCE.md`.
- **Deviation from the original outline, stated plainly:** the outline suggested 15–25 test epitopes with 20–50 positives each, and keeping synthetic negatives out of the test set. The real long tail cannot support that at stable sizes, and real negatives (916 total) cannot support a binary test set at all. The package therefore uses 2 test + 2 validation epitopes at 300–800 positives each, with per-epitope-balanced synthetic negatives **in** the test set, labelled in the data card. Shrinking the held-out set rather than padding it with near-duplicates was the outline's own fallback rule.

## Metrics

- **Primary: macro-averaged AUPRC** across held-out epitopes (each epitope weighted equally), plus pooled AUPRC. AUPRC, not AUROC alone, because AUROC flatters a model that ranks easy negatives correctly.
- **Secondary:** pooled AUPRC, macro-averaged AUROC.
- **Validity checks:** one probability per test row in [0,1], no missing/extra rows (`evaluation/evaluate.py` asserts this); shuffle check (labels shuffled → metric collapses to chance); constant-baseline check (predicting train prevalence scores near chance on held-out epitopes).
- **Health:** generic categorical checks on thresholded predictions (p ≥ 0.5): predicting a single class for every validation row is **blocking**; dominant-class fraction > 0.8 is **recording**.

## Baselines — what the numbers actually are

Both baselines are 3-mer-count features (CDR3α + CDR3β + epitope, plus MHC one-hot) with `baseline/baseline.py`. Full run log: `BASELINE_RESULTS.md`.

| Baseline | Validation AUPRC / AUROC | Test macro AUPRC / AUROC |
|---|---|---|
| Logistic regression | 0.458 / 0.461 | **0.488 / 0.423** |
| k-NN (k=5) | 0.499 / 0.517 | **0.434 / 0.374** |

**Do not quote the pre-decontamination k-NN score of 0.745 macro AUPRC.** It was entirely leakage from 6 near-duplicate epitopes in train (e.g. `CTFEYVSQPFLMDLE` next to test epitope `TFEYVSQPFLMDLE`); after removing them it collapsed to 0.434. That collapse is the package's strongest evidence that the evaluation design works: simple additive features are at or below chance on truly unseen epitopes, and beating them requires learning TCR × epitope interaction. The leakage demo (`evaluation/leakage_demo.py`) shows the same story from the other side: on real data the random-split vs held-out gap is only +0.036, because cross-paired synthetic negatives make CDR3 memorisation useless even under a random split.

## Package layout

```
siderius_package/
  README.md            # this file
  SPLIT_EVIDENCE.md    # required split evidence (independence, membership, counterexample)
  CITATION.md          # IEDB citation + licence
  BASELINE_RESULTS.md  # full run log (sample + real runs)
  task.yaml            # SIDERIUS manifest (checked against canonical yuema137/SIDERIUS)
  workflow.py          # experiment declaration
  llm.json             # LLM profile (no API keys; supply your own)
  declared/            # dataset_profile, task_config, metric declarations, health, k-mer vocab
  plugins/             # data path, metrics, health views, reference MLP
  data/                # train.csv, validation.csv, test_inputs.csv, DATA_CARD.md, heldout_epitopes.txt
  evaluation/          # evaluate.py, leakage_demo.py, test_labels.csv (EVALUATOR-ONLY)
  baseline/            # baseline.py (logreg + k-NN on the same features)
```

## How to run

Prerequisites: Python 3.12, `pip install -r requirements.txt` (from the parent `track1-tcr-epitope/` folder), and a SIDERIUS checkout for a full run — a real SIDERIUS run needs an NVIDIA GPU with CUDA (no CPU mode) and bubblewrap.

Baselines only (CPU, no SIDERIUS needed):
```bash
cd track1-tcr-epitope
python baseline.py --splits siderius_package/data --model logreg --predict-test --out runs_pkg
# note: evaluate.py in this package reads labels from evaluation/test_labels.csv
python siderius_package/evaluation/evaluate.py --splits siderius_package/data --pred runs_pkg/test_predictions_logreg.csv
python siderius_package/evaluation/leakage_demo.py data/normalised_real.csv --test-epitopes SPRWYFYYL TFEYVSQPFLMDLE --val-epitopes MPASWVMRI VLWAHGFEL
```

Full SIDERIUS run (from a SIDERIUS checkout):
```bash
bash scripts/launch/run_chain.sh \
  --mode lilab \
  --workspace /path/to/workspace \
  --run_name tcr_epitope_v1 \
  --task_composition /path/to/siderius_package/task.yaml \
  --data_dir /path/to/siderius_package/data \
  --llm_config /path/to/siderius_package/llm.json \
  --healthgate_mode blocking \
  --result_authority scientific \
  --num_iterations 2 --max_rounds 3
# Start with --dry-run to see the exact child commands before spending anything.
```

Reproduce the data build from the raw export (parent folder):
```bash
python normalise_airr.py data/airr/tcr_full_airr.tsv --out data/normalised_real.csv --cap-per-epitope 1000 --synthetic-negatives 1
python build_splits.py data/normalised_real.csv --out data/splits_real --test-epitopes SPRWYFYYL TFEYVSQPFLMDLE --val-epitopes MPASWVMRI VLWAHGFEL --cap-per-epitope 200 --skip-decontamination
```

## Known limitations

- Only 2 test epitopes and 2 validation epitopes: the long tail cannot support more at stable sizes. Per-epitope scores are therefore noisy; the macro average over 2 epitopes is the primary score by design, not by choice.
- Synthetic negatives are not experimental nonbinders: a cross-paired TCR *could* bind its assigned epitope in reality (it was simply never reported to). The task measures discrimination against plausible nonbinders, not against confirmed ones.
- MHC is often generic (`HLA class I` / `HLA class II`) rather than a specific allele; unseen MHC alleles map to the all-zero one-hot vector.
- The reference model (3-mer MLP) is a floor, not a target: it is expected to score near chance.
- A real SIDERIUS run of this package has **not** been executed end-to-end here (no GPU/CUDA in this sandbox). The manifest was checked field-by-field against the canonical repo and template, the plugins compile, and the data path + metrics were exercised on the real data with stubbed SIDERIUS imports — but a full composed run on GPU is still unverified. See `SUBMISSION_CHECKLIST.md` in the parent folder.
