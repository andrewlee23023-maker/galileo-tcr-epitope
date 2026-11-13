# Split evidence — required by SIDERIUS `docs/guides/define-a-task.md` Step 2

Before using a task package for scientific evaluation, supply and retain a
reproducible, task-owned split check. Scope payloads are opaque to SIDERIUS:
different filenames or hashes, successful execution and a passing Health
check do **not** prove that training and evaluation are independent. This
file is that evidence for the shipped configuration. A proof for these
defaults does not certify another configuration or another dataset.

## 1. Independence definition and identity key

- **Claim:** no epitope in validation or test appears in training, and no
  train epitope is a near-duplicate of a held-out epitope. The identity key
  is the **epitope peptide sequence** (`epitope_sequence` / `epitope_id`),
  not the row ID: the same TCR (receptor group / CDR3) *may* appear in
  training with a different epitope — that is the task (same TCR, different
  epitope, different label), and the synthetic negatives are built on it.
- **What is forbidden:** a held-out epitope (or a near-duplicate of one)
  in the training file; a test TCR–epitope pair already present in train
  under a different assay/publication (prevented by dedup on
  receptor-group + CDR3α + CDR3β + epitope + MHC before splitting).
- **What is shared and justified:** CDR3 sequences and receptor groups
  recur across epitopes by design. A model that memorises "this CDR3 binds"
  scores at chance on held-out epitopes, because the same CDR3 appears as a
  positive for its true epitope and as a synthetic negative for other
  epitopes. Sharing TCR identity across the split is therefore not leakage
  here; sharing epitope identity would be.

## 2. Inputs checked

- Dataset: IEDB AIRR export `tcr_full_airr.tsv`, Nov 3, 2026 build,
  downloaded Nov 8, 2026 (`downloader.php?file_name=doc/airr_full.zip`).
- Split declaration: `build_splits.py` in the parent build folder, run as
  `build_splits.py data/normalised_real.csv --out data/splits_real
  --test-epitopes SPRWYFYYL TFEYVSQPFLMDLE --val-epitopes MPASWVMRI
  VLWAHGFEL --cap-per-epitope 200 --skip-decontamination`.
- Task code revision: this package, Nov 8, 2026 (see `BASELINE_RESULTS.md`).
- Populations: train 61,585 rows / 2,500 epitopes; validation 3,042 rows /
  2 epitopes; test 2,876 rows / 2 epitopes. All ~50% positive, all nonempty
  (an empty set would pass an intersection check vacuously — checked
  explicitly below).

## 3. Real membership check (run Nov 8, 2026)

From `data/` in this package (and the identical `data/splits_real/`):

```
train epitopes:      2,500
validation epitopes: 2  (MPASWVMRI, VLWAHGFEL)
test epitopes:       2  (SPRWYFYYL, TFEYVSQPFLMDLE)
train ∩ validation:  ∅  (exact sequence)
train ∩ test:        ∅  (exact sequence)
validation ∩ test:   ∅  (exact sequence)
```

Epitope-similarity decontamination (build step 4b): any train epitope that
is a substring of, within edit distance ≤ 2 of, or sharing ≥ 0.60 3-mer
Jaccard with a held-out epitope was removed — **6 epitopes, 2,604 rows**
(train 62,197 → 61,585). After that step, no train epitope satisfies any
of those three conditions against any held-out epitope (verified by
re-running the same check on the shipped train file; the check is part of
`build_splits.py`, not a one-off script).

Reproduce:
```bash
cd track1-tcr-epitope
python -c "
import pandas as pd
tr=set(pd.read_csv('siderius_package/data/train.csv', usecols=['epitope_id'])['epitope_id'])
va=set(pd.read_csv('siderius_package/data/validation.csv', usecols=['epitope_id'])['epitope_id'])
te=set(pd.read_csv('siderius_package/data/test_inputs.csv', usecols=['epitope_id'])['epitope_id'])
print('train', len(tr), 'val', va, 'test', te)
print('train&val', tr & va, 'train&test', tr & te, 'val&test', va & te)
assert not (tr & va) and not (tr & te) and not (va & te)
print('disjoint: OK')
"
```

## 4. Failing leakage counterexample

Two deliberate leaks, both demonstrated during the build (see
`BASELINE_RESULTS.md` for the numbers):

1. **Near-duplicate epitopes left in train (pre-decontamination run).**
   Train contained `CTFEYVSQPFLMDLE`, `CTFEYVSQPFLM` and
   `NCTFEYVSQPFLMDL` next to test epitope `TFEYVSQPFLMDLE`, and
   `LSPRWYFYY` / `LPRWYFYYL` next to `SPRWYFYYL`. The k-NN baseline
   scored **macro AUPRC 0.745 / AUROC 0.753** on test — a misleading high
   score produced by epitope similarity, not generalisation. After the
   decontamination step removed those 6 epitopes, the same baseline scored
   **0.434 / 0.374**. The check that catches this leak is the
   similarity-decontamination step itself; without it, the score is
   inflated by ~0.3 AUPRC.
2. **Random split (leakage demo).** `evaluation/leakage_demo.py` trains the
   same logistic-regression baseline on a random row split (same epitopes
   in train and test) and on the epitope-held-out split. On the synthetic
   sample data the gap was +0.304 (0.985 vs 0.681); on the real data it is
   only +0.036 (0.489 vs 0.453), because cross-paired synthetic negatives
   make CDR3 memorisation useless even under a random split. The demo is
   retained because it shows *why* the gap collapsed — the negative
   construction, not a failure of the split.

## 5. Reproducibility and limitations

- Exact commands are in this package's `README.md` ("Reproduce the data
  build") and in `BASELINE_RESULTS.md`.
- **CDR3 decontamination was NOT applied, deliberately:** with
  Hamming ≤ 1 / 3-mer Jaccard ≥ 0.90 it drops 1,599/2,876 test rows
  (55.6%), because synthetic negatives reuse train TCRs by definition.
  Applying it would delete the task's own negatives and bias the test set
  toward rare TCRs. This is a limitation of the negative construction, not
  a gap in the evidence: the epitope-level check above is the operative
  independence claim.
- **What this evidence does not certify:** a different held-out epitope
  set, a different cap, or a re-run with different synthetic-negative
  sampling would need its own check. The two held-out epitopes per split
  are few; per-epitope scores are noisy.
- Synthetic negatives are not confirmed nonbinders: a cross-paired TCR
  could in principle bind its assigned epitope. The split evidence covers
  epitope independence, not the ground truth of synthetic labels.
