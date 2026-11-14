# Task card: TCR–epitope binding prediction, epitope-held-out

## Scientific question

Can a model predict whether a T-cell receptor (TCR) binds a given peptide epitope presented by MHC, for an epitope it has never seen? Each sample is one TCR–epitope pair: CDR3-alpha and/or CDR3-beta sequences, an epitope peptide (8–15 amino acids), and an MHC allele. The answer matters whenever the epitope is new — a novel pathogen, a tumour neoantigen, or a new CAR-T target — because a model that only recognises familiar epitopes is useless exactly when it is needed most.

## Task design

One sample is one TCR–epitope pair, labelled binding (1) or non-binding (0). The model input is 7,180 numeric features per row: 3-mer counts over CDR3-alpha, CDR3-beta and the epitope peptide (7,035 features, vocabulary in `task/declared/kmer_vocab.json`), plus a one-hot MHC allele (145 alleles; an unseen allele maps to the all-zero vector). The output is two logits, class 1 = binding; the deliverable is one binding probability per row in [0,1].

Splits are by epitope, not by row. Train: 61,585 rows, 2,500 epitopes. Validation: 3,042 rows, epitopes `MPASWVMRI` and `VLWAHGFEL` — used for model selection, never for tuning on test. Test: 2,876 rows, epitopes `SPRWYFYYL` and `TFEYVSQPFLMDLE`. Train, validation and test epitopes are disjoint by exact sequence, and train epitopes within edit distance ≤ 2 / 3-mer Jaccard ≥ 0.60 of any held-out epitope were removed (6 epitopes, 2,604 rows). This split matches the scientific use: the test epitopes are genuinely new cases, the way a new pathogen's epitopes would be. `task.yaml` implements it through `task_data_path` (train/validation files), `dataset_profile` (held-out epitope lists in `topology`), and the decontamination step in the build scripts (documented in `SPLIT_EVIDENCE.md`).

## Data

- DOI: 10.5281/zenodo.24212820
- Size: 10,655,356 bytes total on Zenodo, 3 files — train.csv 9,736,175 bytes (61,585 rows), validation.csv 458,980 bytes (3,042 rows), test_inputs.csv 460,201 bytes (2,876 rows, no labels)
- License: cc-by-4.0
- Access: public
- Processing: `fetch_data.py` downloads the three files from Zenodo, verifies each file's size and md5 against `data_source.json`, and writes `train.csv` and `validation.csv` into the run data folder. `test_inputs.csv` is verified but not written there — it is the held-out test input. Test labels are not on Zenodo; they are evaluator-only, held by the task author, and are not included in this submission, on Zenodo, or anywhere public.

Source: Immune Epitope Database (IEDB), AIRR export (`tcr_full_airr.tsv`, Nov 3, 2026 build). Human TRA/TRB only, CDR3 present, peptidic epitopes 8–15 amino acids, deduplicated by receptor group + CDR3 + epitope + MHC, positives capped at 1,000 per epitope. Known quality issues, stated plainly: real IEDB non-binding outcomes are only 916 of 339,842 chain rows (0.4%), so negatives are mostly synthetic — each positive gets one cross-paired negative (a TCR from another epitope, never reported to bind this one), balanced 1:1 per epitope and flagged `is_synthetic_negative` in train/validation (the flag is removed from test inputs, where it would reveal the label). A cross-paired TCR could in principle bind its assigned epitope; the task measures discrimination against plausible non-binders, not confirmed ones. MHC is often generic (`HLA class I` / `HLA class II`) rather than a specific allele. See `README.md` and `data/DATA_CARD.md` in the build folder.

## Evaluation

Primary metric: `macro_auprc` (higher is better) — average precision, macro-averaged across the held-out epitopes in the evaluation split, each epitope weighted equally. AUPRC, not AUROC alone, because positives are rare per epitope and AUROC flatters a model that ranks easy negatives correctly. Secondary metrics: `pooled_auprc` (average precision over all validation rows pooled) and `macro_auroc` (mean per-epitope AUROC). All three are computed by the task's own metric plugins (`task/plugins/_tcr_metrics.py`) from binding probabilities.

Health checks (generic categorical checks on thresholded predictions, p ≥ 0.5): predicting a single class for every validation row is **blocking** — a collapsed classifier's score says nothing. Dominant-class fraction above 0.8 is **recording**. A trivial predictor scores at chance by design: a constant predictor gets macro AUPRC ≈ 0.50 (the validation prevalence), and the shipped logistic-regression baseline gets macro AUPRC 0.437 on validation (0.488 on the held-out test) — the task is hard, and the baseline is the floor, not the target.

## Baseline

Method: scikit-learn logistic regression (`LogisticRegression`, C=1.0, max_iter=1000, lbfgs) on the task's 7,180 3-mer + MHC features; no initial model (seed + hyperparameters reproduce it). Model file: `baseline/final_model.skops` (skops format, loads without unpickling arbitrary objects). Extra packages: `baseline/requirements.txt` pins `scikit-learn==1.7.2` and `skops==0.13.0`. The bonus SIDERIUS model plugin (`task.yaml` `model_plugins`, PyTorch reference MLP `tcr_reference_mlp`) is provided separately and is not the shipped baseline.

Scores on the validation split (from `baseline/baseline_result.json`): macro_auprc 0.437, pooled_auprc 0.458, macro_auroc 0.435. Seed 20261108. Hardware: CPU (exact CPU model recorded in `baseline_result.json`). Reproduction tolerance is 1e-4 (absolute) per metric: retraining with the same seed on other CPUs or BLAS thread counts moved the scores by up to 3e-6 (lbfgs floating-point differences), so 0 would fail on reviewer hardware. Wall time recorded in `baseline_result.json`. No GPU needed. On the held-out test split the same model scored macro AUPRC 0.488, pooled AUPRC 0.531, macro AUROC 0.423 (development evaluation with the build folder's `evaluation/evaluate.py`, not part of this submission; test labels are evaluator-only). A second baseline, k-NN (k=5, distance-weighted, same features), scored macro AUPRC 0.434 on test after decontamination — and 0.745 before it, a leakage artefact from near-duplicate epitopes in train that the decontamination step removes (see `SPLIT_EVIDENCE.md`). Neither baseline beats chance on truly unseen epitopes; beating them requires learning TCR × epitope interaction.

## Novelty

Existing TCR–epitope benchmarks (VDJdb, McPAS-TCR, TITAN) ship random or epitope-agnostic splits and report scores inflated by epitope memorisation; TITAN ships no epitope-held-out split at all. This task's contribution is the evaluation design, not the data: an epitope-held-out split with epitope-similarity decontamination, per-epitope-balanced synthetic negatives disclosed as synthetic, and a demonstrated leakage counterexample (k-NN 0.745 → 0.434 when the near-duplicates are removed). Closest prior work: TITAN (Weber et al., https://github.com/PaccMann/TITAN) and ERGO-II (Springer et al., https://github.com/IdoSpringer/ERGO-II).

## Difficulty

The baseline is at chance on unseen epitopes, and that is the honest result: with cross-paired synthetic negatives, the same TCR appears as a positive for its true epitope and as a negative for other epitopes, so additive sequence features cannot separate the classes — a model must learn the interaction between the TCR and the epitope. Headroom is large but requires interaction models (cross-attention, bilinear pooling, pretrained protein embeddings), not bigger additive models. Additional difficulty: only 2 validation and 2 test epitopes (the long tail cannot support more at stable sizes), so per-epitope scores are noisy; synthetic negatives are plausible rather than confirmed non-binders; MHC is often generic.

## Known shortcuts and leakage risks

1. **Near-duplicate epitopes in train.** Exact-sequence splitting is not enough: `CTFEYVSQPFLMDLE` next to test epitope `TFEYVSQPFLMDLE` inflated k-NN to 0.745 macro AUPRC. Mitigation: epitope-similarity decontamination (substring / edit distance ≤ 2 / 3-mer Jaccard ≥ 0.60), which removed 6 epitopes and collapsed the score to 0.434. The check is part of the build, and the counterexample is documented in `SPLIT_EVIDENCE.md`.
2. **CDR3 reuse across epitopes.** 483 of 1,024 unique test CDR3-alphas and 1,210 of 2,557 unique test CDR3-betas appear exactly in train too. This is by design (synthetic negatives reuse train TCRs), not a leak: the same CDR3 with a different epitope has a different label, so memorising CDR3s scores at chance. CDR3 decontamination was tested and rejected — it drops 55.6% of test rows, mostly the synthetic negatives the task is built on.
3. **`is_synthetic_negative` flag.** A synthetic negative is label 0 by construction, so the flag would reveal test labels. It is kept in train/validation for stratification and removed from test inputs.
4. **Test labels.** Kept evaluator-only, not on Zenodo, not in the run data folder. If they were published, the task would be trivially gameable.
5. **Epitope-prior shortcut.** Predicting an epitope's train prevalence scores at chance on held-out epitopes by construction (each held-out epitope is 50% positive with no train prevalence to copy).

## Contact

Andrew Lee, University of Pennsylvania (Bioengineering '29), lee59@engineering.upenn.edu. GitHub: andrewlee23023. Built for the Galileo Hackathon, Track 1, November 2026. If this task is made public for Track 2, it is published under the MIT License with the author's agreement.
