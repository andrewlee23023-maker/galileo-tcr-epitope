# Track 1 submission checklist — due Nov 18, 2026, 8:00am PST (11:00am ET)

Submission folder: `submission/tcr-epitope-heldout/` — built to the official kit spec
(https://github.com/yuema137/siderius-hackathon, SPEC.md v1.0.0), NOT a zip of
`siderius_package/` (that draft is 11 MB with data inside and violates the kit's
limits: one zip ≤ 20 MiB, any file ≤ 5 MiB, **no data in the zip** — data lives
on Zenodo). Portal: https://submit.siderius.scitra.org — sign in with the
registered GitHub account (andrewlee23023), create a task to get its ID
(`t1-NNN`), upload the zip. Re-uploads replace the current submission.

## Decisions made

- [x] **PUBLIC for Track 1** (Andrew, Nov 8, 9:56pm ET). The task may be publicly released and become a Track 2 challenge, published under MIT. Consequence: his solo Track 2 team **cannot challenge his own task** on Nov 21 — he must pick a different task. The actual permission is granted in the portal at upload time ("Choose whether to permit public release").
- [x] Confirmation email / portal link arrived (Nov 8): https://siderius.scitra.org/hackathon/track-1 — submission is via the portal, not the registration Google Form.

## Done

- [x] IEDB data downloaded, filtered, split (train 61,585 / val 3,042 / test 2,876)
- [x] Epitope-similarity decontamination added and re-run (6 epitopes, 2,604 rows)
- [x] Baselines run on real data (logreg 0.488 / kNN 0.434 macro AUPRC after decontamination)
- [x] `task.yaml` checked field-by-field against canonical `yuema137/SIDERIUS` repo + CSV template + composition reference (Nov 8, 2026)
- [x] Package validated locally: 42/42 checks pass (`validate_package.py`)
- [x] Split evidence written (`SPLIT_EVIDENCE.md`)
- [x] Test labels kept evaluator-only (`evaluation/test_labels.csv`, not in `data/`, not on Zenodo)
- [x] Submission folder restructured to kit layout: `TASK_CARD.md` (9 headings), `data_source.json` (real sizes/md5s, record id placeholder), `fetch_data.py` (verified: downloads/verifies 3 Zenodo files, writes train+validation only), `task/` (task.yaml + declared/ + plugins/), `baseline/run_baseline.py` (sklearn logreg, skops) + `requirements.txt` (scikit-learn==1.7.2, skops==0.13.0)
- [x] `fetch_data.py --source_dir` verified locally: 3 files verified against `data_source.json`, train+validation written, test_inputs held out
- [x] SIDERIUS v1.1.0 checkout at `~/workspace/siderius-checkout` (`uv sync --group dev --frozen` running; needs torch+CUDA, several GB)

## Still to do before submitting

- [ ] **Andrew: create the task in the portal** (https://submit.siderius.scitra.org, sign in with GitHub andrewlee23023, title "TCR–epitope binding prediction, epitope-held-out") and give me the assigned task ID (`t1-NNN`). The checker needs it for step 3a; the zip builder names the zip after it.
- [ ] **Andrew: publish the data on Zenodo** (https://zenodo.org, free account). Upload exactly 3 files — `train.csv`, `validation.csv`, `test_inputs.csv` (from `siderius_package/data/`) — license CC BY 4.0, access public. Then give me the **published version record ID** (not the concept ID) and DOI; I fill them into `data_source.json` (currently placeholder `ZENODO_RECORD_ID`, which fails checker step 2a by design). Do NOT upload `test_labels.csv` to Zenodo.
- [x] **Baseline training + `baseline_result.json` + `final_model.skops`** — done Nov 8 (late): logreg, seed 20261108, train 55.0 s / evaluate 3.4 s CPU; validation macro_auprc 0.4369, pooled_auprc 0.4582, macro_auroc 0.4353 (train and evaluate agree exactly). TASK_CARD fixed: it had mislabelled the test numbers (0.488/0.531/0.423) as validation scores.
- [x] **Deterministic checker run** — done Nov 8 (late): every step PASSES except 2a (placeholder Zenodo id — needs Andrew's publication). 4b retrain reproduces metrics within tolerance 0.0. Zip built: 132,434 bytes, sha256 032b1595082d0365caba8ea1265bacec3a756cbb7377ab169377b7239302a1c5. Full report in `check-report.json`. Re-run after Zenodo with `--download` + `--retrain` and the real task id; note the checker needs `TMPDIR` on the big disk (/tmp is a 512 MB tmpfs) and `uv` on PATH.
- [ ] **Andrew: upload the final zip to the portal** and tick **permit public release**. Until the deadline you can re-upload; the new upload replaces the current one, only the latest is kept.
- [ ] **Freeze test labels.** `evaluation/test_labels.csv` is final as shipped; do not regenerate splits after this point without re-running the split evidence check.

## Notes

- The old "full SIDERIUS GPU run" item is superseded: Track 1 does NOT require running SIDERIUS research loops. The checker uses a SIDERIUS v1.1.0 checkout only to confirm the package loads (steps 3a/3b); no LLM key and no GPU unless the baseline needs one. Our baseline is CPU.
- Registration closes Nov 15, 11:59pm PST — already registered for both tracks (per Andrew, Nov 8).
