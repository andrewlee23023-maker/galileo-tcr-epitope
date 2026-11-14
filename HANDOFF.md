# Handoff — Galileo Track 1: TCR–epitope binding prediction

**For:** Claude (or any agent) continuing this work.
**Owner:** Andrew Lee (UPenn BSE Bioengineering '29, Goodman lab).
**Deadline:** Track 1 submissions due **Nov 18, 2026, 8:00am PST (11:00am ET)**.
**Decisions already made (do not re-ask):** both tracks; Track 2 solo;
**PUBLIC for Track 1** (task may be released under MIT and become a Track 2
challenge — which means Andrew's solo Track 2 team *cannot* challenge his own
task on Nov 21, he must pick a different one). Andrew registered for both
tracks himself; confirmation email + portal link arrived.

## 1. Current state (Nov 12, 2026)

- Submission folder `submission/tcr-epitope-heldout/` is built to the official
  kit spec (yuema137/siderius-hackathon, SPEC.md v1.0.0): `TASK_CARD.md` (9
  headings), `data_source.json`, `fetch_data.py`, `task/` (SIDERIUS package),
  `baseline/` (logreg, skops, `baseline_result.json`).
- Repo: https://github.com/andrewlee23023-maker/galileo-tcr-epitope
  (public). `train.csv` is stored there as 10 shards
  (`siderius_package/data/train.csv.00.part` … `.09.part`) — reassemble after
  cloning with `cat train.csv.?? .part > train.csv` and verify md5
  `cd53b10080cda95dd762b0c77e0bcaea`. `test_labels.csv` is NOT in the repo.
- Deterministic checker: **every step PASSES except 2a** (placeholder
  `ZENODO_RECORD_ID`). Report: `check-report.json`. Zip built: 132,434 bytes,
  sha256 `032b1595082d0365caba8ea1265bacec3a756cbb7377ab169377b7239302a1c5`.
- Baseline trained (seed 20261108, CPU, train 55 s / eval 3.4 s):
  validation macro_auprc **0.4369**, pooled_auprc 0.4582, macro_auroc 0.4353.
  Retrain reproduces metrics exactly (tolerance 0.0).
- Held-out test (development evaluation only, via `evaluation/evaluate.py`):
  logreg macro AUPRC 0.488, pooled 0.531, macro AUROC 0.423; kNN 0.434 / 0.425 /
  0.374. **Never quote the pre-decontamination kNN 0.745 — it was leakage.**
- Test epitopes: `SPRWYFYYL`, `TFEYVSQPFLMDLE`. Validation: `MPASWVMRI`,
  `VLWAHGFEL`. Andrew has not explicitly approved this list — it is frozen
  unless he says otherwise. **Do not regenerate splits without re-running the
  split-evidence check** (`SPLIT_EVIDENCE.md`).

## 2. Open items — in order

### (a) Andrew (human) — portal task ID
Sign in at https://submit.siderius.scitra.org with GitHub `andrewlee23023`
(the registered account; it must have one of his registered emails added and
verified), create a task titled "TCR–epitope binding prediction,
epitope-held-out", and return the assigned ID (`t1-NNN`). The checker takes it
as `--task_id`; `make_submission.py` names the zip after it.

### (b) Andrew (human) — publish the data on Zenodo
https://zenodo.org (free account). Upload exactly three files from
`siderius_package/data/`: `train.csv`, `validation.csv`, `test_inputs.csv`.
License **CC BY 4.0**, access **public**. Return the **published version
record ID** (not the concept ID) and DOI. **Do NOT upload `test_labels.csv`.**

### (c) Agent — fill in the record, re-check, rebuild the zip
1. Put the record ID + DOI (`10.5281/zenodo.<id>`) into
   `submission/tcr-epitope-heldout/data_source.json` (`zenodo_record_id`,
   `doi`) and into the `## Data` section of `TASK_CARD.md` (the DOI must
   appear there literally). Access must stay `"public"`, license `"cc-by-4.0"`.
2. Re-run the checker **with `--download --retrain` and `--task_id t1-NNN`**:
   target verdict PASS (exit 0). See §3 for the exact command and gotchas.
3. Rebuild the zip with `make_submission.py --task_id t1-NNN` and record the
   new sha256. Verify the portal receipt sha256 matches.

### (d) Andrew (human) — upload + public release
Upload the final zip in the portal and tick **permit public release**. Until
the deadline, re-uploads replace the current submission; only the latest is
kept. A failed-validation upload leaves the current submission unchanged.

## 3. Environment + exact commands

- SIDERIUS checkout: `~/workspace/siderius-checkout` (v1.1.0, commit 5ccb436).
  Interpreter: `$SIDERIUS/.venv/bin/python`. It has torch 2.10.0+cu128 (no GPU
  in sandbox — fine, the baseline is CPU), scikit-learn 1.7.2, skops 0.13.0,
  and `uv` 0.13.7 installed in `.venv/bin/`.
- Submission kit: `~/workspace/siderius-hackathon-kit` (SPEC.md, tools/, skills/).
- **Gotcha 1 — TMPDIR:** `/tmp` here is a 512 MB tmpfs and the checker builds
  a ~7.6 GB isolated venv under TMPDIR. Always run the checker with
  `TMPDIR=<dir-on-big-disk>` (e.g. `~/workspace/.tmp-check`, delete after).
- **Gotcha 2 — `uv` on PATH:** step 4.1 SKIPs unless `uv` is findable:
  prefix PATH with `$SIDERIUS/.venv/bin`.
- **Gotcha 3 — `--data_dir` semantics:** it must hold the *raw Zenodo files*
  (all 3 listed in `data_source.json`), not the processed output. The checker
  verifies them, then runs `fetch_data.py --source_dir` itself. For local
  runs use `siderius_package/data/`; the processed `data_check/` folder is
  the wrong input.
- **Gotcha 4 — `bwrap` missing:** the sandbox has no bubblewrap. The checker
  only records this (step 3a passes anyway). Track 2 runs need it; the checker
  does not.

```bash
REPO="$PWD"   # repo root; run the commands below from the repo root
SIDERIUS=~/workspace/siderius-checkout
KIT=~/workspace/siderius-hackathon-kit
SUB=submission/tcr-epitope-heldout

# Train / evaluate the baseline (regenerates final_model.skops + metrics.json)
mkdir -p ~/workspace/.tmp-check
"$SIDERIUS/.venv/bin/python" "$REPO/$SUB/baseline/run_baseline.py" train \
  --data_dir "$REPO/siderius_package/data" --seed 20261108 --output_dir /tmp/tcr-baseline-run
# copy /tmp/tcr-baseline-run/final_model.skops into $REPO/$SUB/baseline/
# (baseline_result.json is already written with these exact numbers — do NOT
#  edit its metrics to make a check pass; retrain instead.)

# Full checker (after Zenodo + task ID exist)
cd "$KIT"
TMPDIR=~/workspace/.tmp-check PATH="$SIDERIUS/.venv/bin:$PATH" \
  "$SIDERIUS/.venv/bin/python" tools/check_submission.py "$REPO/$SUB" \
  --task_id t1-NNN --siderius "$SIDERIUS" --download --retrain \
  --report "$REPO/check-report.json"
# (before Zenodo exists, use --data_dir "$REPO/siderius_package/data" instead
#  of --download; expect 2a FAIL on the placeholder id.)

# Build the zip
mkdir -p /tmp/tcr-uploads
"$SIDERIUS/.venv/bin/python" tools/make_submission.py "$REPO/$SUB" \
  --task_id t1-NNN --output_dir /tmp/tcr-uploads
rm -rf ~/workspace/.tmp-check
```

## 4. Numbers you must not get wrong

- **Validation (the shipped baseline, `baseline_result.json`):** macro_auprc
  0.4369476683173579, pooled_auprc 0.4581749881994889, macro_auroc
  0.4353005117364963. Seed 20261108. These were once mislabelled in TASK_CARD
  as the test numbers — that was fixed Nov 8; do not reintroduce it.
- **Test (development evaluation, evaluator-only labels):** logreg macro AUPRC
  0.488, pooled 0.531, macro AUROC 0.423; kNN 0.434. Pre-decontamination kNN
  0.745 was leakage from 6 near-duplicate train epitopes — never quote it.
- **Leakage demo (real data):** random split AUPRC 0.489 vs epitope-held-out
  0.453, gap +0.036. (The synthetic sample gap +0.304 is NOT the real one.)
- **Splits:** train 61,585 (49.6% pos, ~2,506 epitopes), validation 3,042,
  test 2,876. Decontamination dropped 6 train epitopes / 2,604 rows.
  CDR3 decontamination was deliberately SKIPPED: it would delete 55.6% of
  test rows, because synthetic negatives are defined by reusing train TCRs
  with a new epitope. 483/1,024 test CDR3-alphas and 1,210/2,557 test CDR3-betas
  appear exactly in train — by design, documented in TASK_CARD.
- **Zenodo files (exact):** train.csv 9,736,175 bytes, md5
  `cd53b10080cda95dd762b0c77e0bcaea`; validation.csv 458,980,
  `952126c4dfc34f5e4eff7de898f0a834`; test_inputs.csv 460,201,
  `4a986a3998c9e42c36bd64ebec32395d`. (test_labels.csv 127,313,
  `b0eb8f8fd20d11ed075caa7fbc366c07` — **not** for Zenodo or GitHub.)
- **Features:** 7,180 (7,035 3-mers with train+validation count ≥ 2, + 145 MHC
  one-hots). Hyperparameters: `LogisticRegression(C=1.0, max_iter=1000,
  solver='lbfgs', random_state=seed)`.

## 5. Hard rules

- **Test labels stay evaluator-only.** Never put `test_labels.csv` in the zip,
  on Zenodo, on GitHub, or in the run data folder. Never put the
  `is_synthetic_negative` flag in test inputs. Never invent scores, checksums,
  or successful checks — the checker's verdict is the verdict.
- Verify internship/task facts on the official live source before stating
  them; never pitch eligibility or deadlines from inference or mirrors.
- Kit AGENTS.md: the participant (Andrew) handles registration, Zenodo
  publication, and the final portal upload (including the public-release
  choice). The agent prepares files and explains next actions.
- Track 1 does NOT require running SIDERIUS research loops — the checker only
  confirms the package loads (steps 3a/3b). No LLM key, no GPU for the baseline.

## 6. After Track 1 ships — Track 2 prep (Nov 21, solo, 4–7pm ET)

Tasks are revealed on the day; preparation is a playbook, not a pre-picked
task: ESM embedding script, leakage-safe split function, baseline-first
workflow, report template. Pick Biology/Medicine tasks first; one task only;
baseline-first; Reward Hacking Discovery Award (documenting a reproducible
evaluation loophole) is the fallback path. Track 2 compute: up to 2 h H100 per
team, first-come first-served by registration time; $30 LLM credits are
in-person only; bring-your-own compute/API keys allowed. **Cannot challenge
his own task** (public release) — pick a different one.

## 7. Update — Nov 13, 2026 (moved to Andrew's Windows machine)

- Working copy: `C:\Users\andre\Galileo Hackathon Track 1` (clone of the repo,
  `core.autocrlf=false` + `.gitattributes` `* -text`). **CRLF hazard:** with
  autocrlf on, the CSVs check out with CRLF and every md5 is wrong. Upload to
  Zenodo only the LF files in `siderius_package/data/` (md5s verified there).
- Checker runs in **WSL Ubuntu** (step 3a requires Linux): `~/workspace/`
  holds `siderius-checkout` (v1.1.0, 5ccb436; venv pinned to sklearn 1.7.2,
  numpy 2.4.3, skops 0.13.0, uv 0.13.7), `siderius-hackathon-kit`, and a
  mirror `galileo-tcr-epitope`. Sync the Windows submission folder into the
  mirror with rsync before checking.
- **Tolerance changed 0.0 -> 1e-4** in `baseline_result.json` (metrics
  untouched). Same-seed retrain on this CPU gave macro_auprc 0.4369457 (24
  threads) / 0.4369465 (1 or 4 threads) vs 0.4369477 reported: lbfgs/BLAS
  floating-point drift up to ~3e-6. Tolerance 0 would fail 4a/4b on reviewer
  hardware. SPEC section 4.2 allows tolerance for thread-count drift.
- `final_model.skops` regenerated on this machine (it is gitignored).
- Checker (local data, `--retrain`): all PASS except 2a placeholder. Zip
  113,935 bytes, sha256 5e486ee9... (pre-Zenodo, not final).
- Docs no longer claim test labels are in the submission (`evaluation/` was
  never in the zip). `test_labels.csv` is NOT on this machine or in the repo.
- **Nov 13 final:** portal task `t1-007`; Zenodo record 24212820
  (DOI 10.5281/zenodo.24212820; an earlier record 24212768 mistakenly held the
  synthetic sample and is not used). Checker with `--download --retrain
  --task_id t1-007`: **PASS**, exit 0. Final zip `dist/t1-007.zip` (gitignored),
  113,903 bytes, sha256
  `bdcf585c64786e1018a448c6692c2d65c3576fa39594e99774ecf3de2723fb07`.
  Remaining: Andrew uploads it in the portal, ticks permit public release, and
  checks the receipt sha256 matches.
