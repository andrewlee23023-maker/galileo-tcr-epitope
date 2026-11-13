#!/usr/bin/env python3
"""Filter + deduplicate + epitope-held-out splits for the TCR-epitope task.

Input: canonical-schema CSV (see SCHEMA.md).
Output (in --out dir):
  train.csv, validation.csv — with labels
  test_inputs.csv — WITHOUT labels (solver-visible)
  test_labels.csv — evaluator-only
  heldout_epitopes.txt, DATA_CARD.md

Rules:
- Human alpha-beta TCRs, peptidic epitopes only.
- Dedup by (receptor_group_id, cdr3_alpha, cdr3_beta, epitope_id, mhc_allele).
- 10X-flagged rows excluded from the main splits (counted, not hidden).
- Cap pairs per epitope in train (--cap-per-epitope, default 200).
- Splits by epitope: train/val/test epitopes disjoint.
- CDR3 decontamination: drop any test row whose CDR3-beta is within
  --max-cdr3-distance (default 1 AA substitution, same length) of a train
  CDR3-beta, or shares >= --min-kmer-overlap (default 0.90) of its 3-mers
  with a train CDR3-beta. This is the second line of defence after the
  epitope split; the threshold is reported, not tuned to flatter a baseline.
"""
import argparse
from collections import Counter
from pathlib import Path

import pandas as pd


def kmers(seq, k=3):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def hamming_le(a, b, max_d):
    if len(a) != len(b):
        return False
    return sum(1 for x, y in zip(a, b) if x != y) <= max_d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("--out", default="data/splits")
    ap.add_argument("--test-epitopes", nargs="+", default=["E11", "E12"],
                    help="Epitopes held out for test (must exist in the data).")
    ap.add_argument("--val-epitopes", nargs="+", default=["E09", "E10"])
    ap.add_argument("--cap-per-epitope", type=int, default=200)
    ap.add_argument("--max-cdr3-distance", type=int, default=1)
    ap.add_argument("--min-kmer-overlap", type=float, default=0.90)
    ap.add_argument("--skip-decontamination", action="store_true",
                    help="Skip CDR3 decontamination. Use for the real data: "
                    "synthetic negatives are DEFINED by reusing a train TCR "
                    "with a new epitope, so Hamming<=1 decontamination drops "
                    "~56% of test (1,599/2,876) — mostly the synthetic "
                    "negatives the task is built on. The epitope-disjoint "
                    "split is the leakage defence; cross-epitope CDR3 reuse "
                    "is the signal, not a leak.")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, dtype=str, keep_default_na=False)
    # label / is_10x need real types
    df["label"] = pd.read_csv(args.input, usecols=["label"])["label"].astype(int)
    df["is_10x"] = df["is_10x"].astype(str).str.lower().isin(["true", "1", "yes"])

    counts = {"input_rows": len(df)}
    # 1. human only. Receptor type: paired alpha-beta is the minority in
    # real IEDB AIRR data (only ~33k of 302k cell IDs have both chains;
    # most are TRB-only). Restricting to alpha-beta would discard ~80%
    # of the data, so Beta-only and Alpha-only rows are kept — the
    # baseline uses whichever CDR3(s) are present. Corrected Nov 8, 2026.
    mask_human_ab = (df["species"] == "Homo sapiens")
    counts["dropped_non_human_or_non_ab"] = int((~mask_human_ab).sum())
    df = df[mask_human_ab].copy()
    # 2. need at least one CDR3 and an epitope
    mask_complete = ((df["cdr3_alpha"] != "") | (df["cdr3_beta"] != "")) & (df["epitope_sequence"] != "")
    counts["dropped_missing_cdr3_or_epitope"] = int((~mask_complete).sum())
    df = df[mask_complete].copy()
    # 3. 10X-flagged rows excluded from main splits
    counts["excluded_10x_flagged"] = int(df["is_10x"].sum())
    df_main = df[~df["is_10x"]].copy()
    # 4. dedup
    before = len(df_main)
    df_main = df_main.drop_duplicates(
        subset=["receptor_group_id", "cdr3_alpha", "cdr3_beta", "epitope_id", "mhc_allele"])
    counts["dropped_exact_duplicates"] = before - len(df_main)
    counts["after_filtering"] = len(df_main)

    epitopes = sorted(df_main["epitope_id"].unique())
    test_eps = [e for e in args.test_epitopes if e in epitopes]
    val_eps = [e for e in args.val_epitopes if e in epitopes and e not in test_eps]
    train_eps = [e for e in epitopes if e not in test_eps and e not in val_eps]
    assert test_eps, f"no test epitopes found in data (asked {args.test_epitopes}, have {epitopes})"
    assert val_eps, "no validation epitopes found in data"
    assert train_eps, "no train epitopes left"

    train = df_main[df_main["epitope_id"].isin(train_eps)].copy()
    val = df_main[df_main["epitope_id"].isin(val_eps)].copy()
    test = df_main[df_main["epitope_id"].isin(test_eps)].copy()

    # 4b. Epitope-similarity decontamination of train vs held-out epitopes.
    # Added Nov 8, 2026: exact-sequence splitting is NOT enough. Real data
    # showed test epitope TFEYVSQPFLMDLE has train neighbours CTFEYVSQPFLM,
    # CTFEYVSQPFLMDLE and NCTFEYVSQPFLMDL, and SPRWYFYYL has LSPRWYFYY /
    # LPRWYFYYL — k-NN scored 0.745 macro AUPRC largely by matching these
    # near-duplicates. A train epitope is dropped if it is a substring of,
    # or contains, a held-out epitope, or is within edit distance <= 2 of
    # one, or shares >= 0.60 of its 3-mers (Jaccard) with one. Count is
    # reported, not tuned. Held-out epitopes are never dropped from val/test.
    def lev(a, b, max_d=2):
        if abs(len(a) - len(b)) > max_d:
            return max_d + 1
        prev = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            cur = [i]
            for j, cb in enumerate(b, 1):
                cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
            prev = cur
            if min(prev) > max_d:
                return max_d + 1
        return prev[-1]

    heldout = test_eps + val_eps
    held_kmers = {h: kmers(h) for h in heldout}
    bad_train_eps = set()
    for ep in train_eps:
        ek = kmers(ep)
        for h in heldout:
            hk = held_kmers[h]
            jac = len(ek & hk) / max(len(ek | hk), 1) if ek and hk else 0.0
            if ep in h or h in ep or lev(ep, h) <= 2 or jac >= 0.60:
                bad_train_eps.add(ep)
                break
    counts["dropped_train_epitopes_near_heldout"] = len(bad_train_eps)
    counts["dropped_train_rows_epitope_decontamination"] = int(
        train["epitope_id"].isin(bad_train_eps).sum())
    train = train[~train["epitope_id"].isin(bad_train_eps)].copy()
    train_eps = [e for e in train_eps if e not in bad_train_eps]

    # 5. cap per epitope in train (keeps one dominant epitope from swamping the score)
    capped = []
    for _, g in train.groupby("epitope_id"):
        capped.append(g.sample(n=min(len(g), args.cap_per_epitope), random_state=args.seed))
    train = pd.concat(capped).sample(frac=1, random_state=args.seed)

    # 6. CDR3 decontamination of test vs train. Naive all-pairs is
    # O(test x train) and far too slow on real data (~2.9k x ~100k).
    # Use a 3-mer index: a train CDR3 can only be within Hamming <= 1
    # or Jaccard >= 0.90 of a test CDR3 if it shares at least one 3-mer
    # with it. Corrected Nov 8, 2026 — the naive loop was killed mid-run.
    # Skipped entirely with --skip-decontamination (see flag help).
    dropped_contam = 0
    if not args.skip_decontamination:
        from collections import defaultdict
        train_beta = train["cdr3_beta"].tolist()
        train_kmers = [kmers(s) for s in train_beta]
        kmer_index = defaultdict(set)
        for i, tk in enumerate(train_kmers):
            for k in tk:
                kmer_index[k].add(i)
        keep = []
        for _, row in test.iterrows():
            b = row["cdr3_beta"]
            kb = kmers(b)
            candidates = set()
            for k in kb:
                candidates.update(kmer_index.get(k, ()))
            contaminated = False
            for i in candidates:
                tb, tk = train_beta[i], train_kmers[i]
                if hamming_le(b, tb, args.max_cdr3_distance):
                    contaminated = True
                    break
                if kb and tk and len(kb & tk) / max(len(kb | tk), 1) >= args.min_kmer_overlap:
                    contaminated = True
                    break
            if contaminated:
                dropped_contam += 1
            else:
                keep.append(row)
        test = pd.DataFrame(keep) if keep else test.iloc[0:0]
    counts["dropped_test_cdr3_decontamination"] = dropped_contam

    train.to_csv(out / "train.csv", index=False)
    val.to_csv(out / "validation.csv", index=False)
    # test_inputs must not leak labels: drop label AND is_synthetic_negative
    # (a synthetic negative is label 0 by construction, so revealing the
    # flag would reveal its label).
    test.drop(columns=[c for c in ["label", "is_synthetic_negative"]
                       if c in test.columns]).to_csv(out / "test_inputs.csv", index=False)
    test[["receptor_group_id", "epitope_id", "label"]].to_csv(out / "test_labels.csv", index=False)
    (out / "heldout_epitopes.txt").write_text(
        "test: " + ",".join(test_eps) + "\nvalidation: " + ",".join(val_eps) + "\n")

    per_ep = df_main.groupby("epitope_id").size().to_dict()
    per_mhc = df_main.groupby("mhc_allele").size().to_dict()
    card = ["# Data card — TCR-epitope, epitope-held-out", "",
            f"Input rows: {counts['input_rows']}",
            f"After filtering/dedup: {counts['after_filtering']}",
            f"Dropped non-human / non-alpha-beta: {counts['dropped_non_human_or_non_ab']}",
            f"Dropped missing CDR3/epitope: {counts['dropped_missing_cdr3_or_epitope']}",
            f"Excluded 10X-flagged: {counts['excluded_10x_flagged']}",
            f"Dropped exact duplicates: {counts['dropped_exact_duplicates']}",
            f"Dropped test rows by CDR3 decontamination: {counts['dropped_test_cdr3_decontamination']}",
            f"Dropped train epitopes near held-out (epitope-similarity decontamination): {counts.get('dropped_train_epitopes_near_heldout', 0)}",
            f"Dropped train rows by epitope decontamination: {counts.get('dropped_train_rows_epitope_decontamination', 0)}",
            f"Train epitopes ({len(train_eps)}): {', '.join(train_eps)}",
            f"Validation epitopes ({len(val_eps)}): {', '.join(val_eps)}",
            f"Test epitopes ({len(test_eps)}): {', '.join(test_eps)}",
            f"Train rows: {len(train)} ({train['label'].mean():.1%} positive)",
            f"Validation rows: {len(val)} ({val['label'].mean():.1%} positive)",
            f"Test rows: {len(test)} ({test['label'].mean():.1%} positive)" if len(test) else "Test rows: 0",
            f"Cap per epitope in train: {args.cap_per_epitope}",
            f"CDR3 decontamination: Hamming <= {args.max_cdr3_distance} or 3-mer Jaccard >= {args.min_kmer_overlap}",
            "", "## Pairs per epitope (after filtering)",
            *[f"- {k}: {v}" for k, v in sorted(per_ep.items())],
            "", "## Pairs per MHC allele (after filtering)",
            *[f"- {k or '(blank)'}: {v}" for k, v in sorted(per_mhc.items())],
            "",
            "Licence: IEDB data, CC BY 4.0 — cite IEDB. Sample data (if used) is synthetic.",
            ]
    (out / "DATA_CARD.md").write_text("\n".join(card))
    print(f"train={len(train)} val={len(val)} test={len(test)} "
          f"train_eps={train_eps} val_eps={val_eps} test_eps={test_eps}")
    print(f"decontamination dropped {dropped_contam} test rows; wrote {out}")


if __name__ == "__main__":
    main()
