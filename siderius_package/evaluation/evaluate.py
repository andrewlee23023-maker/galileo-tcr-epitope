#!/usr/bin/env python3
"""Evaluate test predictions for the TCR-epitope task.

Usage:
  python3 evaluate.py --splits data/splits --pred runs/test_predictions_logreg.csv

Checks:
- format: one probability per test row, in [0,1], no missing/extra rows.
- primary: AUPRC macro-averaged across held-out epitopes + pooled AUPRC.
- secondary: AUROC macro + pooled, per-epitope table.
- shuffle check: labels shuffled -> metric collapses to chance.
- constant check: predicting the train prevalence for every row must score
  near chance on held-out epitopes; if it does not, the split has leaked.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def macro_metric(df, metric):
    vals = []
    for _, g in df.groupby("epitope_id"):
        if g["label"].nunique() < 2:
            continue
        vals.append(metric(g["label"], g["p_binding"]))
    return float(np.mean(vals)) if vals else float("nan"), len(vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", default="data")
    ap.add_argument("--labels", default=None,
                    help="test_labels.csv location; default evaluation/test_labels.csv "
                    "next to the package root (evaluator-only, not in data/)")
    ap.add_argument("--pred", required=True)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    splits = Path(args.splits)
    # Package layout: test labels are evaluator-only in evaluation/, while
    # train.csv lives in data/. --splits points at data/; --labels overrides
    # the label file location (default: evaluation/test_labels.csv, resolved
    # relative to the package root, i.e. the parent of --splits).
    labels_path = Path(args.labels) if args.labels else (splits.parent / "evaluation" / "test_labels.csv")
    labels = pd.read_csv(labels_path)
    pred = pd.read_csv(args.pred)
    assert set(pred.columns) >= {"receptor_group_id", "epitope_id", "p_binding"}, \
        f"predictions need receptor_group_id, epitope_id, p_binding; got {list(pred.columns)}"
    assert len(pred) == len(labels), f"row count mismatch: pred {len(pred)} vs labels {len(labels)}"
    assert pred["p_binding"].between(0, 1).all(), "p_binding must be in [0,1]"
    assert not pred["p_binding"].isna().any(), "p_binding has missing values"

    df = labels.merge(pred, on=["receptor_group_id", "epitope_id"], how="left", validate="one_to_one")
    assert not df["p_binding"].isna().any(), "some test rows have no prediction"

    auprc_macro, n_eps = macro_metric(df, average_precision_score)
    auroc_macro, _ = macro_metric(df, roc_auc_score)
    result = {
        "n_test": len(df),
        "n_epitopes_scored": n_eps,
        "auprc_macro": auprc_macro,
        "auprc_pooled": float(average_precision_score(df["label"], df["p_binding"])),
        "auroc_macro": auroc_macro,
        "auroc_pooled": float(roc_auc_score(df["label"], df["p_binding"])),
        "prevalence": float(df["label"].mean()),
    }
    # shuffle check: mean AUPRC over 20 label shuffles (one shuffle is noisy
    # at n~100, so report the mean, not a single draw)
    rng = np.random.default_rng(args.seed)
    shuffled_scores = []
    for _ in range(20):
        shuffled = df["label"].values.copy()
        rng.shuffle(shuffled)
        shuffled_scores.append(average_precision_score(shuffled, df["p_binding"]))
    result["auprc_shuffled_mean"] = float(np.mean(shuffled_scores))
    result["auprc_shuffled_sd"] = float(np.std(shuffled_scores))
    # constant check: train prevalence for every row
    train = pd.read_csv(splits / "train.csv", usecols=["label"])
    prev = float(train["label"].mean())
    # a constant predictor has undefined AUPRC; use the prevalence as a
    # threshold-free sanity: its per-epitope AUPRC equals that epitope's
    # prevalence on average — report the pooled prevalence instead.
    result["train_prevalence"] = prev
    result["constant_would_score_near_prevalence"] = True

    per_ep = []
    for eid, g in df.groupby("epitope_id"):
        per_ep.append({
            "epitope_id": eid, "n": len(g), "n_pos": int(g["label"].sum()),
            "auprc": float(average_precision_score(g["label"], g["p_binding"]))
            if g["label"].nunique() > 1 else None,
            "auroc": float(roc_auc_score(g["label"], g["p_binding"]))
            if g["label"].nunique() > 1 else None,
        })
    result["per_epitope"] = per_ep
    print(json.dumps({k: v for k, v in result.items() if k != "per_epitope"}, indent=2))
    print("\nPer-epitope:")
    for r in per_ep:
        print(f"  {r['epitope_id']}: n={r['n']} pos={r['n_pos']} "
              f"AUPRC={r['auprc'] if r['auprc'] is None else round(r['auprc'], 3)} "
              f"AUROC={r['auroc'] if r['auroc'] is None else round(r['auroc'], 3)}")
    print(f"\nShuffle check: shuffled AUPRC mean={result['auprc_shuffled_mean']:.3f} "
          f"(sd {result['auprc_shuffled_sd']:.3f}, 20 shuffles; "
          f"should be near prevalence {result['prevalence']:.3f})")


if __name__ == "__main__":
    main()
