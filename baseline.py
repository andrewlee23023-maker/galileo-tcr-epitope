#!/usr/bin/env python3
"""Baselines for the TCR-epitope task.

Features: 3-mer counts over CDR3-alpha + CDR3-beta + epitope sequence,
plus MHC allele as a categorical feature. Two models:
  knn  — k-nearest neighbours on the k-mer vectors (TITAN-style k-NN baseline)
  logreg — logistic regression on the same vectors (the honest simple model)

Trains on train.csv, scores validation.csv, and (only via --predict-test)
writes test predictions. Test labels are never used for training or tuning.
Saves a run record JSON with config, counts and validation metrics.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import FunctionTransformer


def seq_text(df):
    return (df["cdr3_alpha"].fillna("") + " " + df["cdr3_beta"].fillna("") + " "
            + df["epitope_sequence"].fillna("") + " " + df["mhc_allele"].fillna("")).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", default="data/splits")
    ap.add_argument("--model", choices=["knn", "logreg"], default="logreg")
    ap.add_argument("--predict-test", action="store_true")
    ap.add_argument("--out", default="runs")
    args = ap.parse_args()

    splits = Path(args.splits)
    train = pd.read_csv(splits / "train.csv", dtype=str, keep_default_na=False)
    val = pd.read_csv(splits / "validation.csv", dtype=str, keep_default_na=False)
    y_train = pd.read_csv(splits / "train.csv", usecols=["label"])["label"].astype(int).values
    y_val = pd.read_csv(splits / "validation.csv", usecols=["label"])["label"].astype(int).values

    # character 3-mers over the concatenated sequences; MHC token kept as a word
    vec = CountVectorizer(analyzer="char_wb", ngram_range=(3, 3), min_df=2)
    if args.model == "knn":
        clf = KNeighborsClassifier(n_neighbors=5, weights="distance")
    else:
        clf = LogisticRegression(max_iter=1000, C=1.0)
    model = make_pipeline(vec, clf)
    model.fit(seq_text(train), y_train)

    def predict_proba(df):
        if args.model == "knn":
            return model.predict_proba(seq_text(df))[:, 1]
        return model.predict_proba(seq_text(df))[:, 1]

    p_val = predict_proba(val)
    rec = {
        "model": args.model,
        "n_train": len(train), "n_val": len(val),
        "val_auprc": float(average_precision_score(y_val, p_val)),
        "val_auroc": float(roc_auc_score(y_val, p_val)),
        "val_mean_p_pos": float(p_val[y_val == 1].mean()) if (y_val == 1).any() else None,
        "val_mean_p_neg": float(p_val[y_val == 0].mean()) if (y_val == 0).any() else None,
    }
    print(f"{args.model}: val AUPRC={rec['val_auprc']:.3f} AUROC={rec['val_auroc']:.3f} "
          f"(n_train={len(train)}, n_val={len(val)})")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.predict_test:
        test = pd.read_csv(splits / "test_inputs.csv", dtype=str, keep_default_na=False)
        p_test = predict_proba(test)
        pd.DataFrame({"receptor_group_id": test["receptor_group_id"],
                      "epitope_id": test["epitope_id"],
                      "p_binding": p_test}).to_csv(out / f"test_predictions_{args.model}.csv", index=False)
        rec["n_test"] = len(test)
        print(f"wrote {out / f'test_predictions_{args.model}.csv'} ({len(test)} rows)")
    (out / f"run_{args.model}.json").write_text(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
