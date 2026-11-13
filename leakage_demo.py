#!/usr/bin/env python3
"""Leakage demo: same baseline, random split vs epitope-held-out split.

Rebuilds splits twice from the same normalised input — once with a random
row split (the leaky way), once with the epitope-held-out split — and trains
the logistic-regression baseline on each. Prints both validation AUPRCs.

If the random-split score is far above the held-out score, the task is
measuring memorisation under a random split, and the epitope-held-out
design is doing its job. Include this output in the Track 1 package.
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.pipeline import make_pipeline

HERE = Path(__file__).parent


def seq_text(df):
    return (df["cdr3_alpha"].fillna("") + " " + df["cdr3_beta"].fillna("") + " "
            + df["epitope_sequence"].fillna("")).tolist()


def score(train, test):
    vec = CountVectorizer(analyzer="char_wb", ngram_range=(3, 3), min_df=2)
    model = make_pipeline(vec, LogisticRegression(max_iter=1000, C=1.0))
    model.fit(seq_text(train), train["label"].astype(int))
    p = model.predict_proba(seq_text(test))[:, 1]
    return average_precision_score(test["label"].astype(int), p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="normalised CSV (canonical schema)")
    ap.add_argument("--test-epitopes", nargs="+", default=["E11", "E12"])
    ap.add_argument("--val-epitopes", nargs="+", default=["E09", "E10"])
    args = ap.parse_args()

    df = pd.read_csv(args.input, dtype=str, keep_default_na=False)
    df["label"] = pd.read_csv(args.input, usecols=["label"])["label"].astype(int)
    df = df[(df["species"] == "Homo sapiens")]
    df = df[((df["cdr3_alpha"] != "") | (df["cdr3_beta"] != "")) & (df["epitope_sequence"] != "")]
    df = df.drop_duplicates(subset=["receptor_group_id", "cdr3_alpha", "cdr3_beta",
                                    "epitope_id", "mhc_allele"])

    # random split (leaky)
    rnd = df.sample(frac=1, random_state=7)
    n = len(rnd)
    r_train, r_test = rnd.iloc[:int(0.8 * n)], rnd.iloc[int(0.8 * n):]
    rnd_score = score(r_train, r_test)

    # epitope-held-out split, via build_splits.py
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([sys.executable, str(HERE / "build_splits.py"), args.input,
                        "--out", td, "--test-epitopes", *args.test_epitopes,
                        "--val-epitopes", *args.val_epitopes,
                        "--skip-decontamination"],
                       check=True, capture_output=True, text=True)
        val = pd.read_csv(Path(td) / "validation.csv")
        train = pd.read_csv(Path(td) / "train.csv")
        held_score = score(train, val)

    print(f"random split      AUPRC = {rnd_score:.3f}  (leaky: same epitopes in train and test)")
    print(f"epitope-held-out AUPRC = {held_score:.3f}  (validation epitopes unseen in train)")
    print(f"gap = {rnd_score - held_score:+.3f}")
    if rnd_score - held_score > 0.10:
        print("Interpretation: a large gap confirms the random split inflates the score — "
              "the held-out design is measuring generalisation.")
    else:
        print("Interpretation: small gap on this data — check whether the sample motifs "
              "are too weak, or the held-out epitopes too similar to train epitopes.")


if __name__ == "__main__":
    main()
