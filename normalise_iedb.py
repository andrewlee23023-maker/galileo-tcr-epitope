#!/usr/bin/env python3
"""Map an IEDB export to the canonical schema (see SCHEMA.md).

WARNING — verified Nov 8, 2026: tcr_full_v3.xlsx is POSITIVE-ONLY
(302,807 rows, no assay-outcome column). Do NOT use this script's
xlsx path for the binary task. Use normalise_airr.py on
tcr_full_airr.tsv (from airr_full.zip) instead — that file has
iedb_assay_outcome, though real negatives are only 1,410 rows.
This script is kept for inspecting the xlsx columns and for
positive-only analyses; its label mapping will fail loudly if no
outcome column is found, rather than inventing labels.

Two inputs were originally supported:
  1. tcr_full_v3.xlsx — the full TCR export from the Database Export page.
  2. airr_full.tsv — the full AIRR export (one chain per row; cell_id
     reconnects alpha and beta rows for the same receptor-epitope pair).

Run --inspect first. It prints the actual columns in the file. If a required
mapping is missing, the script exits with an error listing the available
columns instead of silently guessing. Update COLUMN_MAP below to match the
real export, then re-run.

IEDB assay outcomes: Positive -> label 1, Negative -> label 0. Rows with any
other outcome are dropped and counted.
"""
import argparse
import sys

import pandas as pd

# Column names as documented in the IEDB receptor export. VERIFY with
# --inspect on the real file before trusting this table.
COLUMN_MAP_XLSX = {
    "receptor_group_id": ["Receptor Group ID", "receptor_group_id", "IEDB Receptor ID"],
    "cdr3_alpha": ["CDR3 Alpha", "cdr3_alpha"],
    "cdr3_beta": ["CDR3 Beta", "cdr3_beta"],
    "epitope_sequence": ["Epitope", "epitope_sequence", "Epitope Sequence"],
    "epitope_id": ["Epitope ID", "epitope_id"],
    "mhc_allele": ["MHC Restriction", "mhc_allele", "MHC Allele"],
    "species": ["Species", "species"],
    "receptor_type": ["Receptor Type", "receptor_type"],
    "assay_outcome": ["Assay Outcome", "assay_outcome", "Qualitative Measure"],
    "assay_id": ["Assay ID", "assay_id"],
    "reference_id": ["Reference ID", "reference_id"],
    "assay_type": ["Assay Type", "assay_type"],
}


def pick(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("--out", default="data/normalised.csv")
    ap.add_argument("--inspect", action="store_true",
                    help="print the columns in the file and exit")
    args = ap.parse_args()

    if args.input.endswith((".xlsx", ".xls")):
        df = pd.read_excel(args.input, dtype=str)
    else:
        df = pd.read_csv(args.input, dtype=str, sep=None, engine="python")

    if args.inspect:
        print(f"{args.input}: {len(df)} rows")
        for c in df.columns:
            print(f"  - {c}")
        return

    mapping = {}
    missing = []
    for canonical, candidates in COLUMN_MAP_XLSX.items():
        col = pick(df, candidates)
        if col is None:
            missing.append((canonical, candidates))
        else:
            mapping[canonical] = col
    if missing:
        print("ERROR: could not map these required columns. Available columns:", file=sys.stderr)
        for c in df.columns:
            print(f"  - {c}", file=sys.stderr)
        for canonical, candidates in missing:
            print(f"  missing {canonical}: tried {candidates}", file=sys.stderr)
        print("Update COLUMN_MAP_XLSX in normalise_iedb.py and re-run.", file=sys.stderr)
        sys.exit(2)

    out = pd.DataFrame({
        "receptor_group_id": df[mapping["receptor_group_id"]],
        "cdr3_alpha": df[mapping["cdr3_alpha"]].fillna("") if "cdr3_alpha" in mapping else "",
        "cdr3_beta": df[mapping["cdr3_beta"]].fillna("") if "cdr3_beta" in mapping else "",
        "v_gene_alpha": "",
        "v_gene_beta": "",
        "epitope_sequence": df[mapping["epitope_sequence"]],
        "epitope_id": df[mapping["epitope_id"]],
        "mhc_allele": df[mapping["mhc_allele"]].fillna("") if "mhc_allele" in mapping else "",
        "species": df[mapping["species"]],
        "receptor_type": df[mapping["receptor_type"]],
        "assay_id": df[mapping["assay_id"]] if "assay_id" in mapping else "",
        "reference_id": df[mapping["reference_id"]] if "reference_id" in mapping else "",
        "assay_type": df[mapping["assay_type"]].fillna("") if "assay_type" in mapping else "",
    })
    outcome = df[mapping["assay_outcome"]].astype(str).str.lower()
    out["label"] = outcome.map({"positive": 1, "negative": 0})
    dropped = out["label"].isna().sum()
    out = out.dropna(subset=["label"])
    out["label"] = out["label"].astype(int)
    out["is_10x"] = out["assay_type"].str.contains("10X", case=False, na=False)
    out.to_csv(args.out, index=False)
    print(f"wrote {args.out}: {len(out)} rows ({dropped} dropped for non Positive/Negative outcome), "
          f"{out['label'].mean():.1%} positive, {out['epitope_id'].nunique()} epitopes")


if __name__ == "__main__":
    main()
