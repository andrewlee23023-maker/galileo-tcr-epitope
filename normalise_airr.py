#!/usr/bin/env python3
"""Normalise the IEDB AIRR TCR export to the canonical schema (see SCHEMA.md).

Input: tcr_full_airr.tsv from airr_full.zip
  (https://www.iedb.org/database_export_v3.php -> airr_full.zip,
   direct link: https://www.iedb.org/downloader.php?file_name=doc/airr_full.zip)

AIRR is one chain per row; cell_id reconnects the alpha and beta rows for
the same receptor-epitope pair. Columns used:
  cell_id, iedb_receptor_group, reactivity_ref, iedb_epitope_name,
  iedb_assay_outcome (positive/negative), iedb_mhc_alleles,
  locus (TRA/TRB), locus_species, v_call, cdr3_aa

Verified on the Nov 3, 2026 build (downloaded Nov 8, 2026):
  339,842 chain rows; 302,511 cell_ids; 3,349 epitopes.
  Outcomes: positive 338,432 (99.6%), negative 1,410 (0.4%).
  Human (NCBITAXON:9606): 327,682 rows. Paired TRA+TRB: only 33,621
  cell_ids; most receptors are TRB-only (198,504) or TRA-only (59,161).
  Top epitope LLLLDRLNQLESKMS alone has 60,677 rows.

Consequences (do not hide these):
  1. Real negatives are far too few for a binary task (1,232 in the
     human, pure-peptide 8-15 subset). Synthetic negatives are therefore
     REQUIRED, and must be labelled as synthetic: each positive pair gets
     --synthetic-negatives cross-paired TCRs with a different epitope it
     was never reported to bind. Synthetic negatives are generated AFTER
     the epitope split is fixed in build_splits.py? No — they are generated
     here, but only within the same split: build_splits.py splits by
     epitope, so a synthetic negative whose epitope is a test epitope
     lands in test, and one whose epitope is a train epitope lands in
     train. A synthetic negative never leaks a test epitope into train,
     because its epitope is the split key. The TCR in a synthetic negative
     may appear in train with its true epitope — that is the intended
     task (same TCR, different epitope, different label).
  2. The task is therefore primarily TRB-driven: receptor_type is
     "Alpha beta" only where both chains are present, else "Beta" or
     "Alpha". Baselines use CDR3-beta where present.
  3. Epitope skew is extreme. --cap-per-epitope (default 1000) caps
     positives per epitope before synthetic negatives are made.

Filters:
  - human only (NCBITAXON:9606), TRA/TRB only, CDR3 present
  - epitope: pure amino acids, length 8-15
  - MHC: first specific allele (containing '*') if present, else the
    generic class label ("HLA class I" / "HLA class II"), else blank
  - dedup by (receptor_group, cdr3_alpha, cdr3_beta, epitope, mhc)
"""
import argparse
import json
import random
import re

import pandas as pd

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{8,15}$")


def first_specific_mhc(s):
    if not isinstance(s, str) or not s:
        return ""
    parts = [p.strip() for p in s.split("|") if p.strip()]
    for p in parts:
        if "*" in p:
            return p
    return parts[0] if parts else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="tcr_full_airr.tsv")
    ap.add_argument("--out", default="data/normalised.csv")
    ap.add_argument("--cap-per-epitope", type=int, default=1000)
    ap.add_argument("--synthetic-negatives", type=int, default=1,
                    help="synthetic negatives per positive (0 to disable)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    cols = ["cell_id", "iedb_receptor_group", "reactivity_ref", "iedb_epitope_name",
            "iedb_assay_outcome", "iedb_mhc_alleles", "locus", "locus_species",
            "v_call", "cdr3_aa"]
    df = pd.read_csv(args.input, sep="\t", usecols=cols, dtype=str)
    n_in = len(df)
    df = df[df["locus_species"] == "NCBITAXON:9606"]
    df = df[df["locus"].isin(["TRA", "TRB"])]
    df = df[df["cdr3_aa"].notna() & (df["cdr3_aa"] != "")]
    df = df[df["iedb_epitope_name"].str.fullmatch(AA_RE, na=False)]
    n_filtered = len(df)

    # one pair row per (cell_id, epitope, outcome): combine TRA + TRB
    rows = []
    for (cell, epitope, outcome), g in df.groupby(
            ["cell_id", "iedb_epitope_name", "iedb_assay_outcome"], sort=False):
        tra = g[g["locus"] == "TRA"]
        trb = g[g["locus"] == "TRB"]
        first = g.iloc[0]
        rows.append({
            "receptor_group_id": first["iedb_receptor_group"],
            "cdr3_alpha": tra["cdr3_aa"].iloc[0] if len(tra) else "",
            "cdr3_beta": trb["cdr3_aa"].iloc[0] if len(trb) else "",
            "v_gene_alpha": tra["v_call"].iloc[0] if len(tra) else "",
            "v_gene_beta": trb["v_call"].iloc[0] if len(trb) else "",
            "epitope_sequence": epitope,
            "epitope_id": epitope,
            "mhc_allele": first_specific_mhc(first["iedb_mhc_alleles"]),
            "species": "Homo sapiens",
            "receptor_type": ("Alpha beta" if len(tra) and len(trb)
                              else "Beta" if len(trb) else "Alpha"),
            "label": 1 if outcome == "positive" else 0,
            "assay_id": "",
            "reference_id": first["iedb_receptor_group"],
            "assay_type": "",
            "is_10x": False,
            "iedb_epitope_ref": first["reactivity_ref"],
        })
    pairs = pd.DataFrame(rows)
    n_pairs = len(pairs)
    pairs = pairs.drop_duplicates(
        subset=["receptor_group_id", "cdr3_alpha", "cdr3_beta", "epitope_id", "mhc_allele"])
    n_dedup = len(pairs)

    # cap positives per epitope (keeps the mega-epitopes from dominating)
    pos = pairs[pairs["label"] == 1]
    neg_real = pairs[pairs["label"] == 0]
    capped = []
    for _, g in pos.groupby("epitope_id"):
        capped.append(g.sample(n=min(len(g), args.cap_per_epitope), random_state=args.seed))
    pos_capped = pd.concat(capped) if capped else pos.iloc[0:0]

    # synthetic negatives: PER-EPITOPE balanced. The first version of this
    # script sampled a random epitope per positive TCR, which left each
    # epitope with ~35 synthetic negatives against ~500-1000 positives
    # (test prevalence ~94% — AUPRC would be meaningless). Instead, for
    # each epitope E, sample as many TCRs from OTHER epitopes as E has
    # positives, and assign them to E as synthetic negatives. That gives
    # ~1:1 per epitope in every split (splits are by epitope, so a
    # synthetic negative for E lands in the same split as E's positives).
    # A sampled TCR is skipped if that TCR-epitope pair is a known positive.
    known = set(zip(pairs["cdr3_alpha"], pairs["cdr3_beta"], pairs["epitope_id"]))
    syn_rows = []
    if args.synthetic_negatives > 0:
        # pool of donor TCRs: one row per distinct (cdr3a, cdr3b) with its true epitope
        donors = pos_capped.drop_duplicates(subset=["cdr3_alpha", "cdr3_beta"])
        donor_list = list(donors.itertuples(index=False))
        donor_by_epitope = {}
        for d in donor_list:
            donor_by_epitope.setdefault(d.epitope_id, []).append(d)
        for ep, g in pos_capped.groupby("epitope_id"):
            n_needed = int(len(g) * args.synthetic_negatives)
            # donor TCRs whose true epitope is not ep
            pool = [d for d in donor_list if d.epitope_id != ep]
            rng.shuffle(pool)
            made = 0
            for d in pool:
                if made >= n_needed:
                    break
                if (d.cdr3_alpha, d.cdr3_beta, ep) in known:
                    continue
                syn_rows.append({
                    "receptor_group_id": d.receptor_group_id,
                    "cdr3_alpha": d.cdr3_alpha,
                    "cdr3_beta": d.cdr3_beta,
                    "v_gene_alpha": d.v_gene_alpha,
                    "v_gene_beta": d.v_gene_beta,
                    "epitope_sequence": ep,
                    "epitope_id": ep,
                    "mhc_allele": g["mhc_allele"].iloc[0],
                    "species": "Homo sapiens",
                    "receptor_type": d.receptor_type,
                    "label": 0,
                    "assay_id": "",
                    "reference_id": d.receptor_group_id,
                    "assay_type": "",
                    "is_10x": False,
                    "iedb_epitope_ref": "",
                    "is_synthetic_negative": True,
                })
                made += 1
    syn = pd.DataFrame(syn_rows) if syn_rows else pd.DataFrame(columns=pos_capped.columns)
    out = pd.concat([pos_capped, neg_real, syn], ignore_index=True)
    if "is_synthetic_negative" not in out.columns:
        out["is_synthetic_negative"] = False
    out["is_synthetic_negative"] = out.get("is_synthetic_negative", False).fillna(False).astype(bool)
    # real negatives have no synthetic flag
    out.loc[out["label"] == 1, "is_synthetic_negative"] = False
    out = out.sample(frac=1, random_state=args.seed)
    out.to_csv(args.out, index=False)

    summary = {
        "input_chain_rows": int(n_in),
        "after_human_tra_trb_cdr3_peptide_filter": int(n_filtered),
        "pair_rows_before_dedup": int(n_pairs),
        "pair_rows_after_dedup": int(n_dedup),
        "real_positives_after_cap": int(len(pos_capped)),
        "real_negatives_kept": int(len(neg_real)),
        "synthetic_negatives": int(len(syn)),
        "total_rows": int(len(out)),
        "distinct_epitopes": int(out["epitope_id"].nunique()),
        "positive_rate": float(out["label"].mean()),
        "cap_per_epitope": args.cap_per_epitope,
    }
    print(json.dumps(summary, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
