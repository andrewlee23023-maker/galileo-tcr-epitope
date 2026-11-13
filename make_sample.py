#!/usr/bin/env python3
"""Generate a small synthetic TCR-epitope dataset in the canonical schema.

Design (so the leakage demo shows a real gap):
- 12 epitopes, each with a private 3-mer motif. Binders for an epitope carry
  that motif in CDR3-beta; non-binders are random CDR3s.
- A model that has seen an epitope can learn its motif (random split: high
  score). A model facing an unseen epitope has never seen its motif
  (epitope-held-out: near chance). That gap is the point of the task.
- Includes deliberate duplicates (same CDR3+epitope in 2 publications),
  near-duplicate CDR3s (1 AA change), a few mouse rows, a few TRA-only rows,
  and a few 10X-flagged rows, so filtering/dedup/decontamination have
  something to do. These are counted in DATA_CARD.md, not hidden.

This sample proves the pipeline runs. Its numbers say nothing about how
solvable the real IEDB task is.
"""
import argparse
import random
import string

import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"


def rand_cdr3(rng, length=13):
    return "".join(rng.choice(AA) for _ in range(length))


def with_motif(rng, motif):
    seq = list(rand_cdr3(rng, 13))
    pos = rng.randint(2, 8)
    seq[pos:pos + 3] = list(motif)
    return "".join(seq)


EPITOPES = [
    ("E01", "ALEPGPVTA", "HLA-A*02:01"),
    ("E02", "YLEPAPVTA", "HLA-A*02:01"),
    ("E03", "GILGFVFTL", "HLA-A*02:01"),
    ("E04", "NLVPMVATV", "HLA-A*02:01"),
    ("E05", "KLGGALQAK", "HLA-A*03:01"),
    ("E06", "RLRAEAQVK", "HLA-A*03:01"),
    ("E07", "ATDALMTGF", "HLA-B*07:02"),
    ("E08", "RPRGEVRFL", "HLA-B*07:02"),
    ("E09", "TPQDLNTML", "HLA-B*08:01"),
    ("E10", "FLRGRAYGL", "HLA-B*08:01"),
    ("E11", "AAVVRFQEA", "HLA-A*02:01"),
    ("E12", "QAAVVRFQE", "HLA-A*02:01"),
]
MOTIFS = {eid: "".join(random.Random(eid).choices(AA, k=3)) for eid, _, _ in EPITOPES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="sample_data.csv")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--n-per-epitope", type=int, default=50)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    rows = []
    gid = 1000
    for eid, epitope, mhc in EPITOPES:
        motif = MOTIFS[eid]
        for i in range(args.n_per_epitope):
            is_pos = i % 2 == 0
            cdr3b = with_motif(rng, motif) if is_pos else rand_cdr3(rng)
            cdr3a = rand_cdr3(rng, 12)
            rows.append({
                "receptor_group_id": gid,
                "cdr3_alpha": cdr3a,
                "cdr3_beta": cdr3b,
                "v_gene_alpha": "TRAV1-2",
                "v_gene_beta": "TRBV20-1",
                "epitope_sequence": epitope,
                "epitope_id": eid,
                "mhc_allele": mhc,
                "species": "Homo sapiens",
                "receptor_type": "Alpha beta",
                "label": 1 if is_pos else 0,
                "assay_id": f"A{gid}",
                "reference_id": f"R{rng.randint(1, 20)}",
                "assay_type": "tetramer",
                "is_10x": False,
            })
            gid += 1
            # deliberate exact duplicate in a second publication
            if i % 17 == 0:
                dup = dict(rows[-1])
                dup["assay_id"] = f"A{gid}"
                dup["reference_id"] = f"R{rng.randint(21, 40)}"
                rows.append(dup)
                gid += 1
            # near-duplicate CDR3 (1 AA change) — decontamination target
            if i % 23 == 0:
                near = list(cdr3b)
                near[rng.randint(0, 12)] = rng.choice(AA)
                nd = dict(rows[-1])
                nd["receptor_group_id"] = gid
                nd["cdr3_beta"] = "".join(near)
                nd["assay_id"] = f"A{gid}"
                rows.append(nd)
                gid += 1

    # rows that filtering must drop / flag
    for k in range(6):
        rows.append({
            "receptor_group_id": gid, "cdr3_alpha": rand_cdr3(rng, 12),
            "cdr3_beta": rand_cdr3(rng), "v_gene_alpha": "TRAV1-2",
            "v_gene_beta": "TRBV20-1", "epitope_sequence": "GILGFVFTL",
            "epitope_id": "E03", "mhc_allele": "H-2Kb",
            "species": "Mus musculus", "receptor_type": "Alpha beta",
            "label": 1, "assay_id": f"A{gid}", "reference_id": "R99",
            "assay_type": "tetramer", "is_10x": False,
        })
        gid += 1
    for k in range(4):
        rows.append({
            "receptor_group_id": gid, "cdr3_alpha": rand_cdr3(rng, 12),
            "cdr3_beta": "", "v_gene_alpha": "TRAV1-2", "v_gene_beta": "",
            "epitope_sequence": "NLVPMVATV", "epitope_id": "E04",
            "mhc_allele": "HLA-A*02:01", "species": "Homo sapiens",
            "receptor_type": "Alpha", "label": 1, "assay_id": f"A{gid}",
            "reference_id": "R98", "assay_type": "tetramer", "is_10x": False,
        })
        gid += 1
    for k in range(5):
        rows.append({
            "receptor_group_id": gid, "cdr3_alpha": rand_cdr3(rng, 12),
            "cdr3_beta": rand_cdr3(rng), "v_gene_alpha": "TRAV1-2",
            "v_gene_beta": "TRBV20-1", "epitope_sequence": "KLGGALQAK",
            "epitope_id": "E05", "mhc_allele": "HLA-A*03:01",
            "species": "Homo sapiens", "receptor_type": "Alpha beta",
            "label": 1, "assay_id": f"A{gid}", "reference_id": "R97",
            "assay_type": "10X single-cell", "is_10x": True,
        })
        gid += 1

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out}: {len(df)} rows, {df['epitope_id'].nunique()} epitopes, "
          f"{df['label'].mean():.2%} positive")


if __name__ == "__main__":
    main()
