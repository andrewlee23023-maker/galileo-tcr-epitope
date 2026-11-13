# Canonical schema

All downstream scripts (`build_splits.py`, `baseline.py`, `evaluate.py`) read this schema. `normalise_iedb.py` is the only place that maps IEDB exports to it.

## Columns (one row = one TCR–epitope–MHC assay outcome)

| Column | Type | Meaning |
|---|---|---|
| `receptor_group_id` | str/int | IEDB receptor group ID (deduplication key). |
| `cdr3_alpha` | str | CDR3α amino-acid sequence, or empty if β-only. |
| `cdr3_beta` | str | CDR3β amino-acid sequence, or empty if α-only. |
| `v_gene_alpha` | str | e.g. TRAV1-2 (may be empty). |
| `v_gene_beta` | str | e.g. TRBV20-1 (may be empty). |
| `epitope_sequence` | str | Peptide epitope, amino-acid sequence. |
| `epitope_id` | str/int | IEDB epitope ID. |
| `mhc_allele` | str | e.g. HLA-A*02:01 (may be empty). |
| `species` | str | e.g. Homo sapiens. |
| `receptor_type` | str | e.g. "Alpha beta", "Beta", "Alpha". |
| `label` | int | 1 = binding (Positive assay), 0 = nonbinding (Negative assay). |
| `assay_id` | str/int | IEDB assay ID (may be empty). |
| `reference_id` | str/int | Source publication/reference (may be empty). |
| `assay_type` | str | e.g. tetramer, multimer, 10X single-cell — used for the 10X flag. |
| `is_10x` | bool | True if the assay came from a 10X Genomics single-cell subset. |

## IEDB mapping notes (to be confirmed with `--inspect` on the real export)

- Full export: `tcr_full_v3.xlsx` from https://www.iedb.org/database_export_v3.php (weekly builds). AIRR full export: `airr_full.zip` from the same page.
- AIRR format: one receptor chain per row; `cell_id` reconnects the α and β rows for the same receptor–epitope pair. Custom fields are prefixed `iedb_` and carry epitope, MHC and assay information. AIRR exports Calculated data fields only.
- Identifiers are Compact URIs (CURIEs): `iedb.org/receptor/<group>` for receptor groups, `iedb.org/epitope/<id>` for epitopes.
- Example from the IEDB paper (Fig. 6): receptor group 802 has a positive assay for epitope ALEPGPVTA and a negative assay for YLEPAPVTA — the same receptor, two epitopes, opposite labels. The pipeline must keep both rows until dedup, and must not let one leak into train while the other is in test for the same epitope (the epitope split prevents this; the CDR3 decontamination is the second line of defence).
- **Do not hard-code column names from memory.** Run `normalise_iedb.py --inspect` on the real file, compare with the mapping table in that script, and update the table if IEDB renamed a column. The script exits with a clear error listing the available columns if a required mapping is missing.
