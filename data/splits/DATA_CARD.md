# Data card — TCR-epitope, epitope-held-out

Input rows: 687
After filtering/dedup: 636
Dropped non-human / non-alpha-beta: 10
Dropped missing CDR3/epitope: 0
Excluded 10X-flagged: 5
Dropped exact duplicates: 36
Dropped test rows by CDR3 decontamination: 0
Train epitopes (8): E01, E02, E03, E04, E05, E06, E07, E08
Validation epitopes (2): E09, E10
Test epitopes (2): E11, E12
Train rows: 424 (50.9% positive)
Validation rows: 106 (50.9% positive)
Test rows: 106 (50.9% positive)
Cap per epitope in train: 200
CDR3 decontamination: Hamming <= 1 or 3-mer Jaccard >= 0.9

## Pairs per epitope (after filtering)
- E01: 53
- E02: 53
- E03: 53
- E04: 53
- E05: 53
- E06: 53
- E07: 53
- E08: 53
- E09: 53
- E10: 53
- E11: 53
- E12: 53

## Pairs per MHC allele (after filtering)
- HLA-A*02:01: 318
- HLA-A*03:01: 106
- HLA-B*07:02: 106
- HLA-B*08:01: 106

Licence: IEDB data, CC BY 4.0 — cite IEDB. Sample data (if used) is synthetic.