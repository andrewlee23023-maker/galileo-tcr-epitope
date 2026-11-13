# Citation and licence

## Data — IEDB

All task data come from the Immune Epitope Database (IEDB), https://www.iedb.org.

> Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB):
> 2018 update. *Nucleic Acids Research* 2019;47(D1):D339–D343.
> https://doi.org/10.1093/nar/gky1006

> Scheffer L, et al. The Immune Epitope Database: TCR and BCR sequence data.
> bioRxiv 2026-06-06. (Scope figures used in the data card: 195,808 unique
> αβ TCR–epitope pairs, 181,168 TCRs, 2,849 peptidic epitopes, data as of
> 2026-05-28.)

IEDB data are freely available for non-commercial and commercial use under
**CC BY 4.0** since May 2017, with citation of the IEDB requested. This
package redistributes a filtered subset (human TRA/TRB, peptidic epitopes
8–15 AA, deduplicated, capped) plus synthetic negatives derived from it;
the IEDB attribution above must be kept with any redistribution.

## Software

- SIDERIUS (task framework): https://github.com/yuema137/SIDERIUS — see its
  own LICENSE. This package adapts its CSV classification template
  (`examples/task_templates/csv_classification/`).
- Baselines in this package were written for this task (3-mer logistic
  regression and k-NN, in the style of TITAN's `knn_cv.py`); no TITAN or
  ERGO-II code or pretrained weights are shipped, and no TITAN Box data are
  used. TITAN: https://github.com/PaccMann/TITAN (MIT). ERGO-II:
  https://github.com/IdoSpringer/ERGO-II (MIT).

## Task package

TCR–epitope binding prediction, epitope-held-out. Author: Andrew Lee,
University of Pennsylvania. Built for the Galileo Hackathon, Track 1, Nov
2026. If this task is made public for Track 2, it is published under the
MIT License, only with the author's agreement (see the hackathon page).
