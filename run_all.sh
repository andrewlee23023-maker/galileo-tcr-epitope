#!/bin/bash
# Sample end-to-end run: no IEDB download needed.
#   make sample -> splits -> baselines -> evaluate -> leakage demo
set -e
cd "$(dirname "$0")"
PY="${PY:-python3}"
$PY make_sample.py --out sample_data.csv
$PY build_splits.py sample_data.csv --out data/splits
$PY baseline.py --splits data/splits --model logreg --predict-test --out runs
$PY baseline.py --splits data/splits --model knn --predict-test --out runs
$PY evaluate.py --splits data/splits --pred runs/test_predictions_logreg.csv
$PY evaluate.py --splits data/splits --pred runs/test_predictions_knn.csv
$PY leakage_demo.py sample_data.csv
echo "OK — sample run complete. See runs/ and data/splits/DATA_CARD.md"
