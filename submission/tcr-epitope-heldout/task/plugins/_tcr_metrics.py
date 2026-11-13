"""Scores for TCR-epitope binding prediction.

Primary (declared by ``declared/metric_macro_auprc.json``):
  ``macro_auprc`` — average precision, macro-averaged across the held-out
  epitopes in the validation file (each epitope weighted equally), higher is
  better. AUPRC, not AUROC alone, because positives are rare per epitope and
  AUROC flatters a model that ranks the easy negatives correctly.

Secondaries (observational only; they influence nothing):
  ``pooled_auprc`` — average precision over all validation rows pooled.
  ``macro_auroc`` — AUROC, macro-averaged across held-out epitopes.

All three read probabilities from the deliverable
(``{row_id: {"p_binding": float}}``) and truth labels + epitopes from the
bound task data path's ``truth_by_sample_id``, so the metric and the
datasets parse the file with one set of rules and select the same rows.
A row with no prediction counts as p=0.0 (wrong, not ignored); a prediction
for an ID outside the selected rows is ignored. Epitopes with only one class
in the selected rows are left out of the macro mean (they carry no ranking
signal) — with the shipped validation file both epitopes have both classes.

Average precision and AUROC are implemented here in numpy (no scikit-learn
dependency in the SIDERIUS runtime): AP is the standard step-function
definition, AUROC is the Mann-Whitney U form with tie handling.

Declared ids must match ``IMPLEMENTS`` exactly; the framework refuses a
metric whose implementation rewrites its declared id.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np

from execute_tools.evaluation_metric import EvaluationMetric
from execute_tools.task_data_path import require_bound_task_data_path


def _average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    """Step-function average precision; labels are 0/1, higher score = more binding."""
    order = np.argsort(-scores, kind="stable")
    ranked = labels[order]
    n_pos = int(ranked.sum())
    if n_pos == 0:
        raise ValueError("average precision needs at least one positive")
    hits = np.cumsum(ranked)
    precisions = hits / np.arange(1, len(ranked) + 1)
    return float((precisions[ranked == 1]).sum() / n_pos)


def _auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    """Mann-Whitney U AUROC with average ranks for ties."""
    order = np.argsort(scores, kind="stable")
    ranks = np.empty(len(scores), dtype=float)
    sorted_scores = scores[order]
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUROC needs both classes")
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def _truth_and_scores(
    evaluation_payload: Mapping[str, Mapping[str, Any]],
    task_scope: Any,
    data_dir: str | None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Aligned (labels, p_binding, epitope) for the selected validation rows."""
    if data_dir is None:
        raise ValueError("TCR metrics read labels from the data folder; data_dir is None")
    truth_fn = getattr(require_bound_task_data_path(), "truth_by_sample_id", None)
    if not callable(truth_fn):
        raise TypeError("the bound task data path does not provide truth_by_sample_id")
    truth = truth_fn(task_scope, data_dir)
    if not isinstance(truth, dict) or not truth:
        raise ValueError("TCR metrics need at least one validation row")
    labels, scores, epitopes = [], [], []
    for sample_id in sorted(truth):
        label, epitope = truth[sample_id]
        entry = evaluation_payload.get(sample_id)
        p = float(entry["p_binding"]) if isinstance(entry, Mapping) and "p_binding" in entry else 0.0
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p_binding {p} for {sample_id} is outside [0, 1]")
        labels.append(label)
        scores.append(p)
        epitopes.append(epitope)
    return np.asarray(labels, dtype=np.int64), np.asarray(scores, dtype=float), epitopes


def _macro(labels: np.ndarray, scores: np.ndarray, epitopes: list[str], fn) -> float:
    """Mean of ``fn`` over epitopes that contain both classes."""
    values = []
    for epitope in sorted(set(epitopes)):
        mask = np.asarray([e == epitope for e in epitopes])
        if len(set(labels[mask].tolist())) < 2:
            continue
        values.append(fn(labels[mask], scores[mask]))
    if not values:
        raise ValueError("macro metric needs at least one epitope with both classes")
    return float(np.mean(values))


class MacroAuprcMetric(EvaluationMetric):
    """Primary: mean per-epitope average precision on the validation file."""

    IMPLEMENTS: ClassVar[tuple[str, ...]] = ("macro_auprc",)

    def _compute(
        self,
        deliverables: Mapping[int, str],
        /,
        *,
        evaluation_payload: Mapping[str, Mapping[str, Any]],
        task_scope: Any,
        data_dir: str | None = None,
    ) -> tuple[float, list[float | None] | None, tuple[str, ...]]:
        labels, scores, epitopes = _truth_and_scores(evaluation_payload, task_scope, data_dir)
        return _macro(labels, scores, epitopes, _average_precision), None, ()


class PooledAuprcMetric(EvaluationMetric):
    """Secondary: average precision over all validation rows pooled."""

    IMPLEMENTS: ClassVar[tuple[str, ...]] = ("pooled_auprc",)

    def _compute(
        self,
        deliverables: Mapping[int, str],
        /,
        *,
        evaluation_payload: Mapping[str, Mapping[str, Any]],
        task_scope: Any,
        data_dir: str | None = None,
    ) -> tuple[float, list[float | None] | None, tuple[str, ...]]:
        labels, scores, _ = _truth_and_scores(evaluation_payload, task_scope, data_dir)
        return _average_precision(labels, scores), None, ()


class MacroAurocMetric(EvaluationMetric):
    """Secondary: mean per-epitope AUROC on the validation file."""

    IMPLEMENTS: ClassVar[tuple[str, ...]] = ("macro_auroc",)

    def _compute(
        self,
        deliverables: Mapping[int, str],
        /,
        *,
        evaluation_payload: Mapping[str, Mapping[str, Any]],
        task_scope: Any,
        data_dir: str | None = None,
    ) -> tuple[float, list[float | None] | None, tuple[str, ...]]:
        labels, scores, epitopes = _truth_and_scores(evaluation_payload, task_scope, data_dir)
        return _macro(labels, scores, epitopes, _auroc), None, ()
