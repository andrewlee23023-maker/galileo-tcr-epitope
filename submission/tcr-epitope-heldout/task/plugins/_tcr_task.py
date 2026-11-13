"""Task data path for TCR-epitope binding prediction from IEDB data.

Adapted from the canonical SIDERIUS CSV classification template
(``examples/task_templates/csv_classification/plugins/_csv_task.py`` in
``yuema137/SIDERIUS``): same scope / deliverable / health-coverage structure,
different features. Each row is one TCR-epitope pair; the model sees

  * 3-mer counts (min count 2 in train+validation) over CDR3-alpha,
    CDR3-beta and the epitope peptide — 7,035 features, and
  * a one-hot MHC allele (145 alleles seen in train+validation; an unseen
    allele maps to the all-zero vector) — 145 features,

for 7,180 float32 features total. The vocabulary lives in
``declared/kmer_vocab.json`` and is passed via ``task.yaml``
``task_data_path.config.vocab_file`` (a ``{ref:}`` path, resolved against the
manifest directory), so the package is relocatable and the feature space is
fixed at composition time.

The deliverable is ``{row_id: {"p_binding": float}}`` — one binding
probability per validation row. The metric reads probabilities (not hard
classes), because the primary score is macro-averaged AUPRC. Health views
threshold at 0.5 to get class indices for the generic categorical checks.

Data layout. The run's data folder holds ``train.csv`` and
``validation.csv`` (epitope-disjoint by construction; see
``SPLIT_EVIDENCE.md``). The final test file (``test_inputs.csv`` /
``test_labels.csv``) is kept OUTSIDE the run data folder: nothing in this
package reads it. ``test_labels.csv`` is evaluator-only.

Scope builders receive no data folder, so a scope records only *how* to
select rows (which file, which fraction, which seed); rows are resolved in
the stages that receive the data folder, exactly as in the template.
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from core.seeding import MissingRunSeed
from execute_tools.task_data_path import (
    DeliverableWriteRequest,
    EpochSamplingParams,
    EvalMaterializationParams,
    EvaluationReadRequest,
    HealthCoverageRequest,
    HealthCoverageResult,
    ScopeBuildRequest,
    StorageReadScope,
    TaskOutputArtifactInventory,
    ValidationScopeError,
)

if TYPE_CHECKING:
    from torch.utils.data import Dataset

TASK_DATA_PATH_ID = "tcr_epitope_heldout"

TRAINING_PARTITION = 0
VALIDATION_PARTITION = 1

SCOPE_KIND = "tcr_epitope_rows"

SEQUENCE_COLUMNS = ("cdr3_alpha", "cdr3_beta", "epitope_sequence")


class TcrRowSelection(BaseModel):
    """A scope: which CSV file to read and how to sample its rows.

    ``seed`` is always concrete, so every later stage resolves the identical
    row subset. A sampled scope (``portion < 1`` or ``max_samples``) takes the
    framework's ``request.seed``; the builder refuses a sampled request
    without one rather than inventing a seed. An unsampled scope reads every
    row and stores ``seed=0``, which no stage uses.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    kind: Literal["tcr_epitope_rows"] = SCOPE_KIND
    file: str = Field(min_length=1)
    portion: float = Field(gt=0.0, le=1.0)
    seed: int
    max_samples: int | None = Field(default=None, ge=1)


@dataclass(frozen=True)
class TcrTable:
    """One parsed CSV file: row IDs, features, labels and epitopes, row-aligned."""

    sample_ids: tuple[str, ...]
    features: np.ndarray  # float32 [rows, feature_count]
    labels: np.ndarray  # int64 [rows]
    epitopes: tuple[str, ...]


class TcrEpitopeTaskDataPath:
    """Reads labelled TCR-epitope rows and writes ``{row_id: {p_binding}}`` JSON.

    Constructor arguments come from ``task.yaml`` ``task_data_path.config``.
    Every argument is validated here, so a wrong column name or vocab file
    fails at composition or first read, with a message naming the value.
    """

    task_data_path_id: ClassVar[str] = TASK_DATA_PATH_ID

    def __init__(
        self,
        *,
        train_file: str,
        validation_file: str,
        label_column: str,
        classes: Sequence[Any],
        id_column: str | None = None,
        vocab_file: str | None = None,
    ) -> None:
        for name, value in (("train_file", train_file), ("validation_file", validation_file)):
            if not value or os.path.isabs(value) or os.path.normpath(value).startswith(".."):
                raise ValueError(
                    f"{name} must be a file name inside the data folder; got {value!r}"
                )
        if train_file == validation_file:
            raise ValueError("train_file and validation_file must be different files")
        self._train_file = train_file
        self._validation_file = validation_file
        self._label_column = label_column
        self._id_column = id_column
        self._classes = tuple(str(value) for value in classes)
        if len(self._classes) < 2 or len(set(self._classes)) != len(self._classes):
            raise ValueError(f"classes must list at least two distinct labels; got {classes!r}")
        self._class_index = {label: index for index, label in enumerate(self._classes)}
        if vocab_file is None:
            raise ValueError("vocab_file is required: path to declared/kmer_vocab.json")
        with open(vocab_file, encoding="utf-8") as handle:
            vocab_doc = json.load(handle)
        self._kmer_index = {k: i for i, k in enumerate(vocab_doc["vocab"])}
        self._mhc_index = {a: i for i, a in enumerate(vocab_doc["mhc_vocab"])}
        self._feature_count = len(self._kmer_index) + len(self._mhc_index)

    @property
    def class_count(self) -> int:
        return len(self._classes)

    @property
    def feature_count(self) -> int:
        return self._feature_count

    # -- reading the CSV files -------------------------------------------------

    def _featurise(self, record: Mapping[str, str]) -> np.ndarray:
        """3-mer counts over the three sequence columns + MHC one-hot.

        Returns a float32 array, not a Python list: a dense 7,180-wide row
        as Python floats would cost ~10x the memory and OOMs on the 61,585
        training rows (found in validation, Nov 8, 2026).
        """
        vec = np.zeros(self._feature_count, dtype=np.float32)
        for col in SEQUENCE_COLUMNS:
            seq = (record.get(col) or "").strip()
            for i in range(len(seq) - 2):
                idx = self._kmer_index.get(seq[i:i + 3])
                if idx is not None:
                    vec[idx] += 1.0
        mhc = (record.get("mhc_allele") or "").strip()
        midx = self._mhc_index.get(mhc)
        if midx is not None:
            vec[len(self._kmer_index) + midx] = 1.0
        return vec

    def read_table(self, data_dir: str, file: str) -> TcrTable:
        """Parse one CSV file, refusing missing values, unknown labels and duplicate IDs."""
        path = os.path.join(data_dir, file)
        if not os.path.isfile(path):
            raise ValidationScopeError(f"CSV file not found in the data folder: {path}")
        with open(path, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            header = tuple(reader.fieldnames or ())
            required = [self._label_column, "epitope_sequence"] + list(SEQUENCE_COLUMNS)
            if self._id_column:
                required.append(self._id_column)
            missing = [c for c in required if c not in header]
            if missing:
                raise ValueError(f"{path} has no column(s) {missing}; its header is {list(header)}")
            records = list(reader)
        if not records:
            raise ValueError(f"{path} has a header but no rows")
        # Preallocate the feature matrix once and fill it in place: holding
        # one float32 array per row and then stacking doubles peak memory
        # (61,585 x 7,180), which OOMs a modest machine (found Nov 8, 2026).
        features = np.zeros((len(records), self._feature_count), dtype=np.float32)
        ids: list[str] = []
        labels: list[int] = []
        epitopes: list[str] = []
        for line, record in enumerate(records, start=2):
            label = (record.get(self._label_column) or "").strip()
            if label not in self._class_index:
                raise ValueError(
                    f"{path} line {line}: label {label!r} is not one of classes "
                    f"{list(self._classes)}"
                )
            sample_id = (
                f"row{line - 2:06d}"
                if self._id_column is None
                else (record.get(self._id_column) or "").strip()
            )
            if not sample_id:
                raise ValueError(f"{path} line {line}: empty {self._id_column!r}")
            ids.append(sample_id)
            features[line - 2] = self._featurise(record)
            labels.append(self._class_index[label])
            epitopes.append((record.get("epitope_sequence") or "").strip())
        if len(set(ids)) != len(ids):
            raise ValueError(f"{path} repeats a value in the ID column {self._id_column!r}")
        return TcrTable(
            sample_ids=tuple(ids),
            features=features,
            labels=np.asarray(labels, dtype=np.int64),
            epitopes=tuple(epitopes),
        )

    @staticmethod
    def _subsample(
        count: int, portion: float | None, max_samples: int | None, seed: int | None
    ) -> list[int]:
        """Indices ``0..count-1`` kept by a portion and a ceiling, in file order.

        The kept rows are a seeded random subset, never the first rows of the
        file, so a ceiling on a label-sorted CSV still sees every class.
        """
        keep = count
        if portion is not None and portion < 1.0:
            keep = max(1, math.floor(count * portion))
        if max_samples is not None:
            keep = min(keep, max_samples)
        if keep >= count:
            return list(range(count))
        if seed is None:
            raise MissingRunSeed("subsampling CSV rows needs a seed; none was supplied")
        drawn = np.random.default_rng(seed).permutation(count)[:keep]
        return sorted(int(index) for index in drawn)

    def _select(self, row_count: int, scope: TcrRowSelection) -> list[int]:
        return self._subsample(row_count, scope.portion, scope.max_samples, scope.seed)

    def _rows(self, scope: object, data_dir: str) -> tuple[TcrTable, list[int]]:
        selection = self._require_scope(scope)
        table = self.read_table(data_dir, selection.file)
        return table, self._select(len(table.sample_ids), selection)

    def labels_by_sample_id(self, scope: object, data_dir: str) -> dict[str, int]:
        """Ground-truth class index per selected sample; the metric's truth source."""
        table, positions = self._rows(scope, data_dir)
        return {table.sample_ids[p]: int(table.labels[p]) for p in positions}

    def truth_by_sample_id(self, scope: object, data_dir: str) -> dict[str, tuple[int, str]]:
        """Ground truth per selected sample: (class index, epitope sequence).

        The macro-AUPRC metric needs the epitope per row, not just the label,
        so it can average within each held-out epitope. Labels and epitopes
        come from the same parse as the datasets, so no stage can disagree
        about them.
        """
        table, positions = self._rows(scope, data_dir)
        return {
            table.sample_ids[p]: (int(table.labels[p]), table.epitopes[p])
            for p in positions
        }

    # -- datasets ----------------------------------------------------------------

    @staticmethod
    def _dataset(features: np.ndarray, labels: np.ndarray) -> Dataset[Any]:
        import torch
        from torch.utils.data import TensorDataset

        return TensorDataset(torch.from_numpy(features.copy()), torch.from_numpy(labels.copy()))

    def training_dataset(self, scope: object, params: EpochSamplingParams) -> Dataset[Any]:
        table, positions = self._rows(scope, params.data_dir)
        kept = self._subsample(
            len(positions), params.train_portion, params.max_samples, params.epoch_seed
        )
        positions = [positions[index] for index in kept]
        return self._dataset(table.features[positions], table.labels[positions])

    def validation_dataset(self, scope: object, params: EvalMaterializationParams) -> Dataset[Any]:
        table, positions = self._rows(scope, params.data_dir)
        return self._dataset(table.features[positions], table.labels[positions])

    # -- deliverable: {row_id: {"p_binding": float}} as one JSON file --------------

    @staticmethod
    def deliverable_name(*, model_type: str, run_name: str, exp_id: str) -> str:
        """The file name from the manifest's ``deliverable:`` section."""
        return _deliverable_name(model_type=model_type, run_name=run_name, exp_id=exp_id)

    def write_deliverable(self, outputs: Any, request: DeliverableWriteRequest) -> None:
        """Pair each output with its row ID and write the binding probabilities.

        Generic inference passes one output per validation row, in
        ``validation_dataset`` order, plus the scope and the data folder in
        ``request.source_context``. A ready ``{row_id: {"p_binding": ...}}``
        mapping is also accepted. Logits become probabilities by softmax over
        the two class logits; a bare probability is accepted as-is.
        """
        if isinstance(outputs, Mapping):
            predictions = {str(key): self._prob_of(value) for key, value in outputs.items()}
        else:
            if request.task_scope is None or request.source_context is None:
                raise ValueError(
                    "positional outputs need request.task_scope and request.source_context"
                )
            table, positions = self._rows(request.task_scope, request.source_context.data_dir)
            sequence = list(outputs)
            if len(sequence) != len(positions) or len(positions) != (
                request.source_context.sample_count
            ):
                raise ValueError(
                    f"{len(sequence)} outputs for {len(positions)} validation rows "
                    f"(inference reported {request.source_context.sample_count}); refusing to pair"
                )
            predictions = {
                table.sample_ids[p]: self._prob_of(output)
                for p, output in zip(positions, sequence, strict=True)
            }
        os.makedirs(request.output_dir, exist_ok=True)
        name = self.deliverable_name(
            model_type=request.model_type, run_name=request.run_name, exp_id=request.exp_id
        )
        with open(os.path.join(request.output_dir, name), "w", encoding="utf-8") as handle:
            json.dump(
                {k: {"p_binding": v} for k, v in predictions.items()},
                handle, sort_keys=True, indent=2,
            )

    def _prob_of(self, prediction: object) -> float:
        """Binding probability from model logits (softmax) or a bare probability."""
        import torch

        if isinstance(prediction, torch.Tensor):
            if prediction.ndim == 0:
                value = float(prediction.item())
            else:
                import torch.nn.functional as F

                value = float(F.softmax(prediction.float(), dim=-1)[1].item())
        elif isinstance(prediction, Mapping) and "p_binding" in prediction:
            value = float(prediction["p_binding"])  # type: ignore[arg-type]
        elif isinstance(prediction, (int, float, np.floating)) and not isinstance(prediction, bool):
            value = float(prediction)
        else:
            raise ValueError(f"cannot read a binding probability from {type(prediction).__name__}")
        if not 0.0 <= value <= 1.0 or not math.isfinite(value):
            raise ValueError(f"binding probability {value} is outside [0, 1]")
        return value

    def read_evaluation_payload(self, request: EvaluationReadRequest) -> object:
        """Decode the deliverable back into ``{row_id: {"p_binding": float}}``."""
        name = self.deliverable_name(
            model_type=request.model_type, run_name=request.run_name, exp_id=request.exp_id
        )
        path = os.path.join(request.deliverable_dir, name)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"prediction file not found: {path}")
        with open(path, encoding="utf-8") as handle:
            decoded = json.load(handle)
        if not isinstance(decoded, dict):
            raise ValueError(f"{path} must hold a JSON object")
        return {
            str(key): {"p_binding": float(value["p_binding"])}
            for key, value in decoded.items()
        }

    def enumerate_output_artifacts(
        self, request: EvaluationReadRequest
    ) -> TaskOutputArtifactInventory:
        name = self.deliverable_name(
            model_type=request.model_type, run_name=request.run_name, exp_id=request.exp_id
        )
        path = Path(request.deliverable_dir) / name
        return TaskOutputArtifactInventory(
            run_name=request.run_name,
            exp_id=request.exp_id,
            model_type=request.model_type,
            relative_paths=(name,) if path.exists() or path.is_symlink() else (),
        )

    # -- scopes ------------------------------------------------------------------

    def _require_scope(self, scope: object) -> TcrRowSelection:
        if not isinstance(scope, TcrRowSelection):
            raise ValueError(f"expected a TcrRowSelection scope, got {type(scope).__name__}")
        return scope

    def _selection(self, file: str, request: ScopeBuildRequest) -> TcrRowSelection:
        if request.subset_ref is not None:
            raise ValueError(f"{TASK_DATA_PATH_ID} has no data subsets; got {request.subset_ref!r}")
        seed = request.seed
        if seed is None:
            if request.portion < 1.0 or request.max_samples is not None:
                raise MissingRunSeed(
                    "a sampled CSV scope (portion < 1 or max_samples) requires request.seed"
                )
            seed = 0  # every row is read; the seed is never used
        return TcrRowSelection(
            file=file, portion=request.portion, seed=seed, max_samples=request.max_samples
        )

    def build_training_scope(self, request: ScopeBuildRequest) -> object:
        if request.selection_strategy == "target" and tuple(request.target_partitions) != (
            TRAINING_PARTITION,
        ):
            raise ValueError(
                f"the only training partition is {TRAINING_PARTITION}; "
                f"got target partitions {list(request.target_partitions)}"
            )
        return self._selection(self._train_file, request)

    def build_eval_scope(self, request: ScopeBuildRequest) -> object:
        return self._selection(self._validation_file, request)

    def serialize_scope(self, scope: object) -> str:
        return json.dumps(
            self._require_scope(scope).model_dump(mode="json"),
            sort_keys=True,
            separators=(",", ":"),
        )

    def deserialize_scope(self, payload: str) -> object:
        return TcrRowSelection.model_validate_json(payload)

    # -- optional capabilities the console and Health use --------------------------

    def storage_read_scope(self, data_dir: str, scope: object) -> StorageReadScope:
        """The one CSV file a scope reads, for the Task tab's data volume."""
        path = os.path.abspath(os.path.join(data_dir, self._require_scope(scope).file))
        return StorageReadScope(file_paths=(path,), expected_on_disk_bytes=os.path.getsize(path))

    def validate_health_coverage(self, request: HealthCoverageRequest) -> HealthCoverageResult:
        """Health reads the predictions for the validation file; confirm the scope is it."""
        scope = self._require_scope(request.evaluation_scope)
        override = request.health_gate_files
        if override is not None and tuple(override) != (VALIDATION_PARTITION,):
            reason = (
                f"Health reads only partition {VALIDATION_PARTITION} (the validation "
                f"file); a health_gate_files override of {list(override)} is not supported"
            )
            covered = False
        elif scope.file != self._validation_file:
            reason = f"evaluation scope reads {scope.file}, not {self._validation_file}"
            covered = False
        else:
            reason = f"Health checks read the predictions for {scope.file}"
            covered = True
        return HealthCoverageResult(applicable=True, covered=covered, reason=reason)


def _deliverable_name(*, model_type: str, run_name: str, exp_id: str) -> str:
    from execute_tools.deliverable_spec import resolve_deliverable_naming

    return resolve_deliverable_naming().name(
        model_type=model_type, run_name=run_name, exp_id=exp_id, input_identity=0
    )


def deliverable_name(request: DeliverableWriteRequest | EvaluationReadRequest) -> str:
    """Module-level form of the deliverable name, used by inference reporting."""
    return _deliverable_name(
        model_type=request.model_type, run_name=request.run_name, exp_id=request.exp_id
    )
