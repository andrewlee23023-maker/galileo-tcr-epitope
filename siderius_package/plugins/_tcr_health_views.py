"""Health view: the predicted classes as a categorical stream.

Listed under ``plugins:`` in ``declared/task_health.yaml`` and loaded by the
Health owner, which calls ``register_view_provider`` below. The framework's
generic classification checks (``categorical_distinct_symbols``,
``categorical_dominant_fraction``) read this view; this file only turns the
deliverable's probabilities into class indices (p_binding >= 0.5 -> 1).
"""

from __future__ import annotations

from typing import Any, ClassVar

import numpy as np

from execute_tools.health_checks import (
    CATEGORICAL_PREDICTIONS,
    CategoricalPredictionsPayload,
    HealthView,
    register_view_provider,
)
from execute_tools.health_checks.schemas import HealthCheckContext

#: Must equal ``providers[].provider_id`` in ``declared/task_health.yaml``.
PROVIDER_ID = "tcr_epitope.predicted_classes"


class PredictedClassesView:
    provider_id: ClassVar[str] = PROVIDER_ID
    capabilities: ClassVar[frozenset[str]] = frozenset({CATEGORICAL_PREDICTIONS})

    def materialize(
        self,
        capability_key: str,
        ctx: HealthCheckContext,
        config: dict[str, Any] | None = None,
    ) -> HealthView:
        decoded = ctx.load_evaluation_payload()
        if not isinstance(decoded, dict):
            raise TypeError(f"predictions must be a mapping, got {type(decoded).__name__}")
        symbols = np.asarray(
            [1 if float(decoded[key]["p_binding"]) >= 0.5 else 0 for key in sorted(decoded)],
            dtype=np.int64,
        )
        return HealthView(
            capability_key=capability_key,
            provider_id=self.provider_id,
            payload=CategoricalPredictionsPayload(symbols=symbols),
        )


register_view_provider(PredictedClassesView())
