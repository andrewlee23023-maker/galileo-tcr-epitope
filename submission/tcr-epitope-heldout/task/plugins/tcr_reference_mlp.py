"""Reference model: Linear(features -> hidden) -> ReLU -> Linear(hidden -> 2).

``task.yaml`` ``model_plugins.require`` names this model type, so every run
has one known model that fits the task's input and output. The research loop
proposes its own models; this one only has to run, not to score well — and
on this task it should not be expected to: the shipped baselines (see
``baseline/`` and ``RUN_LOG.md`` in the parent folder) show additive 3-mer
features are at or below chance on truly unseen epitopes. A model that beats
it needs to learn TCR x epitope interaction, not additive sequence features.

Model plugin contract: the module defines ``PLUGIN_MODEL_TYPE``,
``PLUGIN_OUTPUT_TYPE``, ``PLUGIN_CONFIG_CLASS`` and ``PLUGIN_MODEL_CLASS``.
The forward pass maps ``[B, 7180] float32`` (3-mer counts + MHC one-hot) to
``[B, 2]`` logits, matching ``forward_contract.model_io`` in
``declared/task_config.yaml``.
"""

from __future__ import annotations

import torch
from pydantic import BaseModel, Field
from torch import nn

#: 7,035 3-mer features + 145 MHC one-hot features (declared/kmer_vocab.json).
FEATURE_COUNT = 7180
CLASS_COUNT = 2

PLUGIN_MODEL_TYPE = "tcr_reference_mlp"
PLUGIN_OUTPUT_TYPE = "classifier"


class TcrReferenceMlpConfig(BaseModel):
    """Settings the planner may choose. ``extra`` stays permissive on purpose:
    planner-emitted fields this config does not model are ignored, not refused."""

    #: The planner sends ``model_type`` in every plan; a plugin config declares it.
    model_type: str = Field(default="tcr_reference_mlp")
    #: Read by the training engine at start-up; a tabular input has no segments.
    segmentation_size: int = Field(default=FEATURE_COUNT, ge=1)
    batch_size: int = Field(default=32, ge=1)
    hidden_dim: int = Field(default=32, ge=2, le=1024)


PLUGIN_CONFIG_CLASS = TcrReferenceMlpConfig


class TcrReferenceMlp(nn.Module):
    def __init__(self, config: TcrReferenceMlpConfig):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(FEATURE_COUNT, config.hidden_dim),
            nn.ReLU(),
            nn.Linear(config.hidden_dim, CLASS_COUNT),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Raw logits: the cross-entropy loss and the deliverable codec (softmax) expect them."""
        return self.net(x)


PLUGIN_MODEL_CLASS = TcrReferenceMlp
