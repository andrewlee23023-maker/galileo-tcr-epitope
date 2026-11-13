"""TCR-epitope task workflow: edit task paths and review budgets before preparing a run.

Adapted from the canonical SIDERIUS CSV classification template
(yuema137/SIDERIUS, examples/task_templates/csv_classification/workflow.py).
"""

from workflows.definition import (
    ArgumentBinding,
    ExternalRequest,
    ExternalSource,
    ManagedCandidateBinding,
    NodeInstance,
    ProtocolRequest,
    ResultSource,
    ReusableTaskReferences,
    ReusableWorkflowDefinition,
    ReusableWorkflowDocument,
)

document = ReusableWorkflowDocument(
    definition=ReusableWorkflowDefinition(
        # Paths relative to this workflow; data holds train.csv and validation.csv.
        # test_inputs.csv / test_labels.csv are evaluator-only and live outside
        # the run data folder.
        task=ReusableTaskReferences(manifest="task.yaml", data_root="data"),
        externals=(
            ManagedCandidateBinding(id="initial", contract="proposal.input"),
            ManagedCandidateBinding(id="build_store", contract="storage"),
            ManagedCandidateBinding(id="check_store", contract="storage"),
            ManagedCandidateBinding(id="training", contract="native-training.settings"),
        ),
        nodes=(
            NodeInstance(
                id="hypothesis", capability="proposal", request=ExternalRequest(binding="initial")
            ),
            NodeInstance(
                id="build",
                capability="implementation",
                request=ProtocolRequest(
                    protocol="proposal-to-implementation",
                    arguments=(
                        ArgumentBinding(
                            name="storage", source=ExternalSource(binding="build_store")
                        ),
                        ArgumentBinding(name="output", source=ResultSource(node="hypothesis")),
                    ),
                ),
            ),
            NodeInstance(
                id="check",
                capability="validation",
                request=ProtocolRequest(
                    protocol="implementation-to-validation",
                    arguments=(
                        ArgumentBinding(name="output", source=ResultSource(node="build")),
                        ArgumentBinding(
                            name="storage", source=ExternalSource(binding="check_store")
                        ),
                    ),
                ),
            ),
            NodeInstance(
                id="fit",
                capability="tuning",
                request=ProtocolRequest(
                    protocol="validated-candidate-to-tuning",
                    arguments=(
                        ArgumentBinding(name="proposal", source=ResultSource(node="hypothesis")),
                        ArgumentBinding(name="validation", source=ResultSource(node="check")),
                        ArgumentBinding(name="settings", source=ExternalSource(binding="training")),
                    ),
                ),
            ),
        ),
    ),
    training_settings={
        "launch": {
            # Budgets for 61,585 training rows: longer than the tiny CSV
            # template's 2/5 minutes. Review before preparing a run.
            "trial_time_budget_minutes": 15.0,
            "formal_time_budget_minutes": 60.0,
            "healthgate_mode": "blocking",
            "result_authority": "scientific",
            "is_trial": True,
            "formal_portion": 1.0,
            "max_rounds": 3,
        },
    },
)
