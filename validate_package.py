"""Validate the SIDERIUS package WITHOUT a SIDERIUS checkout or GPU.

Stubs the SIDERIUS imports (core.*, execute_tools.*, torch) with minimal
fakes, then exercises the real package code on the real package data:
- YAML/JSON parse of every declared file, task.yaml, llm.json
- data path: read train/validation, scope build/serialize round-trip,
  training/validation dataset shapes, deliverable write/read round-trip
- metrics: macro AUPRC / pooled AUPRC / macro AUROC vs scikit-learn on
  validation predictions from a logreg trained on the package train file
- data hygiene: test_inputs has no label / is_synthetic_negative; row_id unique
"""
import json, sys, types, importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

PKG = Path.home() / "workspace/goals/galileo-challenge/files/track1-tcr-epitope/siderius_package"
ok = []
def check(name, cond, extra=""):
    ok.append((name, bool(cond)))
    print(("PASS " if cond else "FAIL ") + name + (f" — {extra}" if extra else ""))

# 1. parse all YAML/JSON
for f in [PKG/"task.yaml", PKG/"declared/task_config.yaml", PKG/"declared/task_health.yaml"]:
    yaml.safe_load(f.read_text()); check(f"yaml parses: {f.name}", True)
for f in [PKG/"llm.json", PKG/"declared/dataset_profile.json", PKG/"declared/metric_macro_auprc.json",
          PKG/"declared/metric_pooled_auprc.json", PKG/"declared/metric_macro_auroc.json",
          PKG/"declared/kmer_vocab.json"]:
    json.loads(f.read_text()); check(f"json parses: {f.name}", True)

# manifest sections required by the canonical composition reference
man = yaml.safe_load((PKG/"task.yaml").read_text())
for sec in ["task_data_path","dataset_profile","metric","task_config","task_health"]:
    check(f"task.yaml has required section {sec}", sec in man)
check("task.yaml metric id", man["metric"]["declaration"].endswith("metric_macro_auprc.json"))
check("task.yaml secondary ids unique", len({s["declaration"] for s in man["secondary_metrics"]})==2)
check("task.yaml model_plugins require", man["model_plugins"]["require"]==["tcr_reference_mlp"])
check("task.yaml objective ce", man["objective"]["config"]["loss_type"]=="ce")

# 2. stub SIDERIUS + torch imports
def stub_module(name, **attrs):
    m = types.ModuleType(name); m.__dict__.update(attrs); sys.modules[name]=m; return m
class MissingRunSeed(Exception): pass
stub_module("core"); stub_module("core.seeding", MissingRunSeed=MissingRunSeed)
# fake torch: from_numpy -> np array; TensorDataset stores pairs
class FakeTensor:
    def __init__(self, a): self.a=np.asarray(a)
    @property
    def ndim(self): return self.a.ndim
    def float(self): return self
    def item(self): return float(self.a)
    def __getitem__(self, i): return self.a[i]
    def __len__(self): return len(self.a)
class FakeTD:
    def __init__(self, x, y): self.x, self.y = x, y
    def __len__(self): return len(self.x.a)
torch = stub_module("torch")
torch.Tensor = FakeTensor
torch.from_numpy = lambda a: FakeTensor(a)
torch.utils = stub_module("torch.utils")
torch.utils.data = stub_module("torch.utils.data", TensorDataset=FakeTD, Dataset=object)
torch.nn = stub_module("torch.nn")
torch.nn.functional = stub_module("torch.nn.functional")
def _softmax(t, dim=-1):
    a=np.asarray(t.a if isinstance(t, FakeTensor) else t, dtype=float)
    e=np.exp(a-a.max(axis=dim, keepdims=True)); return FakeTensor(e/e.sum(axis=dim, keepdims=True))
torch.nn.functional.softmax=_softmax
torch.argmax=lambda t: FakeTensor(np.argmax(np.asarray(t.a), axis=-1)) if isinstance(t, FakeTensor) else np.argmax(t)
# execute_tools stubs: simple dataclasses for the request types
et = stub_module("execute_tools")
def _req(**kw):
    from dataclasses import dataclass
    @dataclass
    class R:
        pass
    for k,v in kw.items(): setattr(R, k, v)
    return R
class _Base:
    def __init__(self, **kw): self.__dict__.update(kw)
for nm in ["DeliverableWriteRequest","EpochSamplingParams","EvalMaterializationParams","EvaluationReadRequest","HealthCoverageRequest","HealthCoverageResult","ScopeBuildRequest","StorageReadScope","TaskOutputArtifactInventory","ValidationScopeError"]:
    setattr(et, nm, type(nm, (_Base,), {})) if nm!="ValidationScopeError" else None
et.ValidationScopeError = type("ValidationScopeError", (Exception,), {})
tdp = stub_module("execute_tools.task_data_path")
for nm in ["DeliverableWriteRequest","EpochSamplingParams","EvalMaterializationParams","EvaluationReadRequest","HealthCoverageRequest","HealthCoverageResult","ScopeBuildRequest","StorageReadScope","TaskOutputArtifactInventory","ValidationScopeError"]:
    setattr(tdp, nm, getattr(et, nm))
tdp.require_bound_task_data_path = lambda: _BOUND[0]
_BOUND=[None]
dspec = stub_module("execute_tools.deliverable_spec")
class _Naming:
    def name(self, model_type, run_name, exp_id, input_identity=0):
        return f"tcr_predictions_{exp_id}_{input_identity:04d}.json"
dspec.resolve_deliverable_naming = lambda: _Naming()
em = stub_module("execute_tools.evaluation_metric")
class EvaluationMetric:
    def __init__(self): pass
em.EvaluationMetric = EvaluationMetric
hc = stub_module("execute_tools.health_checks")
hc.CATEGORICAL_PREDICTIONS="categorical_predictions"
class _Payload:
    def __init__(self, **kw): self.__dict__.update(kw)
hc.CategoricalPredictionsPayload=_Payload
hc.HealthView=_Payload
hc.register_view_provider=lambda cls: None
hcs = stub_module("execute_tools.health_checks.schemas")
hcs.HealthCheckContext=_Base

def load(name, path):
    spec=importlib.util.spec_from_file_location(name, path)
    mod=importlib.util.module_from_spec(spec); sys.modules[name]=mod; spec.loader.exec_module(mod); return mod

# 3. data path on real package data
task_mod = load("_tcr_task", PKG/"plugins/_tcr_task.py")
vocab = json.loads((PKG/"declared/kmer_vocab.json").read_text())
dp = task_mod.TcrEpitopeTaskDataPath(train_file="train.csv", validation_file="validation.csv",
    label_column="label", classes=[0,1], id_column="row_id", vocab_file=str(PKG/"declared/kmer_vocab.json"))
check("feature_count 7180", dp.feature_count==7180, str(dp.feature_count))
data_dir=str(PKG/"data")
tr = dp.read_table(data_dir, "train.csv"); va = dp.read_table(data_dir, "validation.csv")
check("train rows 61585", len(tr.sample_ids)==61585, str(len(tr.sample_ids)))
check("val rows 3042", len(va.sample_ids)==3042, str(len(va.sample_ids)))
check("features shape", tr.features.shape==(61585,7180), str(tr.features.shape))
check("val epitopes", set(va.epitopes)=={"MPASWVMRI","VLWAHGFEL"}, str(set(va.epitopes)))
check("train/val epitopes disjoint", not (set(tr.epitopes)&set(va.epitopes)))
# free the big train table before the dataset test re-parses it (sandbox RAM)
import gc as _gc
del tr; _gc.collect()
# scope round-trip
from types import SimpleNamespace
scope = dp.build_training_scope(SimpleNamespace(selection_strategy="default", target_partitions=(0,), subset_ref=None, seed=None, portion=1.0, max_samples=None))
s2 = dp.deserialize_scope(dp.serialize_scope(scope))
check("scope serialize round-trip", s2==scope)
# dataset
ds = dp.training_dataset(scope, SimpleNamespace(data_dir=data_dir, train_portion=1.0, max_samples=None, epoch_seed=7))
check("training dataset len", len(ds)==61585)
# deliverable round-trip with fake logits
import tempfile, os
with tempfile.TemporaryDirectory() as td:
    req = SimpleNamespace(task_scope=dp.build_eval_scope(SimpleNamespace(selection_strategy="default", target_partitions=(1,), subset_ref=None, seed=None, portion=1.0, max_samples=None)),
        source_context=SimpleNamespace(data_dir=data_dir, sample_count=3042),
        output_dir=td, model_type="tcr_reference_mlp", run_name="r1", exp_id="e1")
    # fake outputs: logits as FakeTensor per row (class 1 logit = label)
    outs=[FakeTensor(np.array([0.0, float(l)])) for l in va.labels[:3042]]
    dp.write_deliverable(outs, req)
    files=os.listdir(td); check("deliverable written", len(files)==1, str(files))
    rreq=SimpleNamespace(deliverable_dir=td, model_type="tcr_reference_mlp", run_name="r1", exp_id="e1")
    payload=dp.read_evaluation_payload(rreq)
    check("deliverable round-trip", len(payload)==3042 and all(0<=v["p_binding"]<=1 for v in payload.values()))
    _BOUND[0]=dp
    # metrics vs sklearn
    from sklearn.metrics import average_precision_score, roc_auc_score
    met = load("_tcr_metrics", PKG/"plugins/_tcr_metrics.py")
    scope_v = dp.build_eval_scope(SimpleNamespace(selection_strategy="default", target_partitions=(1,), subset_ref=None, seed=None, portion=1.0, max_samples=None))
    for cls_name, sk_fn in [("MacroAuprcMetric", average_precision_score), ("MacroAurocMetric", roc_auc_score)]:
        m=getattr(met, cls_name)()
        val,_,_=m._compute({}, evaluation_payload=payload, task_scope=scope_v, data_dir=data_dir)
        # sklearn macro: mean per epitope
        labs=np.array([va.labels[i] for i in range(3042)]); ps=np.array([payload[va.sample_ids[i]]["p_binding"] for i in range(3042)])
        sk=np.mean([sk_fn(labs[np.array(va.epitopes)==e], ps[np.array(va.epitopes)==e]) for e in set(va.epitopes)])
        check(f"{cls_name} matches sklearn", abs(val-sk)<1e-6, f"ours={val:.4f} sklearn={sk:.4f}")
    m=met.PooledAuprcMetric(); val,_,_=m._compute({}, evaluation_payload=payload, task_scope=scope_v, data_dir=data_dir)
    labs=np.array(va.labels); ps=np.array([payload[va.sample_ids[i]]["p_binding"] for i in range(3042)])
    check("PooledAuprc matches sklearn", abs(val-average_precision_score(labs, ps))<1e-6, f"{val:.4f}")

# 4. data hygiene
ti=pd.read_csv(PKG/"data/test_inputs.csv", dtype=str, keep_default_na=False)
check("test_inputs has no label", "label" not in ti.columns)
check("test_inputs has no is_synthetic_negative", "is_synthetic_negative" not in ti.columns)
check("test_inputs row_id unique", ti["row_id"].is_unique)
check("test_inputs rows 2876", len(ti)==2876)
check("test_labels in evaluation/, not data/", (PKG/"evaluation/test_labels.csv").exists() and not (PKG/"data/test_labels.csv").exists())
tl=pd.read_csv(PKG/"evaluation/test_labels.csv")
check("test labels match test inputs row_id", set(tl["row_id"])==set(ti["row_id"]))

# 5. syntax compile all plugins + workflow
import py_compile
for f in list((PKG/"plugins").glob("*.py"))+[PKG/"workflow.py"]:
    py_compile.compile(str(f), doraise=True); check(f"compiles: {f.name}", True)

print()
fails=[n for n,c in ok if not c]
print(f"{len(ok)-len(fails)}/{len(ok)} checks passed")
if fails: print("FAILURES:", fails); sys.exit(1)
