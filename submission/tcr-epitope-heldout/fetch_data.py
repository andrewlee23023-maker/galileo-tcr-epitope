"""Download, verify and prepare the TCR-epitope data.

Usage (any Python 3.12+; standard library only)::

    python fetch_data.py --output_dir DIR                 # download from Zenodo
    python fetch_data.py --output_dir DIR --source_dir S  # use files already in S

Reads ``data_source.json`` next to this script. Every file listed there is
checked for size and md5 before use. The Zenodo record holds three files:

    train.csv        61,585 rows, 2,500 epitopes, 49.6% positive
    validation.csv    3,042 rows, 2 epitopes (MPASWVMRI, VLWAHGFEL), 50.0% positive
    test_inputs.csv   2,876 rows, 2 epitopes (SPRWYFYYL, TFEYVSQPFLMDLE), no labels

This script writes the layout that ``task/task.yaml`` reads:

    DIR/train.csv
    DIR/validation.csv

``test_inputs.csv`` is verified (checksum) but NOT written to the data
folder: it is the held-out test input, and nothing in the research loop may
read it. Test labels are not on Zenodo at all; they are evaluator-only
(``evaluation/test_labels.csv`` in the submission) and are never published.

Re-running is safe: an identical output file is kept; a different one is
refused and left unchanged. Exit status 0 on success, 2 on any refusal.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
API = "https://zenodo.org/api/records/{record}/files/{name}/content"

#: The two files the run data folder needs. test_inputs.csv is verified
#: below but deliberately not written (see module docstring).
OUTPUT_FILES = ("train.csv", "validation.csv")


class Refused(RuntimeError):
    pass


def md5_of(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path: Path, size: int, md5: str) -> None:
    if not path.is_file():
        raise Refused(f"{path} is missing")
    if path.stat().st_size != size:
        raise Refused(f"{path}: {path.stat().st_size} bytes, data_source.json says {size}")
    if md5_of(path) != md5:
        raise Refused(f"{path}: md5 differs from data_source.json")


def download(source: dict, name: str, target: Path) -> None:
    record = source["zenodo_record_id"]
    if not record.isdigit():
        raise Refused(f"zenodo_record_id {record!r} is a placeholder; use --source_dir")
    url = API.format(record=record, name=urllib.parse.quote(name))
    if source.get("share_link"):
        token = urllib.parse.parse_qs(urllib.parse.urlparse(source["share_link"]).query).get("token")
        if token:
            url += "?" + urllib.parse.urlencode({"token": token[0]})
    print(f"downloading {name}", flush=True)
    request = urllib.request.Request(url, headers={"User-Agent": "siderius-hackathon-fetch/1"})
    with urllib.request.urlopen(request, timeout=120) as response, target.open("wb") as sink:
        shutil.copyfileobj(response, sink)


def write_once(path: Path, payload: bytes) -> None:
    if path.exists():
        if path.read_bytes() == payload:
            print(f"kept {path} (already identical)")
            return
        raise Refused(f"{path} exists with different content; left unchanged")
    temporary = path.with_name(path.name + ".partial")
    temporary.write_bytes(payload)
    os.replace(temporary, path)
    print(f"wrote {path}")


def main(argv: list[str] | None = None) -> int:
    cli = argparse.ArgumentParser(description="Download, verify and prepare the task data.")
    cli.add_argument("--output_dir", type=Path, required=True, help="the data folder to write")
    cli.add_argument("--source_dir", type=Path, help="folder holding the Zenodo files already")
    args = cli.parse_args(argv)
    source = json.loads((HERE / "data_source.json").read_text(encoding="utf-8"))
    files = {entry["name"]: entry for entry in source["files"]}
    try:
        with tempfile.TemporaryDirectory() as scratch:
            raw = Path(args.source_dir) if args.source_dir else Path(scratch)
            for name, entry in files.items():
                if args.source_dir is None:
                    download(source, name, raw / name)
                verify(raw / name, entry["size"], entry["md5"])
            print(f"verified {len(files)} file(s) against data_source.json")
            args.output_dir.mkdir(parents=True, exist_ok=True)
            for name in OUTPUT_FILES:
                write_once(args.output_dir / name, (raw / name).read_bytes())
            print("test_inputs.csv verified but not written to the data folder (held-out test)")
    except (Refused, OSError) as error:
        print(f"fetch_data refused: {error}", file=sys.stderr)
        return 2
    print(f"Data ready in {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
