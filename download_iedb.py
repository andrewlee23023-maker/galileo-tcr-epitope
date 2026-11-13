#!/usr/bin/env python3
"""Attempt the documented IEDB full-export download, and validate the result.

Documented source (verified Nov 8, 2026):
  https://www.iedb.org/database_export_v3.php
  Weekly builds named there: tcr_full_v3.xlsx, bcr_full_v3.xlsx, airr_full.zip

Direct download (verified Nov 8, 2026 — the missing piece is the `doc/`
prefix; without it the URL returns HTML):
  https://www.iedb.org/downloader.php?file_name=doc/tcr_full_v3.xlsx
  https://www.iedb.org/downloader.php?file_name=doc/airr_full.zip
Both were fetched successfully by curl in this sandbox (tcr_full_v3.xlsx
31.1 MB via the browser, airr_full.zip 9.9 MB via curl).

IMPORTANT: tcr_full_v3.xlsx is positive-only (302,807 rows, no outcome
column). For binding/nonbinding outcomes use airr_full.zip instead —
its tcr_full_airr.tsv has iedb_assay_outcome, but real negatives are
only 1,410 of 339,842 rows, so synthetic negatives are required
(see normalise_airr.py).

Usage:
  python3 download_iedb.py --out data/
  # then, after a browser download:
  python3 normalise_iedb.py --inspect data/tcr_full_v3.xlsx
"""
import argparse
import re
import sys
import urllib.request

EXPORT_PAGE = "https://www.iedb.org/database_export_v3.php"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data")
    args = ap.parse_args()
    req = urllib.request.Request(EXPORT_PAGE, headers={"User-Agent": "Mozilla/5.0"})
    try:
        html = urllib.request.urlopen(req, timeout=30).read().decode("utf-8", "ignore")
    except Exception as e:
        print(f"ERROR: could not fetch {EXPORT_PAGE}: {e}", file=sys.stderr)
        print("Download in a browser instead: open the export page, download "
              "'tcr_full_v3.xlsx' (or 'airr_full.zip'), save it under data/.", file=sys.stderr)
        sys.exit(1)
    links = sorted(set(re.findall(r'href="([^"]*(?:tcr_full|airr_full)[^"]*)"', html)))
    if not links:
        print("No direct download link for tcr_full_v3.xlsx / airr_full.zip was "
              "found on the export page (it is JS-driven).")
        print("Download in a browser:")
        print(f"  1. Open {EXPORT_PAGE}")
        print("  2. Download 'tcr_full_v3.xlsx' (full TCR export) — or 'airr_full.zip'")
        print("  3. Save it as data/tcr_full_v3.xlsx in this folder")
        print("  4. Run: python3 normalise_iedb.py --inspect data/tcr_full_v3.xlsx")
        sys.exit(3)
    print("Found candidate links:")
    for link in links:
        print(f"  - {link}")
    print("Download one with a browser if the link is not a direct file URL, "
          "then run normalise_iedb.py --inspect on it.")


if __name__ == "__main__":
    main()
