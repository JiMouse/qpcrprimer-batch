#!/usr/bin/env python3
"""Write genes.txt with every protein-coding gene that has a MANE Select transcript.

The list comes from NCBI's MANE summary file (one MANE Select transcript per
protein-coding gene), so every symbol in it is one the design resolves through
the same MANE Select choice. Lines of the current genes.txt are kept as they
are when their symbol is in MANE (transcript overrides survive); symbols that
are not in MANE are reported, and kept only with --keep-unknown.

Usage:
  python scripts/protein_coding_genes.py [--genes genes.txt] [--out genes.txt]
      [--summary MANE.GRCh38.vX.summary.txt.gz] [--keep-unknown]
"""
from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from pathlib import Path

import requests

MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"


def fetch_summary() -> tuple[str, bytes]:
    listing = requests.get(MANE_DIR, timeout=60)
    listing.raise_for_status()
    names = sorted(set(re.findall(r'MANE\.GRCh38\.v[\d.]+\.summary\.txt\.gz', listing.text)))
    if not names:
        sys.exit(f"no MANE summary file found at {MANE_DIR}")
    name = names[-1]
    r = requests.get(MANE_DIR + name, timeout=300)
    r.raise_for_status()
    return name, r.content


def mane_select_symbols(raw: bytes) -> list[str]:
    text = gzip.decompress(raw).decode()
    lines = text.splitlines()
    head = lines[0].lstrip("#").split("\t")
    col = {k: i for i, k in enumerate(head)}
    sym, status = col["symbol"], col["MANE_status"]
    out = set()
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) > status and f[status] == "MANE Select" and f[sym]:
            out.add(f[sym])
    return sorted(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", default="genes.txt", help="current gene list (its header and overrides are kept)")
    ap.add_argument("--out", default="genes.txt")
    ap.add_argument("--summary", help="a local MANE summary .txt.gz instead of downloading it")
    ap.add_argument("--keep-unknown", action="store_true", help="keep current symbols that are not in MANE")
    a = ap.parse_args(argv)

    if a.summary:
        name, raw = Path(a.summary).name, Path(a.summary).read_bytes()
    else:
        name, raw = fetch_summary()
    mane = mane_select_symbols(raw)
    mane_set = set(mane)

    header, current = [], {}
    if Path(a.genes).exists():
        for line in Path(a.genes).read_text().splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                if not current:
                    header.append(line)
                continue
            current[s.split()[0]] = line.rstrip()

    unknown = sorted(set(current) - mane_set)
    symbols = sorted(mane_set | (set(unknown) if a.keep_unknown else set()))
    body = [current.get(s, s) for s in symbols]
    Path(a.out).write_text("\n".join(header + body) + "\n")

    added = len(mane_set - set(current))
    print(f"{name}: {len(mane)} MANE Select genes; {len(current)} in the current list, "
          f"{added} added, {len(unknown)} not in MANE"
          + (" (kept)" if a.keep_unknown else " (dropped)") + f" -> {len(symbols)} genes in {a.out}")
    if unknown:
        print("not in MANE: " + " ".join(unknown))


if __name__ == "__main__":
    main()
