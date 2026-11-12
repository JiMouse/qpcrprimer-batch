#!/usr/bin/env python3
"""Merge the per-chunk outputs of a sharded qpcrprimer_batch run.

The gene list is too large for one job, so CI splits it into chunks, designs
each chunk separately and uploads the results. This script stitches them back
into one site/ directory and one Excel workbook.

  python scripts/merge_parts.py parts/ site/

`parts/` holds one directory per chunk, each as produced by
`qpcrprimer_batch.py --json-dir <chunk>/data -o <chunk>/part.xlsx`
(the custom regions are one more such directory).
The workbook is written with openpyxl's write-only writer, so memory stays
flat no matter how many amplicons there are.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Alignment, Font, PatternFill

FONT = Font(name="Arial", size=10)
MONO = Font(name="Consolas", size=10)
BOLD = Font(name="Arial", size=10, bold=True, color="FFFFFF")
HEAD_FILL = PatternFill("solid", start_color="3B2F5C")
STATUS_FILL = {"OK": PatternFill("solid", start_color="D5EFD9"),
               "CHECK": PatternFill("solid", start_color="FDEBC8"),
               "FAIL": PatternFill("solid", start_color="F7D3D0")}
# Sequences read better in a monospaced font; column indices are per sheet.
SEQ_COLS = {"Primers": {7, 14}, "Alternatives": {4, 7}, "Amplicons": {5}}  # = qpcrprimer_batch.SEQ_COLS


def gene_bounds(part: Path, name: str) -> tuple[int, int] | None:
    """Transcript span of a gene, so the site can find a gene from a coordinate."""
    f = part / "data" / "genes" / f"{name}.json"
    if not f.exists():
        return None
    ex = json.loads(f.read_text()).get("exons") or []
    if not ex:
        return None
    return min(e["start"] for e in ex), max(e["end"] for e in ex)


def merge_index(parts: list[Path]) -> dict:
    """Concatenate the per-chunk index.json files, keeping run-wide settings."""
    merged: dict | None = None
    genes: list[dict] = []
    errors: list[dict] = []
    for p in parts:
        f = p / "data" / "index.json"
        if not f.exists():
            print(f"  warning: {f} missing, chunk skipped", file=sys.stderr)
            continue
        d = json.loads(f.read_text())
        # Each gene's span lets the site's region search pick genes without
        # fetching every gene file; new builds write it, older parts get it here.
        for g in d.get("genes", []):
            b = None if "start" in g else gene_bounds(p, g["gene"])
            if b:
                g["start"], g["end"] = b
        if merged is None:
            merged = {k: v for k, v in d.items() if k not in ("genes", "errors")}
        elif not d.get("ispcr", False):
            # One chunk without a specificity check makes the whole run unchecked.
            merged["ispcr"] = False
        genes += d.get("genes", [])
        errors += d.get("errors", [])
    if merged is None:
        sys.exit("no chunk produced an index.json")
    # genes first, custom regions after, each alphabetically
    merged["genes"] = sorted(genes, key=lambda g: (g.get("kind") == "region", g["gene"]))
    merged["errors"] = sorted(errors, key=lambda e: e["gene"])
    return merged


def merge_excel(parts: list[Path], out: Path):
    """Append every chunk's rows into one write-only workbook."""
    books = [(p, load_workbook(x, read_only=True))
             for p in parts for x in [p / "part.xlsx"] if x.exists()]
    if not books:
        sys.exit("no chunk produced a part.xlsx")
    wb = Workbook(write_only=True)
    sheets = books[0][1].sheetnames
    counts: dict[str, int] = {}
    for name in sheets:
        ws = wb.create_sheet(name)
        seq = SEQ_COLS.get(name, set())
        rows = 0
        for i, (_, src) in enumerate(books):
            if name not in src.sheetnames:
                continue
            for j, row in enumerate(src[name].iter_rows(values_only=True)):
                # The header, and the single-run Parameters sheet, come from
                # the first chunk only; everything else is per-gene data.
                if j == 0:
                    if i:
                        continue
                    ws.append([head(ws, v) for v in row])
                    continue
                if name == "Parameters" and i:
                    continue
                ws.append([body(ws, v, k + 1, row, seq) for k, v in enumerate(row)])
                rows += 1
        counts[name] = rows
        ws.freeze_panes = "A2"
    for _, src in books:
        src.close()
    wb.save(out)
    return counts


def head(ws, value):
    c = WriteOnlyCell(ws, value=value)
    c.font, c.fill = BOLD, HEAD_FILL
    c.alignment = Alignment(vertical="center", wrap_text=True)
    return c


def body(ws, value, col, row, seq):
    c = WriteOnlyCell(ws, value=value)
    c.font = MONO if col in seq else FONT
    # Colour the row by its status column, as the single-run writer does.
    if len(row) > 2 and row[2] in STATUS_FILL:
        c.fill = STATUS_FILL[row[2]]
    return c


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, site = Path(sys.argv[1]), Path(sys.argv[2])
    parts = sorted(d for d in src.iterdir() if d.is_dir())
    if not parts:
        sys.exit(f"no chunk directories under {src}")
    print(f"merging {len(parts)} chunks")

    genes_dir = site / "data" / "genes"
    genes_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for p in parts:
        for f in sorted((p / "data" / "genes").glob("*.json")):
            shutil.copy2(f, genes_dir / f.name)
            n += 1
    print(f"  {n} gene files -> {genes_dir}")

    index = merge_index(parts)
    (site / "data" / "index.json").write_text(json.dumps(index, indent=1))
    print(f"  index.json: {len(index['genes'])} genes, "
          f"{len(index['errors'])} errors, ispcr={index['ispcr']}")

    counts = merge_excel(parts, site / "qpcr_primers.xlsx")
    print("  qpcr_primers.xlsx: " + ", ".join(f"{k} {v}" for k, v in counts.items()))


if __name__ == "__main__":
    main()
