#!/usr/bin/env python3
"""Compare primers designed in Primer Express with this tool's rules and Tm.

Primer Express does not publish its Tm formula or its penalty score, so the
only way to line the two up is with primers it actually designed. Give this
script a tab-separated file, one primer per line:

    name  sequence  [Primer Express Tm]  [Primer Express penalty]

(a header line is fine, and the last two columns are optional), and it prints,
for each primer, the Tm this tool computes, the difference from the Primer
Express Tm, and each Primer Express rule as this tool applies it. The summary
gives the mean and spread of the Tm difference: a steady offset means the
Tm window (--tm) should be shifted by that much to select what Primer Express
would select.

    python scripts/compare_primer_express.py pe_primers.tsv
    python scripts/compare_primer_express.py pe_primers.tsv --salt 50 3 0.8 200
"""
from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import primer3  # noqa: E402

import qpcrprimer_batch as qpb  # noqa: E402


def read_rows(path: str) -> list[dict]:
    rows = []
    for line in Path(path).read_text().splitlines():
        f = [x.strip() for x in line.split("\t")]
        if len(f) < 2 or not f[1] or not set(f[1].upper()) <= set("ACGT"):
            continue                                  # header, blank or not a sequence
        num = lambda i: float(f[i].replace(",", ".")) if len(f) > i and f[i] else None  # noqa: E731
        rows.append({"name": f[0], "seq": f[1].upper(), "pe_tm": num(2), "pe_penalty": num(3)})
    return rows


def check(seq: str, p: qpb.Params) -> dict:
    tm = round(primer3.calc_tm(seq, mv_conc=p.salt_mono, dv_conc=p.salt_div,
                               dntp_conc=p.dntp, dna_conc=p.dna_conc), 1)
    gc = 100 * sum(c in "GC" for c in seq) / len(seq)
    return {"tm": tm, "gc": round(gc, 1), "end_gc": qpb.end_gc(seq), "run": qpb.max_run(seq),
            "rules": {"Tm": p.tm_min <= tm <= p.tm_max,
                      "GC": p.gc_min <= gc <= p.gc_max,
                      "3'G/C": qpb.end_gc(seq) <= p.max_end_gc,
                      "run": qpb.max_run(seq) <= p.max_poly_x}}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv")
    d = qpb.Params()
    ap.add_argument("--salt", type=float, nargs=4, metavar=("MONO_MM", "DIV_MM", "DNTP_MM", "DNA_NM"),
                    default=[d.salt_mono, d.salt_div, d.dntp, d.dna_conc])
    a = ap.parse_args(argv)
    p = qpb.Params(salt_mono=a.salt[0], salt_div=a.salt[1], dntp=a.salt[2], dna_conc=a.salt[3])
    rows = read_rows(a.tsv)
    if not rows:
        ap.error("no primer found (expected: name <TAB> sequence [<TAB> PE Tm [<TAB> PE penalty]])")
    print("name\tsequence\tlen\tTm_here\tTm_PE\tdiff\tGC%\t3'G/C\tlongest_run\tfails")
    diffs = []
    for r in rows:
        c = check(r["seq"], p)
        diff = None if r["pe_tm"] is None else round(c["tm"] - r["pe_tm"], 1)
        if diff is not None:
            diffs.append(diff)
        fails = ",".join(k for k, ok in c["rules"].items() if not ok) or "-"
        print(f'{r["name"]}\t{r["seq"]}\t{len(r["seq"])}\t{c["tm"]}\t'
              f'{"" if r["pe_tm"] is None else r["pe_tm"]}\t{"" if diff is None else diff}\t'
              f'{c["gc"]}\t{c["end_gc"]}\t{c["run"]}\t{fails}')
    print(f"\n{len(rows)} primers; conditions {a.salt[0]} mM Na+, {a.salt[1]} mM Mg2+, "
          f"{a.salt[2]} mM dNTP, {a.salt[3]} nM oligo", file=sys.stderr)
    if len(diffs) >= 2:
        print(f"Tm here - Tm Primer Express: mean {statistics.mean(diffs):+.2f} C, "
              f"SD {statistics.stdev(diffs):.2f} C, range {min(diffs):+.1f} to {max(diffs):+.1f} "
              f"(n={len(diffs)})", file=sys.stderr)
        print("A steady offset (small SD) means: shift --tm by the mean to select what Primer Express "
              "would select.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
