#!/usr/bin/env python3
"""Check this tool against primers designed in Primer Express.

Give it either

  * a Primer Express export - the "Primers" table saved as text, one primer
    set per line (Start, Length, Tm, %GC, Primer, [probe columns], ...,
    amplicon Length, Tm, %GC, Ta, Penalty); the probe is ignored; or
  * a plain tab-separated list:  name  sequence  [Primer Express Tm]

and it prints, for every primer, the Tm this tool computes next to Primer
Express's, and every rule this tool would apply to it. For an export it also
recomputes each set's penalty. The summary says how many Tm and penalties
match. Spaces inside a sequence (a pasting artefact) are ignored.

    python scripts/compare_primer_express.py "My gene_Primers.TXT"
    python scripts/compare_primer_express.py pe_primers.tsv --tm-formula santalucia
"""
from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import primer3  # noqa: E402

import qpcrprimer_batch as qpb  # noqa: E402


def _num(x: str) -> float | None:
    x = x.strip().replace(",", ".")
    try:
        return float(x)
    except ValueError:
        return None


def _seq(x: str) -> str | None:
    s = x.replace(" ", "").upper()
    return s if s and set(s) <= set("ACGT") else None


def read_file(path: str) -> tuple[dict, list]:
    """-> ({primer: Primer Express Tm or None}, [primer sets of an export])."""
    primers, sets = {}, []
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        f = line.rstrip("\r").split("\t")
        # Primer Express export: 20 columns (with probe) or 15 (primers only)
        if len(f) >= 15 and _seq(f[4]) and _num(f[2]) is not None:
            probe = len(f) >= 20
            r0 = 10 if probe else 5
            a0 = 15 if probe else 10
            F, R = _seq(f[4]), _seq(f[r0 + 4])
            if R is None:
                continue
            primers[F], primers[R] = _num(f[2]), _num(f[r0 + 2])
            sets.append({"F": F, "R": R, "product": int(_num(f[a0])), "penalty": _num(f[a0 + 4])})
        elif len(f) >= 2 and _seq(f[1]):
            primers[_seq(f[1])] = _num(f[2]) if len(f) > 2 else None
    return primers, sets


def check(seq: str, p: qpb.Params) -> dict:
    tm = round(primer3.calc_tm(seq, **p.tm_args()), 1)
    gc = 100 * sum(c in "GC" for c in seq) / len(seq)
    return {"tm": tm, "gc": round(gc, 1), "end_gc": qpb.end_gc(seq), "run": qpb.max_run(seq),
            "g_run": qpb.max_run(seq, "G"),
            "rules": {"Tm": p.tm_min <= tm <= p.tm_max,
                      "GC": p.gc_min <= gc <= p.gc_max,
                      "3'G/C": qpb.end_gc(seq) <= p.max_end_gc,
                      "run": qpb.max_run(seq) <= p.max_poly_x and qpb.max_run(seq, "G") <= p.max_poly_g,
                      "length": p.size_min <= len(seq) <= p.size_max}}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--tm-formula", choices=["primer-express", "santalucia"], default="primer-express")
    ap.add_argument("--salt", type=float, nargs=4, metavar=("MONO_MM", "DIV_MM", "DNTP_MM", "DNA_NM"))
    a = ap.parse_args(argv)
    p = qpb.Params(tm_formula=a.tm_formula, **qpb.salt_params(a.salt, a.tm_formula))
    primers, sets = read_file(a.file)
    if not primers:
        ap.error("no primer found: expected a Primer Express export, or name <TAB> sequence [<TAB> Tm]")
    print("sequence\tlen\tTm_here\tTm_PE\tdiff\tGC%\t3'G/C\tlongest_run\tlongest_G\tfails")
    diffs = []
    for s, pe in primers.items():
        c = check(s, p)
        diff = None if pe is None else round(c["tm"] - pe, 1)
        if diff is not None:
            diffs.append(diff)
        fails = ",".join(k for k, ok in c["rules"].items() if not ok) or "-"
        print(f'{s}\t{len(s)}\t{c["tm"]}\t{"" if pe is None else f"{pe:g}"}\t{"" if diff is None else diff}\t'
              f'{c["gc"]}\t{c["end_gc"]}\t{c["run"]}\t{c["g_run"]}\t{fails}')
    out = sys.stderr
    print(f"\n{len(primers)} distinct primers; Tm: {p.tm_formula}, {p.salt_mono} mM Na+, {p.salt_div} mM Mg2+, "
          f"{p.dntp} mM dNTP, {p.dna_conc} nM oligo", file=out)
    if diffs:
        same = sum(abs(round(t + d) - t) < 1e-9 for t, d in zip([pe for pe in primers.values() if pe is not None], diffs))
        line = f"Tm here - Tm Primer Express: mean {statistics.mean(diffs):+.2f} C"
        if len(diffs) >= 2:
            line += f", SD {statistics.stdev(diffs):.2f} C, range {min(diffs):+.1f} to {max(diffs):+.1f}"
        print(line + f"; identical once rounded: {same}/{len(diffs)}", file=out)
    if sets:
        ok = sum(p.penalty(x["product"], len(x["F"]), len(x["R"])) == x["penalty"] for x in sets)
        print(f"Penalty recomputed as 5 x (amplicon - 50) + |F len - 20| + |R len - 20|: "
              f"{ok}/{len(sets)} sets identical", file=out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
