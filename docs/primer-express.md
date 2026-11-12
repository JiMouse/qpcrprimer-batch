# Primer Express

What the manual workflow used, what this tool reproduces, and what it cannot.
[← back to the README](../README.md)

## What Primer Express does

Primer Express (Applied Biosystems) designs real-time PCR primers, with or without a TaqMan probe,
from a sequence you paste in. It lists every primer set that meets the parameters of the chosen
document and sorts them by a **penalty score**, lowest first. The usual practice is the one this
repository automates: paste the target, take the lowest-penalty sets, then check SNPs and
uniqueness elsewhere.

## The rules, as applied here

The Primer Express 3.0 *Getting Started Guide* gives these guidelines for quantification primers;
each is a hard constraint here:

| Primer Express guideline | Here |
|---|---|
| Tm 58–60 °C, optimum 59 °C | Primer3 `PRIMER_MIN/OPT/MAX_TM` = 58 / 59 / 60 |
| GC 30–80 % | `PRIMER_MIN/MAX_GC` = 30 / 80 |
| optimal length 20 bases | `PRIMER_OPT_SIZE` = 20, range 18–28 |
| no more than two G/C among the last five 3′ bases | `PRIMER_MAX_END_GC` = 2 |
| avoid runs of identical bases; fewer than four consecutive G | `PRIMER_MAX_POLY_X` = 3, for every base, not only G |
| amplicon 50–150 bp, low amplicon length preferred | product 50–150, optimum 80 with a small penalty per bp away |

Primers only: this is the SYBR Green workflow. The TaqMan probe rules (Tm 68–70 °C, no G at the
5′ end, …) are not applied, because no probe is designed.

## The penalty score cannot be reproduced

Primer Express's penalty is **proprietary**. Thermo Fisher's own description is that it reflects
how many of the preferred parameters on the Parameters page each design meets — and, in the same
FAQ, that sets chosen by the TaqMan-specific document already meet the criteria, so *the score does
not and will not determine the primer/probe performance*. The formula is not published, so the
"penalty < 100" cut-off has no equivalent here.

What this tool does instead:

- every reported pair **already meets every rule above** — the role a low Primer Express penalty
  plays in practice;
- within that, pairs are ranked by **Primer3's penalty**: distance from the optimum Tm and length,
  the pair's ΔTm and the distance from the preferred amplicon size. It is on a different scale
  — on 435 exons of the panel the best pair scored a median of 1.6 and 95 % scored under 7 — so
  do not compare the two numbers;
- the ranking that matters most for genomic qPCR comes *after* the penalty: SNP-free and
  repeat-free primers first, then specific pairs first.

## Tm: calibrating against Primer Express

Primer Express does not publish its Tm formula either, and the 58–60 °C window was set against
*its* Tm. This tool computes the SantaLucia 1998 nearest-neighbour Tm under the conditions in
[Tm](tm.md) — the same as exonprimer-batch. If the two formulas differ by a steady offset, the
window here selects slightly different primers than Primer Express would.

That is measurable with primers Primer Express actually designed. Put them in a tab-separated
file — name, sequence, the Tm Primer Express showed, and optionally its penalty — and run:

```bash
python scripts/compare_primer_express.py pe_primers.tsv
```

It prints this tool's Tm next to Primer Express's for every primer, which rules each primer would
fail here, and the mean and spread of the Tm difference. A steady offset (small SD) means the
window should move by the mean, e.g. `--tm 57.5 58.5 59.5`, or the same through the workflow's
`extra_args`. Fifteen to twenty primers are enough to see it.
