# Primer Express

What the manual workflow used, and how this tool reproduces it.
[← back to the README](../README.md)

## What Primer Express does

Primer Express (Applied Biosystems) designs real-time PCR primers, with or without a TaqMan probe,
from a sequence you paste in. It lists every primer set that meets the parameters of the chosen
document and sorts them by a **penalty score**, lowest first. The manual workflow this repository
automates: paste the target, keep the sets with a penalty under 100, then check SNPs and
uniqueness elsewhere.

## Measured, not guessed

Applied Biosystems does not publish the penalty formula or the Tm formula. Both were measured on a
Primer Express 2.0 export of 200 primer sets (a TaqMan MGB document, 47 distinct primers), and
`scripts/compare_primer_express.py` repeats the measurement on any export:

| | formula | fit |
|---|---|---|
| **Penalty** | **5 × (amplicon − 50) + \|F length − 20\| + \|R length − 20\|** | 200 / 200 sets identical |
| **Tm** | Breslauer 1986 nearest-neighbour, Schildkraut–Lifson salt term, 50 mM Na⁺, 5 nM oligo, no Mg²⁺ | 44 / 47 primers identical after rounding; the other 3 within 0.5 °C (rounding ties) |

Within the rules, Tm, GC and the probe add nothing to the penalty: **it is almost entirely
amplicon length**. Five points per base pair above 50 means the usual "penalty < 100" cut-off
amounts to an amplicon of about 69 bp or less with 20-mer primers.

Both formulas are what this tool uses: Primer3 selects primers on the Primer Express Tm, and its
penalty weights are set so that Primer3's penalty *is* the Primer Express penalty (only amplicon
size and primer length weigh; Tm, GC and ΔTm are rules, not scores). The penalty on the site and in
Excel is the number Primer Express would show for the same pair.

## The rules

| Primer Express | Here |
|---|---|
| Tm 58–60 °C | Primer3 `PRIMER_MIN/MAX_TM` = 58 / 60, on the Primer Express Tm |
| max Tm difference 2 °C | `PRIMER_PAIR_MAX_DIFF_TM` = 2 |
| GC 30–80 % | `PRIMER_MIN/MAX_GC` = 30 / 80 |
| no more than two G/C among the last five 3′ bases | `PRIMER_MAX_END_GC` = 2 |
| fewer than four consecutive G | no run of 4 G (checked on each candidate) |
| runs of other bases | up to 5 (`PRIMER_MAX_POLY_X` = 5): the export contains runs of five C and five A |
| optimal length 20 | penalty 1 per nt away from 20; lengths 16–30 (the export has 16 to 27) |
| amplicon 50–150 bp | product 50–150; penalty 5 per bp above 50 |

Every one of the export's 47 primers passes all of these rules as applied here.

## Differences that remain, on purpose

- **No probe.** This is the SYBR Green workflow. A TaqMan document must also fit a probe between
  the primers, which forces longer amplicons: on TBX1 exon 2 the best set of that export scored 30
  (55 bp), while the primer-only design here finds a pair scoring 4 (50 bp). Primer Express's own
  primer-only document should find such pairs too.
- **Hairpins and primer-dimers** are judged as Primer3 judges them — the structure's Tm must stay
  under 47 °C — in a PCR buffer (50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer). Primer
  Express's own secondary-structure test is not published; in SYBR Green any dimer fluoresces, so
  this check stays strict.
- **Ranking around the penalty.** SNP- and repeat-free primers come first and pairs specific in
  in-silico PCR are preferred; the penalty orders the rest. The search covers the exon and 100 bp
  of intron on each side, as if that whole stretch had been pasted into Primer Express, and the
  lowest penalty wins wherever it sits.
- **The cut-off is a flag, not a filter.** When no pair under 100 passes the checks, the best one
  is still given, marked Check with its penalty.

## Checking another export

Save the "Primers" table from Primer Express as text and run:

```bash
python scripts/compare_primer_express.py "MyGene_Primers.TXT"
```

It prints, for every primer, this tool's Tm next to Primer Express's and any rule the primer would
fail here, then how many Tm and penalties are identical. A plain list (`name <TAB> sequence <TAB>
Tm`) works too. If another version or document of Primer Express ever disagrees, the summary shows
it at once.
