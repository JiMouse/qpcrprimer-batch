# Tm

Three temperatures are reported; each has a single use.
[← back to the README](../README.md)

| column | what it is | use it for |
|---|---|---|
| **Tm** (`F Tm`, `R Tm`) | **Primer Express's Tm**: Breslauer 1986 nearest-neighbour with the Schildkraut–Lifson salt term, 50 mM Na⁺, 5 nM oligo, no Mg²⁺ — the value the design was selected on | the 58–60 °C rule, and comparing with Primer Express |
| Tm UCSC | what UCSC In-Silico PCR prints for the same oligo: the same formula at 50 nM | recognising the pair when you paste it into UCSC; it reads 3–5 °C higher |
| **Product Tm** | Primer3's estimate of the amplicon's melting temperature | what to expect of the melt-curve peak |

## Where the Primer Express formula comes from

It is not published; it was measured on a Primer Express export (see
[Primer Express](primer-express.md)). The Breslauer tables with the Schildkraut salt term are the
older formula UCSC still uses, which is why the two columns move together. At 5 nM oligo it gives
Primer Express's number for 44 of 47 primers once rounded, and is within 0.5 °C for the rest. Any
other export can be checked with `scripts/compare_primer_express.py`.

These are not the temperatures at which the primers anneal in your tube — every formula assumes
its own conditions, and a qPCR master mix has more salt and magnesium. They are a scale to select
primers on, and keeping Primer Express's scale means the 58–60 °C window selects what it selected.

To design on another scale instead — SantaLucia 1998 in a PCR buffer, as exonprimer-batch does —
use `--tm-formula santalucia` (default buffer 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo,
or your own with `--salt`). On that scale the same primers read about 1 °C higher than in Primer
Express.

## At the bench

The usual SYBR Green programme — 95 °C, then 40 cycles of 95 °C / 60 °C — suits primers in the
58–60 °C Primer Express window; that is the pairing Applied Biosystems designed it for.

## Hairpins and dimers

Secondary structures are not judged on the Primer Express scale but in a PCR buffer (50 mM Na⁺,
1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer), where they are more stable: a hairpin, self-dimer or
primer-dimer — anywhere, or anchored on a 3′ end — whose Tm exceeds 47 °C rejects the pair, as in
Primer3.

## Product Tm and the melt curve

Product Tm is an estimate from the amplicon's length and GC content in Primer3's model, not a
measurement: the peak on your instrument depends on the master mix and the dye. Use it to recognise
the right product — one peak, in the expected region — and to tell it from a primer-dimer peak,
which sits well below. A SNP between the primers (listed in Notes) can give a heterozygous sample
a slightly broader or shouldered peak.
