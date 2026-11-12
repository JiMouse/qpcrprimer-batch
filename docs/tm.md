# Tm

Three temperatures are reported; each has a single use.
[← back to the README](../README.md)

| column | what it is | use it for |
|---|---|---|
| **Tm** (`F Tm`, `R Tm`) | SantaLucia 1998 nearest-neighbour Tm at 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo — the value the design was selected on | comparing primers, and the 58–60 °C rule |
| Tm UCSC | the number UCSC In-Silico PCR prints for the same oligo (Breslauer 1986 tables, 50 mM Na⁺, no Mg²⁺) | recognising the pair when you paste it into UCSC; nothing else |
| **Product Tm** | Primer3's estimate of the amplicon's melting temperature | what to expect of the melt-curve peak |

The conditions are those of exonprimer-batch, so a primer reads the same Tm in both tools; its
[Tm page](https://github.com/JiMouse/exonprimer-batch/blob/main/docs/tm.md) explains why UCSC
prints 5–8 °C more for a 20-mer, and why a calculator that models Na⁺ alone prints about 6 °C less.
Change them with `--salt MONO DIV DNTP OLIGO`.

## At the bench

The usual SYBR Green programme — 95 °C, then 40 cycles of 95 °C / 60 °C — suits primers in the
58–60 °C window. Primer Express's own Tm formula is not published (see
[Primer Express](primer-express.md)), so the same primer can read a degree or two differently
there; `scripts/compare_primer_express.py` measures the offset on primers it designed.

## Product Tm and the melt curve

Product Tm is an estimate from the amplicon's length and GC content in Primer3's model, not a
measurement: the peak on your instrument depends on the master mix and the dye. Use it to recognise
the right product — one peak, in the expected region — and to tell it from a primer-dimer peak,
which sits well below. A SNP between the primers (listed in Notes) can give a heterozygous sample
a slightly broader or shouldered peak.
