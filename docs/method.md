# Design method

How each assay is chosen, and what is checked before it is called OK.
[← back to the README](../README.md)

## Targets

1. **Exon structure** — the **MANE Select** transcript (RefSeq NM reported), from the UCSC `mane`
   track, one request per chromosome; Ensembl REST only for a transcript you name in `genes.txt`,
   a gene with no MANE transcript (Ensembl canonical is then used), or a row the track cannot
   provide. This is the exonprimer-batch code, unchanged, so both tools number exons identically.
2. **One target per coding exon** — its coding part by default (`--target coding`), or the whole
   exon (`--target full`). Non-coding exons are not targeted. Exons are numbered by their rank in
   the transcript.
3. **Custom regions** — each line of `regions.txt` (`NAME chr:start-end`, up to 50 kb) is one more
   target, handled exactly like an exon: the region is both the target and the "exon".

## Where the amplicon may sit

The design is meant for **genomic DNA** — confirming an exon deletion or duplication, for instance
— so what matters is that the amplicon interrogates the target:

- the amplicon must cover at least **40 target bases** (`--min-overlap`), or the **whole target**
  when it is shorter than that;
- it is first sought **entirely inside the exon** (placement `exonic`);
- only if no pair fits there may the primers reach up to **100 bp into the flanking intron**
  (`--flank`; placement `spanning`). An exon shorter than the minimum amplicon, 50 bp, is always
  spanned, and the amplicon then contains all of it.

A spanning amplicon is not a warning — it is noted, because a deletion that removes only the exon
and not the primer sites is still detected, but a breakpoint inside the flank would be read
differently.

## Primer rules

Primer3 is run under the Primer Express quantification guidelines, as hard constraints:

| | value | flag |
|---|---|---|
| primer Tm | 58–60 °C, opt 59 | `--tm` |
| primer length | 18–28 nt, opt 20 | `--size` |
| primer GC | 30–80 % | `--gc` |
| 3′ end | ≤ 2 G/C among the last five bases | `--max-end-gc` |
| runs | no run of 4 identical bases | `--max-poly-x` |
| pair ΔTm | ≤ 2 °C | `--max-tm-diff` |
| amplicon | 50–150 bp, 80 preferred | `--product MIN OPT MAX` |

Primer3's thermodynamic checks for hairpins, self-dimers and primer-dimers stay at its defaults;
they matter more in SYBR Green than anywhere, since any dimer fluoresces. Primer Express sets no
GC clamp, and a clamp would fight the 3′ G/C rule, so none is applied.

Pairs come back ranked by **Primer3's penalty** — how far each primer is from the optimum Tm and
length, the pair's ΔTm, and the distance from the preferred amplicon size. Alternatives must differ
from each other by at least 3 bases at their 3′ ends, so they are genuinely different pairs. See
[Primer Express](primer-express.md) for why this is not Primer Express's penalty.

If nothing is found, the design is retried once with relaxed settings (Tm 57–61 °C, ≤ 3 G/C at
the 3′ end, runs of 4, amplicon up to 200 bp), and the result is flagged.

## Safeguards

**1. Common SNPs under primers.** Common variants from UCSC `dbSnp155Common` (MAF ≥ 1 % in
1000 Genomes; `--snp-maf 0.05` to raise the cut-off) are masked before Primer3. On genomic DNA
this is the safeguard that matters most for copy-number work: a heterozygous SNP under a primer
amplifies one allele less well and reads like a heterozygous deletion.

**2. Repeats.** Ensembl soft-masking (RepeatMasker, TRF, Dust) plus UCSC `rmsk`.

Design runs in tiers and stops at the first that works:

| Tier | Accepted | Status |
|---|---|---|
| `clean` | no common SNP and no repeat anywhere under either primer | OK |
| `repeat-5p` | a repeat touches only the 5′ part of a primer (3′ end in unique sequence) | Check |
| `snp-5p` | a common SNP lies under a primer, but never in its last 5 bases (`--snp-3prime`) | Check |

The full order of preference is: Primer Express settings before relaxed; `clean` before
`repeat-5p` before `snp-5p`; inside the exon before spanning. So a clean pair that reaches into
the intron is preferred to an exonic pair with a SNP under it. `--strict` keeps only `clean`.

SNPs **between** the primers are listed as a note, not a warning: they do not affect
amplification, only — slightly — the melt curve.

**3. High homology.**
- Targets overlapping a segmental duplication (UCSC `genomicSuperDups`, ≥ 90 % identity, ≥ 1 kb)
  are flagged with the identity and the duplicated locus. For copy-number work this is the case
  where a specific in-silico PCR result is least reassuring: a mismatched paralog can still
  amplify.
- With `--ispcr` and `--genome-2bit`, every candidate pair (best and alternatives) is run through
  UCSC **isPcr** genome-wide, with the settings of the UCSC In-Silico PCR web tool (15 perfect
  3′ bases, products up to 4 kb). Specific pairs are moved first; off-target products are listed.
  Hits on alt/fix haplotype contigs are ignored. The GitHub workflow downloads isPcr and the
  genome automatically.

## Status

An assay is **OK** only if it is `clean`, designed within the Primer Express rules, outside
segmental duplications and — when checked — gives exactly one in-silico PCR product, at the
target. Otherwise it is **Check**, with the reasons listed. A target with no pair at all is
**FAIL**, with Primer3's own explanation of what it rejected.
