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
- it is first sought **entirely inside the exon** (placement `exonic`), and kept there while its
  penalty stays under 100;
- otherwise the primers may reach up to **100 bp into the flanking intron** (`--flank`; placement
  `spanning`), and the lower-penalty pair of the two searches wins. An exon shorter than the
  minimum amplicon, 50 bp, is always spanned, and the amplicon then contains all of it.

A spanning amplicon is not a warning — it is noted, because a deletion that removes only the exon
and not the primer sites is still detected, but a breakpoint inside the flank would be read
differently.

## Primer rules and penalty

Primer3 is run to reproduce Primer Express — its rules, its Tm and its penalty, measured on its
output (see [Primer Express](primer-express.md)):

| | value | flag |
|---|---|---|
| primer Tm (Primer Express formula) | 58–60 °C | `--tm`, `--tm-formula` |
| primer length | 16–30 nt, 20 costs nothing | `--size` |
| primer GC | 30–80 % | `--gc` |
| 3′ end | ≤ 2 G/C among the last five bases | `--max-end-gc` |
| runs | no run of 4 G; other runs up to 5 | `--max-poly-g`, `--max-poly-x` |
| pair ΔTm | ≤ 2 °C | `--max-tm-diff` |
| amplicon | 50–150 bp | `--product MIN OPT MAX` |
| **penalty** | **5 × (amplicon − 50) + \|F length − 20\| + \|R length − 20\|**; a best pair at 100 or more is flagged | `--penalty-max` |

The penalty is Primer3's own, with weights set to give exactly Primer Express's number, so the
shortest amplicon wins: with these rules most exons get a 50–55 bp amplicon. Alternatives must
differ from each other by at least 3 bases at the 3′ end of each primer, so they are genuinely
different pairs.

**Hairpins and primer-dimers** are refused as Primer3 refuses them — any hairpin, self-dimer or
primer-dimer, anywhere or anchored on a 3′ end, whose Tm exceeds 47 °C — computed in a PCR buffer
(50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer). In SYBR Green any dimer fluoresces, so this
is checked on every candidate. Primer3 searches pairs without these checks, one amplicon-size band
at a time from the shortest, and they are applied to its candidates afterwards with the same
thermodynamic model — many times faster than inside its search. Primer3 may reuse a primer across
the pairs it returns, so a good primer whose first partner fails a check is still paired with
others; the "3 bases apart" rule for alternatives is applied after the checks. Primer Express sets no GC clamp, and a clamp would fight the 3′ G/C rule, so none is
applied.

If nothing is found, the design is retried once with relaxed settings (Tm 57–61 °C, ≤ 3 G/C at
the 3′ end, runs one base longer, amplicon up to 200 bp, primers up to 150 bp into the intron), and
the result is flagged. The extra reach is for short exons between AT-rich introns — PTEN exon 3,
45 bp, has no pair within 100 bp on either side.

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
`repeat-5p` before `snp-5p`; then inside the exon before spanning, as long as the exonic pair's
penalty stays under 100. So a clean pair that reaches into the intron is preferred to an exonic
pair with a SNP under it. `--strict` keeps only `clean`.

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

An assay is **OK** only if it is `clean`, designed within the Primer Express rules with a penalty
under 100, outside segmental duplications and — when checked — gives exactly one in-silico PCR
product, at the target. Otherwise it is **Check**, with the reasons listed. A target with no pair at all is
**FAIL**, with Primer3's own explanation of what it rejected.
