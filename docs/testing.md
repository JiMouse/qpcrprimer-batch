# Tests

What each suite pins down. [← back to the README](../README.md) Both run in the workflow before
anything is designed.

`python tests/test_offline.py` builds a synthetic genome with genes on both strands and fake
Ensembl / UCSC / isPcr responses (the UCSC `mane` track and the batched sequence endpoint
included, so the fast paths are the ones under test). It checks:

- that the defaults **are** the Primer Express rules, and that every primer of every candidate
  pair obeys them — Tm 58–60 °C, GC 30–80 %, ≤ 2 G/C among the last five bases, no run of 4 G,
  ΔTm ≤ 2 °C, amplicon 50–150 bp;
- the Primer Express penalty, pinned on sets of a real export, and equal to Primer3's penalty for
  every candidate; the Primer Express Tm formula, equal to the Tm of every candidate; the amplicon
  size bands; the penalty cut-off flag;
- hairpin and primer-dimer refusal, and that no candidate carries one;
- that each primer matches the genome at its reported coordinates, on both strands, and that the
  published amplicon is that genomic sequence;
- placement: a long exon gets an amplicon inside it; a 31-bp exon is spanned and covered entirely;
  an exon with a long UTR is covered on its coding part, not in the UTR; a short exon whose exonic
  pair would cross the penalty cut-off is spanned instead; the declared placement and coverage
  match the coordinates;
- SNPs: a SNP under the best pair moves the design to another pair; dense SNPs fall back to the
  `snp-5p` tier with no SNP in the last 5 bases and a Check status; `--strict` refuses instead; the
  MAF cut-off; a SNP between the primers is a note, not a warning;
- repeats avoided, or accepted only in a primer's 5′ part when nothing else fits;
- segmental-duplication alerts; an off-target in-silico PCR product moves a pair down, and a
  target with no specific pair is Check;
- custom regions: parsing (commas, reversed ends, chrM, the 50-kb limit, duplicate names), design
  inside the region, a short region spanned like a short exon, a name clashing with a gene refused;
- the outputs: Excel sheets and columns, JSON index with kind and span, the command line end to
  end, the MANE-track model identical to the Ensembl one, one batched sequence request, the merge
  of gene chunks and the regions job, the Tm UCSC values pinned against what UCSC prints, both Tm formulas
  reaching Primer3, and the Primer Express comparison script on an export and on a plain list.

`node tests/test_page.mjs` runs the whole `<script>` of `site/index.html` against a stubbed DOM and
a stubbed `fetch` serving synthetic genes and a custom region. It covers the exon chips and their
status, the assay panel (names, both primers, product Tm, penalty, placement, 3′ G/C, alternatives,
UCSC link), that warnings are HTML-escaped, the region search — a range inside one exon, over a
whole gene, a point in an intron (the exons either side), a minus-strand gene, a custom region, a
range outside the build — and that genes on other chromosomes are never fetched, routing of
symbols, region names and coordinates, and the exports: order rows (`GENE-qE<n>-F`, seven
tab-separated columns in the order of the synthesis order sheet) and the TSV.
