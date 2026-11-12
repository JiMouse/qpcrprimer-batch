# qpcrprimer-batch

### → **[jimouse.github.io/qpcrprimer-batch](https://jimouse.github.io/qpcrprimer-batch)**

qPCR (SYBR Green) primer pairs for **every coding exon of a list of genes**, and for any
genomic region you list, designed under the **Primer Express** rules, SNP-masked, checked
genome-wide by **in-silico PCR**, exported to **Excel** and published as a **GitHub Pages site**.
Search a gene (`PTEN`), a custom region by name, or paste a region (`chr10:87,863,000-87,971,000`,
for instance a CNV call) to list every assay inside it.

Companion to [exonprimer-batch](https://github.com/JiMouse/exonprimer-batch), which designs
the Sanger primers; same gene list, same transcripts, same architecture.

## From the manual workflow to this tool

| Step | By hand | Here |
|---|---|---|
| 1. Define the target | pick a region | each **coding exon** of the MANE Select transcript, or a `chr:start-end` from [`regions.txt`](regions.txt) |
| 2. Extract the sequence | copy from a browser | Ensembl, soft-masked, batched and cached |
| 3. Design | Primer Express | Primer3 with the **Primer Express quantification rules as hard constraints** |
| 4. Rank by penalty | PE penalty < 100 | Primer3 penalty, specific pairs first — see [Primer Express](docs/primer-express.md) for why the two scales differ |
| 5. Check SNPs | dbSNP by eye | common SNPs (dbSNP 155, MAF ≥ 1 %) **masked before design**; any left under a primer is flagged |
| 6. Check uniqueness | UCSC In-Silico PCR | **isPcr genome-wide** on every candidate pair, plus segmental-duplication flags |

## Documentation

| | |
|---|---|
| **[Design method](docs/method.md)** | the rules, where the amplicon may sit, the SNP / repeat / homology tiers, and what makes an assay OK |
| **[Primer Express](docs/primer-express.md)** | which Primer Express rules are applied and how, why its penalty score cannot be reproduced, and how to calibrate the Tm against primers it designed |
| **[Tm](docs/tm.md)** | the three temperatures reported (primer Tm, Tm UCSC, product Tm) and which to use |
| **[Building and publishing](docs/build-and-deploy.md)** | local options, custom regions, the sharded GitHub workflow and the data cache |
| **[Tests](docs/testing.md)** | what each suite pins down |

## In short

For each gene in `genes.txt` (the exonprimer-batch panel, 6,437 genes) and each region in
`regions.txt`:

- **one target per coding exon** (the coding part by default, `--target full` for the whole exon);
- primers **18–28 nt (opt 20), Tm 58–60 °C (opt 59), GC 30–80 %, ≤ 2 G/C among the last five
  3′ bases, no run of 4 identical bases, ΔTm ≤ 2 °C**; amplicon **50–150 bp**, shorter preferred;
- the amplicon covers at least 40 target bases, or the whole exon when it is shorter; it is
  placed **inside the exon** whenever possible and reaches into the intron (≤ 100 bp) only when it
  must — a 31-bp exon is spanned, so all of it is interrogated;
- common SNPs and repeats are masked **before** design: on genomic DNA, a SNP under a primer can
  drop one allele and **mimic a heterozygous deletion**;
- every candidate pair goes through genome-wide in-silico PCR; specific pairs are preferred.

An assay is **OK** only if it is SNP- and repeat-free under both primers, designed within the
Primer Express rules, outside segmental duplications and specific. Otherwise it is **Check**, with
the reasons listed. Every rule has a flag that changes it — see [Design method](docs/method.md).

## Run locally

```bash
pip install -r requirements.txt
python qpcrprimer_batch.py genes.txt -o qpcr_primers.xlsx                      # Excel only
python qpcrprimer_batch.py --regions regions.txt -o my_regions.xlsx             # custom regions only
python qpcrprimer_batch.py genes.txt --regions regions.txt \
       -o site/qpcr_primers.xlsx --json-dir site/data                          # Excel + website data
python -m http.server -d site 8000                                             # preview the site
python qpcrprimer_batch.py --help                                              # all options
```

A small gene list is the quickest way to try it: `echo PTEN > mine.txt` and pass `mine.txt`.
Ensembl and UCSC responses are cached in `.cache/`, so reruns with new parameters are fast.
Genome-wide in-silico PCR and the GitHub workflow are in
[Building and publishing](docs/build-and-deploy.md).

## Excel output

| Sheet | Content |
|---|---|
| Primers | One row per exon or region: status (OK / CHECK / FAIL), names (`GENE_qE3F/R`), sequences, length, Tm, Tm UCSC, GC, **3′ G/C**, ΔTm, **amplicon size and product Tm**, Primer3 penalty, placement (inside the exon or spanning), target bases covered, coordinates, design tier, SNPs under F, under R and between the primers, repeats, segmental duplication, in-silico PCR result, warnings, notes, UCSC link |
| Alternatives | The best pair plus 2 alternatives per target, each with its penalty, SNPs and in-silico PCR result |
| Amplicons | Amplicon sequences, exon UPPER / intron lower |
| Genes | Transcript, RefSeq, targets, OK / to check / failed per gene or region |
| Problems | Genes not found, targets without primers |
| Parameters | All settings used |

Primer names carry a `q` — `PTEN_qE5F`, ordered as `PTEN-qE5-F` — so a qPCR primer never shares
a name with the Sanger primer of the same exon (`PTEN_E5F`).

## Before ordering primers

- Review every **Check** row: a SNP under a primer, a repeat, a segmental duplication or an
  off-target product each has its own reason line.
- A SNP **between** the primers (listed in Notes) does not affect amplification, but it can shift
  the melt curve of a heterozygous sample.
- For genes with known pseudogenes, a specific in-silico PCR result does not rule out mismatched
  priming on the paralog — compare with a second assay, or design against a paralog-specific
  difference.
- Check a new assay's efficiency (standard curve) and its melt curve (one peak, near the product
  Tm shown) before using it.

## Development

Developed with AI assistance (Claude Code). The method, parameters and clinical interpretation
are the author's; every primer in the output should be reviewed before ordering, as above.
