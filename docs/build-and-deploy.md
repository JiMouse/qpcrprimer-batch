# Building and publishing

Local options beyond the quickstart, custom regions, and how the sharded GitHub workflow publishes
the site. [← back to the README](../README.md)

## Local options

Ensembl and UCSC responses are cached in `.cache/`, so reruns with new parameters are fast. UCSC
is queried at most once per second (`--ucsc-delay`), the rate UCSC asks for. `--assembly GRCh37`
uses the GRCh37 Ensembl server and hg19 annotations (no MANE there, canonical is used).

## Custom regions

`regions.txt` holds any target that is not a coding exon — an intronic or intergenic region, a UTR,
a CNV breakpoint:

```
MY_REGION   chr7:117,559,590-117,559,800
```

The name becomes the primer name (`MY_REGION_qF` / `_qR`, ordered as `MY_REGION-q-F`) and the
search key on the site; it may not repeat a gene symbol. Push the file and the workflow designs the
regions in a job of their own; locally, `--regions regions.txt` works with or without a gene list.
Coordinates are on the build's assembly (GRCh38 by default).

## Genome-wide in-silico PCR locally

isPcr lives in the `blat/` subdirectory of the UCSC binaries; replace `linux.x86_64` with
`macOSX.arm64` (Apple Silicon) or `macOSX.x86_64` (Intel):

```bash
mkdir -p tools
curl -o tools/isPcr https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/isPcr && chmod +x tools/isPcr
curl -o tools/hg38.2bit https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit   # ~800 MB
python qpcrprimer_batch.py genes.txt --ispcr tools/isPcr --genome-2bit tools/hg38.2bit
```

isPcr indexes the whole genome in memory (several GB of RAM). If it fails, the run continues and
the output says "in-silico PCR not run". On macOS the downloaded binary is unsigned, so Gatekeeper
may need `xattr -d com.apple.quarantine tools/isPcr`.

## How long a build takes

The requests to Ensembl and UCSC are exactly those of exonprimer-batch — the MANE track, batched
sequence and three UCSC annotation calls per gene — so the waiting is the same; its
[build page](https://github.com/JiMouse/exonprimer-batch/blob/main/docs/build-and-deploy.md)
has the measurements. The design itself costs more here, because each exon is tried inside the
exon before the flanks and the 58–60 °C window is narrow: about **0.7 s per gene** measured on
150 random genes of the panel, so roughly 80 minutes of computation for the whole list, spread
over the parallel jobs.

## Publish the website on GitHub

1. Create the repository `qpcrprimer-batch` and push this folder to its `main` branch.
2. **Settings → Pages → Build and deployment → Source: GitHub Actions.**
3. Run the workflow from the **Actions** tab (or push a change to `genes.txt` / `regions.txt`).
4. The site appears at `https://<user>.github.io/qpcrprimer-batch/`. Links like `…/#PTEN/5` open
   a gene and exon directly, and `…/#q=chr10:87863000-87971000` opens a region search. The full
   Excel is downloadable from the site and from each workflow run's artifacts.

The workflow runs the offline and website tests, designs primers, deploys `site/` and commits the
built data back to `main` as `Built data: N files changed`, so `git log -p site/data/genes/PTEN.json`
is the history of that gene's primers. The workbook is not committed — it is a build artifact and
a download on the site.

### Large gene lists

A few thousand genes do not fit in one GitHub job (6 h limit), so the build is **sharded**:

| Job | What it does |
|---|---|
| `plan` | runs the tests, strips comments from `genes.txt`, `split`s it into chunks, and adds one job for `regions.txt` when it lists any region |
| `design` | a matrix job per chunk, each designing its own genes (or the regions) and uploading `part.xlsx` + `data/` |
| `merge` | `scripts/merge_parts.py` stitches the chunks into one `site/` and one workbook |
| `deploy` | publishes the merged `site/` to Pages |

Tune it from the **Actions** tab: `chunk_size` (genes per job, default 150) and `max_parallel`
(concurrent jobs, default 4). **Keep `max_parallel` low** — every job queries Ensembl and UCSC,
and parallel jobs multiply that load on public servers.

### Changing the website without rebuilding the data

The built data is also published as an asset on the `data-cache` prerelease, named after a
fingerprint of the inputs that determine it — `genes.txt`, `regions.txt`, `qpcrprimer_batch.py`,
`scripts/merge_parts.py`, the assembly, the isPcr flag and any extra options. The website is
deliberately **not** part of that fingerprint: editing `site/index.html` redeploys in a couple of
minutes from the published data, without designing anything. Tick **force_rebuild** in the
Actions tab to rebuild regardless. The three newest fingerprints are kept; a partial build is never
published or deployed.

Expect retired gene symbols to fail: Ensembl does not follow HGNC aliases, so they land in the
Problems sheet. Use current symbols.
