#!/usr/bin/env python3
"""
qpcrprimer_batch.py - qPCR (SYBR Green) primer pairs for every coding exon of
a list of genes, and for any genomic region, with Excel and JSON (website)
output.

The workflow it automates is the manual one:
  1. define the target region          -> each coding exon (MANE Select), or a
                                          chr:start-end from --regions
  2. extract the DNA sequence          -> Ensembl, soft-masked, batched
  3. design in Primer Express          -> Primer3 under Primer Express's
                                          quantification rules, Tm formula and
                                          penalty (below)
  4. keep the lowest-penalty pairs     -> ranked on Primer Express's penalty;
                                          a best pair at 100 or more is flagged
  5. check for SNPs                    -> dbSNP 155 Common masked before design
  6. check uniqueness (In-Silico PCR)  -> UCSC isPcr genome-wide on every
                                          candidate; segmental duplications
                                          flagged

Primer Express quantification guidelines, enforced as hard rules:
  primer Tm 58-60 C, GC 30-80 %, no more than two G/C among the last five 3'
  bases, no run of four G (runs of up to five other bases, as Primer Express
  accepts), amplicon 50-150 bp.
Measured on Primer Express output and reproduced here:
  Tm      = Breslauer 1986 nearest-neighbour, Schildkraut salt term,
            50 mM Na+, 5 nM oligo
  penalty = 5 x (amplicon - 50) + |F length - 20| + |R length - 20|

Placement (genomic DNA, e.g. exon copy-number confirmation): the amplicon
must cover at least --min-overlap bases of the target, or the whole target
when it is shorter. It is first sought entirely inside the exon; only when
that fails may the primers reach up to --flank bp into the flanking intron.

Usage:
  python qpcrprimer_batch.py genes.txt -o qpcr_primers.xlsx --json-dir site/data \
      [--regions regions.txt] [--ispcr tools/isPcr --genome-2bit tools/hg38.2bit]

genes.txt: one HGNC symbol per line; optionally a transcript ID after it
(e.g. "BRCA1  ENST00000357654") to override the MANE Select choice.
regions.txt: one region per line, "NAME chr7:117559590-117559800" (the name
is optional).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlencode

from functools import lru_cache

import primer3
import requests

__version__ = "1.0.0"

SERVERS = {
    "GRCh38": "https://rest.ensembl.org",
    "GRCh37": "https://grch37.rest.ensembl.org",
}
UCSC_DB = {"GRCh38": "hg38", "GRCh37": "hg19"}
UCSC_API = "https://api.genome.ucsc.edu"
SNP_TRACK, RMSK_TRACK, SEGDUP_TRACK = "dbSnp155Common", "rmsk", "genomicSuperDups"
MANE_TRACK = "mane"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
TIERS = {
    "clean": "no common SNP and no repeat under either primer",
    "repeat-5p": "a repeat overlaps the 5' part of a primer (3' end outside repeats)",
    "snp-5p": "a common SNP lies under a primer, outside its 3'-terminal bases",
}
PLACEMENTS = {
    "exonic": "amplicon entirely inside the exon (or the region)",
    "spanning": "primers reach into the flanking sequence, where the lowest-penalty pair was",
}
MAX_REGION = 50_000      # bp; a custom region larger than this is refused


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# Secondary structures are judged as Primer3 judges them by default - the Tm
# of the structure, at most 47 C - but in a PCR buffer, whatever Tm formula
# selects the primers: 50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM oligo.
STRUCT_MAX_TM = 47.0
STRUCT_COND = {"mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.6, "dna_conc": 50.0}


@lru_cache(maxsize=200_000)
def oligo_ok(seq: str) -> bool:
    """No hairpin, self-dimer or 3'-end self-dimer above STRUCT_MAX_TM."""
    return (primer3.calc_hairpin(seq, **STRUCT_COND).tm <= STRUCT_MAX_TM
            and primer3.calc_homodimer(seq, **STRUCT_COND).tm <= STRUCT_MAX_TM
            and primer3.calc_end_stability(seq, seq, **STRUCT_COND).tm <= STRUCT_MAX_TM)


@lru_cache(maxsize=200_000)
def pair_ok(f: str, r: str) -> bool:
    """Both primers clean, and no primer-dimer (any, or 3'-anchored) above STRUCT_MAX_TM."""
    return (oligo_ok(f) and oligo_ok(r)
            and primer3.calc_heterodimer(f, r, **STRUCT_COND).tm <= STRUCT_MAX_TM
            and primer3.calc_end_stability(f, r, **STRUCT_COND).tm <= STRUCT_MAX_TM
            and primer3.calc_end_stability(r, f, **STRUCT_COND).tm <= STRUCT_MAX_TM)


def end_gc(seq: str, n: int = 5) -> int:
    """G/C count among the n 3'-terminal bases (Primer Express allows <= 2)."""
    return sum(c in "GCgc" for c in seq[-n:])


def max_run(seq: str, base: str | None = None) -> int:
    """Longest run of one identical base (of `base` only, if given)."""
    best = run = 0
    prev = ""
    for c in seq.upper():
        run = run + 1 if c == prev else 1
        prev = c
        if base is None or c == base:
            best = max(best, run)
    return best


def tm_ucsc(seq: str) -> float:
    """The Tm UCSC in-silico PCR prints, reproduced exactly.

    hgPcr uses the Breslauer 1986 tables with the Schildkraut-Lifson salt term
    at 50 mM Na+, 50 nM oligo and no magnesium (Primer3's
    PRIMER_TM_FORMULA=0 / PRIMER_SALT_CORRECTIONS=0). It reads 5-8 C above
    the design Tm on a 20-mer; it is reported only so that a pair pasted into
    UCSC In-Silico PCR shows a number you recognise. Same function and same
    pinned values as exonprimer-batch.
    """
    return round(primer3.calc_tm(seq.upper(), mv_conc=50.0, dv_conc=0.0,
                                 dntp_conc=0.0, dna_conc=50.0,
                                 tm_method="breslauer",
                                 salt_corrections_method="schildkraut"), 1)


def ucsc_chrom(chrom: str) -> str:
    return "chrM" if chrom in ("MT", "M") else f"chr{chrom}"


class _Http:
    def __init__(self, cache_dir: str | None, delay: float, name: str):
        self.session = requests.Session()
        self.session.headers["User-Agent"] = f"qpcrprimer-batch/{__version__}"
        self.cache_dir = Path(cache_dir) / name if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.delay = delay

    def post_json(self, path_or_url: str, body: dict, key: str, params=None):
        url = path_or_url if path_or_url.startswith("http") else getattr(self, "base", "") + path_or_url
        return self._request("POST", url, key, params, json.dumps(body))

    def get_json(self, url: str, key: str, params=None, headers=None):
        return self._request("GET", url, key, params, None, headers)

    def _request(self, method: str, url: str, key: str, params=None, body=None, headers=None):
        cache = None
        if self.cache_dir:
            cache = self.cache_dir / (re.sub(r"[^A-Za-z0-9._-]+", "_", key)[:200] + ".json")
            if cache.exists():
                return json.loads(cache.read_text())
        for attempt in range(8):
            try:
                hdr = dict(headers or {})
                if body is not None:
                    hdr.setdefault("Content-Type", "application/json")
                    hdr.setdefault("Accept", "application/json")
                r = self.session.request(method, url, params=params, headers=hdr,
                                         data=body, timeout=600 if body else 120)
            except requests.RequestException:
                time.sleep(2 ** attempt)
                continue
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 2 ** attempt)))
                continue
            if r.status_code in (500, 502, 503, 504):
                time.sleep(2 ** attempt)
                continue
            if r.status_code in (400, 404):
                try:
                    msg = r.json().get("error", r.text)
                except ValueError:
                    msg = r.text
                raise LookupError(msg)
            r.raise_for_status()
            data = r.json()
            if cache:
                cache.write_text(json.dumps(data))
            time.sleep(self.delay)
            return data
        raise RuntimeError(f"request failed repeatedly: {url}")


# --------------------------------------------------------------------------
# Ensembl (exon structure + soft-masked sequence)
# --------------------------------------------------------------------------
class Ensembl(_Http):
    def __init__(self, assembly: str = "GRCh38", cache_dir: str | None = None):
        super().__init__(cache_dir, 0.07, "ensembl")  # < 15 requests/s
        self.base = SERVERS[assembly]

    def _get(self, path, params=None):
        return self.get_json(self.base + path, path + "_" + urlencode(params or {}), params,
                             {"Content-Type": "application/json"})

    def gene(self, symbol):
        return self._get(f"/lookup/symbol/homo_sapiens/{symbol}", {"expand": 1, "mane": 1})

    def transcript(self, transcript_id):
        return self._get(f"/lookup/id/{transcript_id}", {"expand": 1, "mane": 1})

    def sequence(self, chrom, start, end, strand):
        return self._get(f"/sequence/region/human/{chrom}:{start}..{end}:{strand}", {"mask": "soft"})["seq"]

    def sequences(self, regions: list[str]) -> dict[str, str]:
        """Soft-masked sequence for up to 50 regions in one POST, keyed by the
        `query` Ensembl echoes, so a reordered reply can never misalign."""
        if not regions:
            return {}
        data = self.post_json("/sequence/region/human", {"regions": regions},
                              "seqbatch_" + hashlib.sha1("|".join(regions).encode()).hexdigest(),
                              {"mask": "soft"})
        return {it["query"]: it["seq"] for it in data if it.get("query") and it.get("seq")}


# --------------------------------------------------------------------------
# UCSC (common SNPs, RepeatMasker, segmental duplications, MANE)
# --------------------------------------------------------------------------
class UCSC(_Http):
    def __init__(self, assembly: str = "GRCh38", cache_dir: str | None = None, delay: float = 1.0):
        super().__init__(cache_dir, delay, "ucsc")
        self.db = UCSC_DB[assembly]

    def track(self, track: str, chrom: str, start: int, end: int) -> list[dict]:
        """Items of a track overlapping chrom:start-end (1-based inclusive)."""
        q = f"genome={self.db};track={track};chrom={chrom};start={start - 1};end={end};maxItemsOutput=-1"
        data = self.get_json(f"{UCSC_API}/getData/track?{q}", f"{self.db}_{track}_{chrom}_{start}_{end}")
        items = data.get(track, [])
        if isinstance(items, dict):
            items = [x for v in items.values() for x in v]
        return items

    def mane(self, chrom: str) -> list[dict]:
        q = f"genome={self.db};track={MANE_TRACK};chrom={chrom};maxItemsOutput=-1"
        data = self.get_json(f"{UCSC_API}/getData/track?{q}", f"{self.db}_{MANE_TRACK}_{chrom}")
        items = data.get(MANE_TRACK, [])
        if isinstance(items, dict):
            items = [x for v in items.values() for x in v]
        return items

    def snps(self, chrom, start, end, min_maf: float) -> list[dict]:
        out = []
        for it in self.track(SNP_TRACK, chrom, start, end):
            maf = max_maf(it.get("minorAlleleFreq"))
            if maf is not None and maf < min_maf:
                continue
            s0, e0 = int(it["chromStart"]), int(it["chromEnd"])
            out.append({"rs": it.get("name", ""), "start": s0 + 1, "end": max(e0, s0 + 1),
                        "maf": None if maf is None else round(maf, 4),
                        "class": it.get("class", "")})
        return out

    def repeats(self, chrom, start, end) -> list[dict]:
        return [{"start": int(it["genoStart"]) + 1, "end": int(it["genoEnd"]),
                 "name": it.get("repName", ""), "class": it.get("repClass", "")}
                for it in self.track(RMSK_TRACK, chrom, start, end)]

    def segdups(self, chrom, start, end) -> list[dict]:
        out = [{"start": int(it["chromStart"]) + 1, "end": int(it["chromEnd"]),
                "other": f'{it.get("otherChrom", "")}:{int(it.get("otherStart", 0)) + 1}-{it.get("otherEnd", "")}',
                "identity": round(100 * float(it.get("fracMatch", 0)), 1)}
               for it in self.track(SEGDUP_TRACK, chrom, start, end)]
        return sorted(out, key=lambda d: (-d["identity"], d["other"], d["start"]))


def max_maf(v) -> float | None:
    """Highest minor-allele frequency across dbSNP frequency sources."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        vals = [float(v)]
    else:
        vals = []
        for x in (v.split(",") if isinstance(v, str) else v):
            try:
                vals.append(float(x))
            except (TypeError, ValueError):
                pass
    vals = [x for x in vals if 0 <= x <= 1]
    return max(vals) if vals else None


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------
@dataclass
class Params:
    # placement
    target: str = "coding"      # coding = coding part of each exon; full = whole exon
    min_overlap: int = 40       # target bases the amplicon must cover (or all, if fewer)
    flank: int = 100            # how far primers may reach beyond the exon, when they must
    # amplicon - Primer Express: 50-150 bp
    product_min: int = 50
    product_max: int = 150
    # penalty - Primer Express's own, measured on its output (docs/primer-express.md):
    #   5 x (amplicon - 50) + |F length - 20| + |R length - 20|
    product_opt: int = 50
    wt_product: float = 5.0     # penalty per bp above product_opt
    wt_size: float = 1.0        # penalty per nt away from size_opt, per primer
    penalty_max: float = 100.0  # the cut-off used with Primer Express; above it the pair is flagged
    # primers - Primer Express quantification guidelines
    tm_min: float = 58.0
    tm_opt: float = 59.0
    tm_max: float = 60.0
    size_min: int = 16
    size_opt: int = 20
    size_max: int = 30
    gc_min: float = 30.0
    gc_max: float = 80.0
    max_end_gc: int = 2         # G/C among the last five 3' bases
    max_poly_g: int = 3         # longest run of G (Primer Express: fewer than four)
    max_poly_x: int = 5         # longest run of any other base (Primer Express accepts five)
    gc_clamp: int = 0           # PE has none; a clamp would fight max_end_gc
    max_tm_diff: float = 2.0
    # Tm - Primer Express's formula, measured on its output: Breslauer 1986 with
    # the Schildkraut-Lifson salt term, 50 mM Na+, 5 nM oligo, no Mg2+ (docs/tm.md).
    # tm_formula "santalucia" uses SantaLucia 1998 and the salt values below as given.
    tm_formula: str = "primer-express"
    salt_mono: float = 50.0
    salt_div: float = 0.0
    dntp: float = 0.0
    dna_conc: float = 5.0
    num_return: int = 3
    relax: bool = True
    # safeguards
    snp_maf: float = 0.01
    snp_3prime: int = 5
    fallback: bool = True
    ispcr_max_size: int = 4000

    def __post_init__(self):
        if not 0 < self.product_min <= self.product_max:
            raise ValueError("product range must be MIN <= MAX (bp)")
        if self.product_min < 2 * self.size_min:
            raise ValueError(f"a {self.product_min}-bp product cannot hold two {self.size_min}-nt primers")
        if not 0 < self.min_overlap <= self.product_min:
            raise ValueError("min_overlap must be between 1 and the minimum product size")
        if self.target not in ("coding", "full"):
            raise ValueError("target is 'coding' or 'full'")
        if not self.tm_min <= self.tm_opt <= self.tm_max:
            raise ValueError("Tm must be MIN <= OPT <= MAX")
        if self.tm_formula not in ("primer-express", "santalucia"):
            raise ValueError("tm_formula is 'primer-express' or 'santalucia'")

    def tm_args(self) -> dict:
        """The Tm settings for primer3.calc_tm, as the design uses them."""
        pe = self.tm_formula == "primer-express"
        return {"mv_conc": self.salt_mono, "dv_conc": self.salt_div, "dntp_conc": self.dntp,
                "dna_conc": self.dna_conc, "tm_method": "breslauer" if pe else "santalucia",
                "salt_corrections_method": "schildkraut" if pe else "santalucia"}

    def penalty(self, product: int, f_len: int, r_len: int) -> float:
        """Primer Express's penalty, which the Primer3 weights below reproduce."""
        return (self.wt_product * abs(product - self.product_opt)
                + self.wt_size * (abs(f_len - self.size_opt) + abs(r_len - self.size_opt)))

    def primer3_global(self, relaxed: bool = False) -> dict:
        g = {
            "PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,
            "PRIMER_OPT_SIZE": self.size_opt, "PRIMER_MIN_SIZE": self.size_min, "PRIMER_MAX_SIZE": self.size_max,
            "PRIMER_OPT_TM": self.tm_opt, "PRIMER_MIN_TM": self.tm_min, "PRIMER_MAX_TM": self.tm_max,
            "PRIMER_MIN_GC": self.gc_min, "PRIMER_MAX_GC": self.gc_max,
            "PRIMER_GC_CLAMP": self.gc_clamp, "PRIMER_MAX_END_GC": self.max_end_gc,
            "PRIMER_PAIR_MAX_DIFF_TM": self.max_tm_diff, "PRIMER_MAX_POLY_X": self.max_poly_x,
            "PRIMER_MAX_NS_ACCEPTED": 0, "PRIMER_EXPLAIN_FLAG": 1,
            # Penalty = Primer Express's: only amplicon size and primer length
            # count; Tm, GC and Tm difference are rules, not scores.
            "PRIMER_WT_TM_LT": 0.0, "PRIMER_WT_TM_GT": 0.0,
            "PRIMER_WT_SIZE_LT": self.wt_size, "PRIMER_WT_SIZE_GT": self.wt_size,
            "PRIMER_PAIR_WT_DIFF_TM": 0.0,
            "PRIMER_PRODUCT_OPT_SIZE": self.product_opt,
            "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": self.wt_product,
            "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": self.wt_product,
            "PRIMER_TM_FORMULA": 0 if self.tm_formula == "primer-express" else 1,
            "PRIMER_SALT_CORRECTIONS": 0 if self.tm_formula == "primer-express" else 1,
            # Primer3 may reuse a primer across the pairs it returns: a good primer
            # whose first partner fails the structure checks must stay available.
            # Distinct alternatives are chosen afterwards (distinct_pairs).
            "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": -1, "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": -1,
            "PRIMER_SALT_MONOVALENT": self.salt_mono, "PRIMER_SALT_DIVALENT": self.salt_div,
            "PRIMER_DNTP_CONC": self.dntp, "PRIMER_DNA_CONC": self.dna_conc,
        }
        if relaxed:
            g.update({
                "PRIMER_MIN_TM": self.tm_min - 1, "PRIMER_MAX_TM": self.tm_max + 1,
                "PRIMER_MIN_SIZE": max(16, self.size_min - 2), "PRIMER_MAX_SIZE": self.size_max + 2,
                "PRIMER_MIN_GC": max(20.0, self.gc_min - 5), "PRIMER_MAX_GC": min(85.0, self.gc_max + 5),
                "PRIMER_MAX_END_GC": self.max_end_gc + 1, "PRIMER_MAX_POLY_X": self.max_poly_x + 1,
                "PRIMER_PAIR_MAX_DIFF_TM": self.max_tm_diff + 1,
            })
        return g

    def product_range(self, relaxed: bool = False) -> tuple[int, int]:
        return self.product_min, self.product_max + (50 if relaxed else 0)

    def product_bands(self, relaxed: bool = False, cap: int | None = None) -> list[list[int]]:
        """Amplicon sizes in bands, shortest first. Primer3 fills the first band
        before trying the next, which matches a penalty driven by amplicon
        size - and is ~25 times faster than one wide range, where primers that
        all score 0 leave Primer3 nothing to prune its pair search with."""
        lo, hi = self.product_range(relaxed)
        if cap is not None:
            hi = min(hi, cap)
        edges = sorted({lo, *(e for e in (lo + 11, lo + 26, lo + 51, lo + 101) if lo < e <= hi), hi + 1})
        return [[a, b - 1] for a, b in zip(edges, edges[1:])]

    def max_g_run(self, relaxed: bool = False) -> int:
        return self.max_poly_g + (1 if relaxed else 0)

    def flank_for(self, relaxed: bool = False) -> int:
        """How far primers may reach beyond the exon: 50 bp more when relaxed,
        for short exons between AT-rich introns."""
        return self.flank + (50 if relaxed else 0)


@dataclass
class Exon:
    number: int
    start: int
    end: int
    cds_start: int | None
    cds_end: int | None
    ensembl_id: str = ""

    @property
    def coding(self) -> bool:
        return self.cds_start is not None


@dataclass
class Target:
    exon: int
    lo: int          # target (coding part, or whole exon), local coordinates
    hi: int
    elo: int         # the whole exon, local coordinates
    ehi: int
    label: str
    name: str        # primer name stem, e.g. MTHFR_qE5


@dataclass
class GeneModel:
    gene: str
    gene_id: str
    transcript_id: str
    transcript_version: str
    refseq: str
    selection: str
    chrom: str
    strand: int
    exons: list[Exon] = field(default_factory=list)
    kind: str = "gene"                              # gene | region
    region_start: int = 0
    region_end: int = 0
    seq: str = ""                                   # transcript-oriented, soft-masked
    snps: list[dict] = field(default_factory=list)
    repeats: list[dict] = field(default_factory=list)
    segdups: list[dict] = field(default_factory=list)
    annotated: dict = field(default_factory=dict)
    rep: bytearray = field(default_factory=bytearray)
    snp_at: dict = field(default_factory=dict)

    def to_local(self, g: int) -> int:
        return g - self.region_start if self.strand == 1 else self.region_end - g

    def to_genomic(self, i: int) -> int:
        return self.region_start + i if self.strand == 1 else self.region_end - i

    def local_span(self, a: int, b: int) -> tuple[int, int]:
        x, y = self.to_local(a), self.to_local(b)
        return (x, y) if x <= y else (y, x)

    def clip_span(self, a: int, b: int) -> tuple[int, int] | None:
        a, b = max(a, self.region_start), min(b, self.region_end)
        return self.local_span(a, b) if a <= b else None

    def build_masks(self):
        self.rep = bytearray(1 if c.islower() else 0 for c in self.seq)
        for r in self.repeats:
            sp = self.clip_span(r["start"], r["end"])
            if sp:
                self.rep[sp[0]:sp[1] + 1] = b"\x01" * (sp[1] - sp[0] + 1)
        self.snp_at = {}
        for s in self.snps:
            sp = self.clip_span(s["start"], s["end"])
            if sp:
                for i in range(sp[0], sp[1] + 1):
                    self.snp_at.setdefault(i, []).append(s)


# --------------------------------------------------------------------------
# Gene models (identical to exonprimer-batch) and custom regions
# --------------------------------------------------------------------------
def pick_transcript(gene: dict) -> tuple[dict, str]:
    txs = gene.get("Transcript", [])
    for t in txs:
        if not t.get("Translation"):
            continue
        for m in t.get("MANE", []) or []:
            if m.get("type") == "MANE_Select":
                return t, "MANE Select"
    for t in txs:
        if t.get("is_canonical") and t.get("Translation"):
            return t, "Ensembl canonical"
    coding = [t for t in txs if t.get("Translation")]
    if coding:
        return max(coding, key=lambda t: t["Translation"]["length"]), "longest CDS"
    raise LookupError("no protein-coding transcript")


class ManeTable:
    """MANE Select transcripts for the whole genome from the UCSC `mane`
    track: one request per chromosome instead of one Ensembl lookup per gene."""

    def __init__(self, ucsc, log=print):
        self.by_symbol: dict[str, dict] = {}
        if ucsc is None:
            return
        for chrom in CHROMS:
            try:
                rows = ucsc.mane(f"chr{chrom}")
            except Exception as ex:
                log(f"warning: MANE track for chr{chrom} not loaded ({ex}); those genes use Ensembl")
                continue
            for r in rows:
                if r.get("maneStat") != "MANE Select":
                    continue
                sym = r.get("geneName2") or ""
                if sym:
                    self.by_symbol.setdefault(sym.upper(), r)
        log(f"MANE table: {len(self.by_symbol)} genes from {len(CHROMS)} UCSC requests")

    def model(self, symbol: str, pad: int) -> GeneModel | None:
        r = self.by_symbol.get(symbol.upper())
        if not r:
            return None
        strand = -1 if r.get("strand") == "-" else 1
        cds_lo, cds_hi = r["thickStart"] + 1, r["thickEnd"]
        if cds_hi < cds_lo:
            return None
        st = r["chromStart"]
        sizes = [int(v) for v in str(r["blockSizes"]).rstrip(",").split(",") if v != ""]
        offs = [int(v) for v in str(r["chromStarts"]).rstrip(",").split(",") if v != ""]
        if not sizes or len(sizes) != len(offs):
            return None
        blocks = sorted(((st + o + 1, st + o + sz) for o, sz in zip(offs, sizes)),
                        reverse=(strand == -1))
        exons = []
        for i, (lo, hi) in enumerate(blocks, 1):
            cs, ce = max(lo, cds_lo), min(hi, cds_hi)
            exons.append(Exon(i, lo, hi, cs if cs <= ce else None, ce if cs <= ce else None))
        if not any(e.coding for e in exons):
            return None
        m = GeneModel(gene=symbol, gene_id=str(r.get("geneName", "")), transcript_id=str(r["name"]).split(".")[0],
                      transcript_version=str(r["name"]), refseq=str(r.get("ncbiId", "")),
                      selection="MANE Select", chrom=str(r["chrom"]).replace("chr", ""),
                      strand=strand, exons=exons)
        m.region_start = max(1, min(e.start for e in exons) - pad)
        m.region_end = max(e.end for e in exons) + pad
        return m


def fetch_sequences(ens, models: list[GeneModel], batch: int = 50, log=print):
    """Soft-masked sequence for many models in POST batches; a model whose
    sequence does not come back keeps seq == "" and is reported by the caller."""
    todo = [m for m in models if not m.seq]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        regions = [f"{m.chrom}:{m.region_start}..{m.region_end}:{m.strand}" for m in chunk]
        try:
            got = ens.sequences(regions)
        except Exception as ex:
            log(f"warning: batch sequence fetch failed ({ex}); falling back to one request per region")
            got = {}
            for m, q in zip(chunk, regions):
                try:
                    got[q] = ens.sequence(m.chrom, m.region_start, m.region_end, m.strand)
                except Exception as ex2:
                    log(f"    warning: {m.gene} sequence not loaded ({ex2})")
        for m, q in zip(chunk, regions):
            seq = got.get(q)
            if not seq:
                continue
            m.seq = seq
            if len(m.seq) != m.region_end - m.region_start + 1:
                m.region_end = m.region_start + len(m.seq) - 1
        log(f"  sequences {min(i + batch, len(todo))}/{len(todo)}")


def build_model(ens, symbol: str, transcript_id: str | None, pad: int,
                with_seq: bool = True) -> GeneModel:
    """Exon structure from Ensembl (genes the MANE track cannot answer)."""
    gene = ens.gene(symbol)
    if transcript_id:
        tx = ens.transcript(transcript_id.split(".")[0])
        selection = "user"
        if not tx.get("Translation"):
            raise LookupError(f"{transcript_id} is not protein-coding")
    else:
        tx, selection = pick_transcript(gene)
    refseq = next((m.get("refseq_match", "") for m in tx.get("MANE", []) or []
                   if m.get("type") == "MANE_Select"), "")
    strand = tx["strand"]
    tl = tx["Translation"]
    exons = []
    for i, e in enumerate(sorted(tx["Exon"], key=lambda e: e["start"], reverse=(strand == -1)), 1):
        cs, ce = max(e["start"], tl["start"]), min(e["end"], tl["end"])
        exons.append(Exon(i, e["start"], e["end"], cs if cs <= ce else None, ce if cs <= ce else None, e.get("id", "")))
    m = GeneModel(gene=gene.get("display_name", symbol), gene_id=gene["id"], transcript_id=tx["id"],
                  transcript_version=f'{tx["id"]}.{tx.get("version", "")}'.rstrip("."), refseq=refseq,
                  selection=selection, chrom=str(tx["seq_region_name"]), strand=strand, exons=exons)
    m.region_start = max(1, min(e.start for e in exons) - pad)
    m.region_end = max(e.end for e in exons) + pad
    if with_seq:
        m.seq = ens.sequence(m.chrom, m.region_start, m.region_end, strand)
        if len(m.seq) != m.region_end - m.region_start + 1:
            m.region_end = m.region_start + len(m.seq) - 1
    return m


REGION_RE = re.compile(r"^(?:chr)?(2[0-2]|1[0-9]|[1-9]|[XYM]|MT):([0-9,]+)-([0-9,]+)$", re.I)


def parse_region(text: str) -> tuple[str, int, int]:
    m = REGION_RE.match(text.strip())
    if not m:
        raise ValueError(f"not a region: {text!r} (expected chr7:117559590-117559800)")
    chrom = m.group(1).upper()
    a, b = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if a > b:
        a, b = b, a
    if b - a + 1 > MAX_REGION:
        raise ValueError(f"{text}: {b - a + 1:,} bp is larger than the {MAX_REGION:,}-bp limit")
    return ("MT" if chrom == "M" else chrom), a, b


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "region"


def read_region_list(path: str) -> list[tuple[str, str, int, int]]:
    """'NAME chr:start-end' (either order, name optional) -> (name, chrom, start, end)."""
    out, seen = [], set()
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = re.split(r"[\s;\t]+", line)
        loc = next((x for x in parts if REGION_RE.match(x)), None)
        if loc is None:
            raise ValueError(f"{path}:{n}: no chr:start-end in {line!r}")
        chrom, a, b = parse_region(loc)
        rest = [x for x in parts if x != loc]
        name = safe_name(rest[0]) if rest else f"chr{chrom}_{a}_{b}"
        if name.upper() in seen:
            raise ValueError(f"{path}:{n}: region name {name} used twice")
        seen.add(name.upper())
        out.append((name, chrom, a, b))
    return out


def region_model(name: str, chrom: str, start: int, end: int, pad: int) -> GeneModel:
    """A custom region is a one-'exon' model on the + strand: the region is
    the target, and the same placement rules apply to it."""
    m = GeneModel(gene=name, gene_id="", transcript_id="", transcript_version="", refseq="",
                  selection="custom region", chrom=chrom, strand=1,
                  exons=[Exon(1, start, end, start, end)], kind="region")
    m.region_start = max(1, start - pad)
    m.region_end = end + pad
    return m


def annotate(m: GeneModel, ucsc, p: Params, log=print):
    """Load common SNPs, repeats and segmental duplications for the region."""
    chrom = ucsc_chrom(m.chrom)
    for kind, fn in (("snps", lambda: ucsc.snps(chrom, m.region_start, m.region_end, p.snp_maf)),
                     ("repeats", lambda: ucsc.repeats(chrom, m.region_start, m.region_end)),
                     ("segdups", lambda: ucsc.segdups(chrom, m.region_start, m.region_end))):
        if ucsc is None:
            m.annotated[kind] = False
            continue
        try:
            setattr(m, kind, fn())
            m.annotated[kind] = True
        except Exception as ex:
            m.annotated[kind] = False
            log(f"    warning: {m.gene} {kind} not loaded from UCSC ({ex})")
    m.build_masks()


# --------------------------------------------------------------------------
# Targets: one per coding exon (or the region)
# --------------------------------------------------------------------------
def make_targets(m: GeneModel, p: Params) -> list[Target]:
    out = []
    for e in m.exons:
        if not e.coding:
            continue
        a, b = (e.cds_start, e.cds_end) if p.target == "coding" else (e.start, e.end)
        lo, hi = m.local_span(a, b)
        elo, ehi = m.local_span(e.start, e.end)
        if m.kind == "region":
            out.append(Target(e.number, lo, hi, elo, ehi, "", f"{m.gene}_q"))
        else:
            out.append(Target(e.number, lo, hi, elo, ehi, str(e.number), f"{m.gene}_qE{e.number}"))
    return out


# --------------------------------------------------------------------------
# Primer design for one target
# --------------------------------------------------------------------------
def masked_template(m: GeneModel, w0: int, w1: int, tier: str, snp_3prime: int = 0) -> str:
    """'N' = forbidden under a primer; lowercase = a primer's 3' end may not
    sit here (Primer3 rejects a primer whose 3'-terminal base is lowercase).

    clean      SNPs and repeats are N.
    repeat-5p  SNPs are N; repeats are lowercase, so only a primer's 5' part
               may lie in one.
    snp-5p     SNPs stay as the reference base, but every base within
               snp_3prime - 1 of a SNP is lowercase: no primer can then end
               with a SNP among its snp_3prime terminal bases. Both strands
               share one template, so the mask is symmetric and slightly
               stricter than needed; the exact rule is checked on each pair.
    """
    near = set()
    if tier == "snp-5p" and snp_3prime > 0:
        for i in m.snp_at:
            if w0 - snp_3prime <= i <= w1 + snp_3prime:
                near.update(range(i - snp_3prime + 1, i + snp_3prime))
    out = []
    for i in range(w0, w1 + 1):
        c = m.seq[i].upper()
        rep, snp = m.rep[i], i in m.snp_at
        if snp and tier in ("clean", "repeat-5p"):
            c = "N"
        elif rep:
            c = "N" if tier == "clean" else c.lower()
        if i in near:
            c = c.lower()
        out.append(c)
    return "".join(out)


def needed_overlap(t: Target, p: Params) -> int:
    return min(t.hi - t.lo + 1, p.min_overlap)


def design_target(m: GeneModel, t: Target, p: Params) -> dict:
    """Try, in order of preference, until a pair is found:
    normal settings before relaxed; clean before repeat-5p before snp-5p.
    Within a tier the lowest Primer Express penalty wins, wherever the pair
    sits in the exon and its flanks (up to --flank bp into the intron) - as
    Primer Express ranks a pasted exon with its flanking intron. On equal
    penalty, a pair inside the exon comes first. A custom region is what the
    user pasted: its pair stays inside it, and reaches out only when none fits."""
    need = needed_overlap(t, p)
    n_cand = max(6, 2 * p.num_return)
    result = {"pairs": [], "relaxed": False, "tier": "", "placement": "", "explain": ""}
    tiers = list(TIERS) if p.fallback else ["clean"]

    def search(tier: str, relaxed: bool, extra: int) -> tuple[list, dict]:
        g = p.primer3_global(relaxed)
        pmin, _ = p.product_range(relaxed)
        smax = g["PRIMER_MAX_SIZE"]
        w0 = max(0, t.elo - extra)
        w1 = min(len(m.seq) - 1, t.ehi + extra)
        if w1 - w0 + 1 < pmin:
            return [], {}
        # The left primer may not start past the point where the amplicon
        # could still cover `need` target bases, nor the right primer end
        # before it; exact overlap is checked again on every pair below.
        left_end = min(w1, t.hi - need + smax)
        right_start = max(w0, t.lo + need - smax)
        if left_end - w0 + 1 < g["PRIMER_MIN_SIZE"] or w1 - right_start + 1 < g["PRIMER_MIN_SIZE"]:
            return [], {}
        # Primer3's thermodynamic structure checks, inside its pair search,
        # cost ~50x the search itself once the penalty no longer ranks primers
        # by Tm. The search runs without them and they are applied below, by
        # the same thermodynamic model, to the candidates only. The amplicon
        # size bands are walked here, shortest first, until enough pairs pass.
        g_max = p.max_g_run(relaxed)
        template = masked_template(m, w0, w1, tier, p.snp_3prime)
        ok_region = [[0, left_end - w0 + 1, right_start - w0, w1 - right_start + 1]]
        pairs, res = [], {}
        for band in p.product_bands(relaxed, w1 - w0 + 1):
            gg = dict(g, PRIMER_PRODUCT_SIZE_RANGE=[band],
                      PRIMER_LOWERCASE_MASKING=1, PRIMER_NUM_RETURN=400,
                      PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0,
                      PRIMER_MAX_SELF_ANY=9999, PRIMER_MAX_SELF_END=9999,
                      PRIMER_PAIR_MAX_COMPL_ANY=9999, PRIMER_PAIR_MAX_COMPL_END=9999)
            r = primer3.design_primers({"SEQUENCE_ID": t.name, "SEQUENCE_TEMPLATE": template,
                                        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": ok_region}, gg)
            res = res or r
            got = [pair_record(m, t, r, i, w0) for i in range(r.get("PRIMER_PAIR_NUM_RETURNED", 0))]
            got = [b for b in got if b["target_cov"] >= need
                   and max_run(b["f_seq"], "G") <= g_max and max_run(b["r_seq"], "G") <= g_max
                   and (tier != "snp-5p" or all(s["pos3"] > p.snp_3prime for s in b["f_snps"] + b["r_snps"]))
                   and pair_ok(b["f_seq"], b["r_seq"])]
            pairs = distinct_pairs(sorted(pairs + got, key=lambda b: (b["penalty"], b["placement"] != "exonic")))
            if len(pairs) >= n_cand:
                break
        return pairs, res

    for relaxed in ([False, True] if p.relax else [False]):
        for tier in tiers:
            if m.kind == "region":
                pairs, res = search(tier, relaxed, 0)
                if not pairs:
                    pairs, res = search(tier, relaxed, p.flank_for(relaxed))
            else:
                pairs, res = search(tier, relaxed, p.flank_for(relaxed))
            if not pairs and tier == "clean" and not relaxed and res:
                result["explain"] = "; ".join(
                    f"{k.split('_')[1].lower()}: {res.get(k, '')}"
                    for k in ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN"))
            if pairs:
                pairs = pairs[:n_cand]
                for k, b in enumerate(pairs, 1):
                    b["rank"] = k
                result.update(pairs=pairs, relaxed=relaxed, tier=tier, placement=pairs[0]["placement"])
                return result
    if not result["explain"]:
        result["explain"] = "not enough sequence around the target"
    return result


def distinct_pairs(pairs: list[dict], d: int = 3) -> list[dict]:
    """Keep pairs in order, dropping one whose forward or reverse primer ends
    within d bases of the 3' end of a primer already kept, as Primer3's
    PRIMER_MIN_*_THREE_PRIME_DISTANCE does - but only among pairs that passed
    every check."""
    kept, f_ends, r_ends = [], [], []
    for b in pairs:
        fe, re_ = b["f_end3"], b["r_end3"]
        if any(abs(fe - x) < d for x in f_ends) or any(abs(re_ - x) < d for x in r_ends):
            continue
        kept.append(b)
        f_ends.append(fe)
        r_ends.append(re_)
    return kept


def primer_features(m: GeneModel, lo: int, hi: int, three_prime_hi: bool) -> tuple[list, int, list]:
    """SNPs (distance from the 3' end, 1 = terminal base), repeat bases and
    repeat names under a primer occupying local [lo, hi]."""
    snps, seen = [], set()
    for i in range(lo, hi + 1):
        for s in m.snp_at.get(i, []):
            if s["rs"] in seen:
                continue
            seen.add(s["rs"])
            snps.append({"rs": s["rs"], "maf": s["maf"], "pos3": (hi - i + 1) if three_prime_hi else (i - lo + 1)})
    rep_bases = sum(m.rep[lo:hi + 1])
    names = []
    if rep_bases:
        g1, g2 = sorted((m.to_genomic(lo), m.to_genomic(hi)))
        names = sorted({f'{r["name"]} ({r["class"]})' for r in m.repeats if r["start"] <= g2 and r["end"] >= g1})
    return sorted(snps, key=lambda s: s["pos3"]), rep_bases, names


def pair_record(m: GeneModel, t: Target, res: dict, i: int, w0: int) -> dict:
    ls, llen = res[f"PRIMER_LEFT_{i}"]
    rs, rlen = res[f"PRIMER_RIGHT_{i}"]  # rs = 3'-most base of the right primer
    a_lo, a_hi = w0 + ls, w0 + rs
    g1, g2 = sorted((m.to_genomic(a_lo), m.to_genomic(a_hi)))
    exon_mask = bytearray(a_hi - a_lo + 1)
    for e in m.exons:
        lo, hi = m.local_span(e.start, e.end)
        for j in range(max(lo, a_lo), min(hi, a_hi) + 1):
            exon_mask[j - a_lo] = 1
    amp = "".join(c.upper() if exon_mask[j] else c.lower() for j, c in enumerate(m.seq[a_lo:a_hi + 1]))
    f_lo, f_hi = a_lo, a_lo + llen - 1
    r_lo, r_hi = a_hi - rlen + 1, a_hi
    f_snps, f_rep, f_rnames = primer_features(m, f_lo, f_hi, True)
    r_snps, r_rep, r_rnames = primer_features(m, r_lo, r_hi, False)
    # SNPs between the primers: harmless for amplification, reported for the melt curve
    between, seen = [], set()
    for j in range(f_hi + 1, r_lo):
        for s in m.snp_at.get(j, []):
            if s["rs"] not in seen:
                seen.add(s["rs"])
                between.append({"rs": s["rs"], "maf": s["maf"]})
    f_seq = m.seq[f_lo:f_hi + 1].upper()
    r_seq = rc(m.seq[r_lo:r_hi + 1]).upper()
    f_tm, r_tm = round(res[f"PRIMER_LEFT_{i}_TM"], 1), round(res[f"PRIMER_RIGHT_{i}_TM"], 1)
    return {
        "rank": i + 1,
        "f_seq": f_seq, "r_seq": r_seq,
        "f_end3": f_hi, "r_end3": r_lo,     # local 3'-end positions, for distinct_pairs
        "f_tm": f_tm, "r_tm": r_tm,
        "f_tm_ucsc": tm_ucsc(f_seq), "r_tm_ucsc": tm_ucsc(r_seq),
        "d_tm": round(abs(f_tm - r_tm), 1),
        "f_gc": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1), "r_gc": round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 1),
        "f_end_gc": end_gc(f_seq), "r_end_gc": end_gc(r_seq),
        "f_len": llen, "r_len": rlen,
        "product": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
        "product_tm": round(res.get(f"PRIMER_PAIR_{i}_PRODUCT_TM", 0.0), 1),
        "penalty": round(res[f"PRIMER_PAIR_{i}_PENALTY"], 2),
        "placement": "exonic" if (a_lo >= t.elo and a_hi <= t.ehi) else "spanning",
        "target_cov": max(0, min(a_hi, t.hi) - max(a_lo, t.lo) + 1),
        "chrom": m.chrom, "start": g1, "end": g2,
        "amplicon": amp,
        "f_snps": f_snps, "r_snps": r_snps, "amp_snps": between,
        "f_repeat": f_rep, "r_repeat": r_rep, "repeat_names": sorted(set(f_rnames + r_rnames)),
        "ispcr": {"checked": False},
    }


# --------------------------------------------------------------------------
# Whole-gene driver
# --------------------------------------------------------------------------
def design_gene(m: GeneModel, p: Params, assembly: str) -> dict:
    amplicons, failures = [], []
    for t in make_targets(m, p):
        r = design_target(m, t, p)
        tg1, tg2 = sorted((m.to_genomic(t.lo), m.to_genomic(t.hi)))
        sd = [s for s in m.segdups if s["start"] <= tg2 and s["end"] >= tg1]
        rec = {"label": t.label, "name": t.name, "exons": [t.exon],
               "target_start": tg1, "target_end": tg2, "target_len": t.hi - t.lo + 1,
               "exon_len": t.ehi - t.elo + 1, "min_cov": needed_overlap(t, p), "penalty_max": p.penalty_max,
               "relaxed": r["relaxed"], "tier": r["tier"], "placement": r["placement"],
               "candidates": r["pairs"],
               "segdups": [{"other": o, "identity": i} for o, i in
                           sorted({(s["other"], s["identity"]) for s in sd}, key=lambda x: -x[1])]}
        if r["pairs"]:
            amplicons.append(rec)
        else:
            rec["reason"] = r["explain"]
            failures.append(rec)
    g = {"gene": m.gene, "kind": m.kind, "gene_id": m.gene_id, "transcript": m.transcript_version,
         "refseq": m.refseq, "selection": m.selection, "assembly": assembly, "chrom": m.chrom,
         "strand": m.strand, "exon_count": len(m.exons), "annotated": m.annotated,
         "exons": [{"n": e.number, "id": e.ensembl_id, "start": e.start, "end": e.end,
                    "cds_start": e.cds_start, "cds_end": e.cds_end, "coding": e.coding,
                    "len": e.end - e.start + 1} for e in m.exons],
         "amplicons": amplicons, "failures": failures}
    finalize_gene(g, p)
    return g


def finalize_gene(g: dict, p: Params):
    """Rank pairs (specific first, then Primer3 penalty), keep num_return,
    write status, warnings and notes."""
    for a in g["amplicons"]:
        cands = sorted(a["candidates"], key=lambda b: (ispcr_bad(b["ispcr"]), b["rank"]))
        a["pairs"] = [dict(b, rank=k) for k, b in enumerate(cands[:p.num_return], 1)]
        a["warnings"] = amplicon_warnings(g, a)
        a["notes"] = amplicon_notes(g, a)
        a["status"] = "OK" if not a["warnings"] else "CHECK"


def ispcr_bad(ip: dict) -> bool:
    return bool(ip.get("checked") and (ip.get("off_target") or not ip.get("on_target")))


def amplicon_warnings(g: dict, a: dict) -> list[str]:
    """Anything here demotes the amplicon to CHECK."""
    b, w = a["pairs"][0], []
    for side in ("f", "r"):
        for s in b[f"{side}_snps"]:
            maf = f'MAF {s["maf"]:.2f}' if s["maf"] is not None else "common"
            w.append(f'{side.upper()} primer: SNP {s["rs"]} ({maf}) {s["pos3"]} nt from 3\' end - '
                     f'can mimic a heterozygous deletion')
        if b[f"{side}_repeat"]:
            w.append(f'{side.upper()} primer: {b[f"{side}_repeat"]} nt in repeat (5\' part)')
    if b["repeat_names"]:
        w.append("Repeats under primers: " + ", ".join(b["repeat_names"]))
    for s in a["segdups"][:3]:
        w.append(f'Target in segmental duplication: {s["identity"]}% identity with {s["other"]}')
    ip = b["ispcr"]
    if ip.get("checked"):
        if ip.get("off_target"):
            w.append(f'In-silico PCR: {len(ip["off_target"])} off-target product(s): ' + "; ".join(ip["off_target"][:3]))
        if not ip.get("on_target"):
            w.append("In-silico PCR: expected product not found")
    if b["penalty"] >= a.get("penalty_max", 100):
        w.append(f'Primer Express penalty {b["penalty"]:g} (cut-off {a.get("penalty_max", 100):g}): '
                 f'no pair with a shorter amplicon passed the checks')
    if a["relaxed"]:
        w.append("Relaxed settings (outside the Primer Express rules)")
    for kind in ("snps", "repeats", "segdups"):
        if not g["annotated"].get(kind):
            w.append(f"{kind} not checked (UCSC annotation unavailable)")
    return w


def amplicon_notes(g: dict, a: dict) -> list[str]:
    """Information that does not demote the amplicon."""
    b, n = a["pairs"][0], []
    if b["placement"] == "spanning":
        what = "region" if g["kind"] == "region" else f'exon ({a["exon_len"]} bp)'
        n.append("Amplicon reaches into the flanking sequence: no pair fits inside the region"
                 if g["kind"] == "region" else
                 f"Amplicon reaches into the flanking intron: the lowest-penalty pair for this {what}")
    if b["amp_snps"]:
        n.append(f'{len(b["amp_snps"])} common SNP(s) between the primers (not under them): '
                 + ", ".join(s["rs"] for s in b["amp_snps"][:5]))
    return n


# --------------------------------------------------------------------------
# Genome-wide specificity: UCSC isPcr
# --------------------------------------------------------------------------
def run_ispcr(results: list[dict], ispcr: str, genome: str, p: Params, log=print) -> bool:
    queries, where = [], {}
    for gi, g in enumerate(results):
        for ai, a in enumerate(g["amplicons"]):
            for pi, b in enumerate(a["candidates"]):
                pid = f"q{gi}_{ai}_{pi}"
                queries.append(f"{pid}\t{b['f_seq']}\t{b['r_seq']}")
                where[pid] = b
    if not queries:
        return True
    with tempfile.TemporaryDirectory() as td:
        q, out = Path(td) / "pairs.txt", Path(td) / "out.bed"
        q.write_text("\n".join(queries) + "\n")
        cmd = [ispcr, genome, str(q), str(out), "-out=bed", f"-maxSize={p.ispcr_max_size}",
               "-minPerfect=15", "-minGood=15"]
        log(f"Running in-silico PCR on {len(queries)} primer pairs ...")
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except (OSError, subprocess.CalledProcessError) as ex:
            log(f"warning: isPcr failed, specificity not checked ({getattr(ex, 'stderr', '') or ex})")
            return False
        apply_ispcr(where, out.read_text().splitlines())
    return True


def apply_ispcr(where: dict, bed_lines: list[str]):
    hits: dict[str, list] = {}
    for line in bed_lines:
        f = line.split("\t")
        if len(f) >= 4 and not line.startswith(("#", "track")):
            hits.setdefault(f[3], []).append((f[0], int(f[1]) + 1, int(f[2])))
    for pid, b in where.items():
        chrom = ucsc_chrom(b["chrom"])
        on, off = False, []
        for c, s, e in hits.get(pid, []):
            if c == chrom and s <= b["end"] and e >= b["start"]:
                on = True
            elif "_alt" not in c and "_fix" not in c:  # alt haplotypes duplicate the locus
                off.append(f"{c}:{s}-{e} ({e - s + 1} bp)")
        b["ispcr"] = {"checked": True, "on_target": on, "off_target": off}


# --------------------------------------------------------------------------
# Output helpers
# --------------------------------------------------------------------------
def ucsc_pcr_url(f: str, r: str, assembly: str) -> str:
    return "https://genome.ucsc.edu/cgi-bin/hgPcr?" + urlencode(
        {"db": UCSC_DB[assembly], "wp_target": "genome", "wp_f": f, "wp_r": r,
         "wp_size": 4000, "wp_perfect": 15, "wp_good": 15, "Submit": "submit"})


def snp_text(snps: list[dict]) -> str:
    return "; ".join(f'{s["rs"]} ({s["pos3"]} nt from 3\'' + (f', MAF {s["maf"]:.2f})' if s["maf"] is not None else ")")
                     for s in snps)


def ispcr_text(ip: dict) -> str:
    if not ip.get("checked"):
        return "not run"
    if ip.get("off_target"):
        n = len(ip["off_target"]) + (1 if ip.get("on_target") else 0)
        return f"{n} products: " + "; ".join(ip["off_target"])
    return "specific (1 product)" if ip.get("on_target") else "no product found"


def read_gene_list(path: str) -> list[tuple[str, str | None]]:
    out, seen = [], set()
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = re.split(r"[\s,;]+", line)
        if parts[0].upper() not in seen:
            seen.add(parts[0].upper())
            out.append((parts[0], parts[1] if len(parts) > 1 else None))
    return out


def target_label(g: dict, a: dict) -> str:
    return "region" if g["kind"] == "region" else a["label"]


# --------------------------------------------------------------------------
# Excel output
# --------------------------------------------------------------------------
PRIMER_COLS = ["Gene / region", "Exon", "Status", "Transcript", "RefSeq (MANE)",
               "Forward name", "Forward primer (5'→3')", "F len", "F Tm (°C)", "F Tm UCSC (°C)", "F GC %",
               "F 3' G/C (last 5)",
               "Reverse name", "Reverse primer (5'→3')", "R len", "R Tm (°C)", "R Tm UCSC (°C)", "R GC %",
               "R 3' G/C (last 5)",
               "ΔTm F-R (°C)", "Product (bp)", "Product Tm (°C)", "Penalty (Primer Express)", "Placement",
               "Target bases covered", "Chr", "Amplicon start", "Amplicon end",
               "Design tier", "Common SNPs under F", "Common SNPs under R", "Common SNPs between primers",
               "Repeat bases F/R", "Segmental duplication", "In-silico PCR", "Warnings", "Notes", "UCSC check"]
SEQ_COLS = {"Primers": {7, 14}, "Alternatives": {4, 7}, "Amplicons": {5}}


def write_excel(results: list[dict], errors: list[dict], p: Params, path: str, assembly: str, ispcr_ran: bool):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    font, mono = Font(name="Arial", size=10), Font(name="Consolas", size=10)
    bold = Font(name="Arial", size=10, bold=True, color="FFFFFF")
    link = Font(name="Arial", size=10, color="0563C1", underline="single")
    head_fill = PatternFill("solid", start_color="3B2F5C")
    band = PatternFill("solid", start_color="F2F0F7")
    status_fill = {"OK": PatternFill("solid", start_color="D5EFD9"),
                   "CHECK": PatternFill("solid", start_color="FDEBC8"),
                   "FAIL": PatternFill("solid", start_color="F7D3D0")}
    thin = Border(bottom=Side(style="thin", color="D6D2E0"))

    def header(ws, cols, widths, freeze="C2"):
        ws.append(cols)
        for i, w in enumerate(widths, 1):
            c = ws.cell(row=1, column=i)
            c.font, c.fill = bold, head_fill
            c.alignment = Alignment(vertical="center", wrap_text=True)
            ws.column_dimensions[get_column_letter(i)].width = w
        ws.row_dimensions[1].height = 30
        ws.freeze_panes = freeze

    def style_rows(ws, mono_cols=(), link_cols=(), status_col=None, wrap_cols=()):
        prev, shade = None, False
        for row in ws.iter_rows(min_row=2):
            if row[0].value != prev:
                shade, prev = not shade, row[0].value
            for c in row:
                c.font = mono if c.column in mono_cols else font
                c.border = thin
                if shade:
                    c.fill = band
                if c.column in wrap_cols:
                    c.alignment = Alignment(wrap_text=True, vertical="top")
            for col in link_cols:
                c = row[col - 1]
                if c.value:
                    c.hyperlink, c.value, c.font = c.value, "in-silico PCR", link
            if status_col and row[status_col - 1].value in status_fill:
                row[status_col - 1].fill = status_fill[row[status_col - 1].value]

    ws = wb.active
    ws.title = "Primers"
    header(ws, PRIMER_COLS, [12, 7, 8, 19, 15,
                             16, 30, 6, 9, 11, 7, 9,
                             16, 30, 6, 9, 11, 7, 9,
                             9, 9, 9, 9, 10, 9, 5, 12, 12,
                             10, 26, 26, 26, 9, 30, 30, 60, 50, 13])
    for g in results:
        for a in g["amplicons"]:
            b = a["pairs"][0]
            ws.append([g["gene"], target_label(g, a), a["status"], g["transcript"], g["refseq"],
                       a["name"] + "F", b["f_seq"], b["f_len"], b["f_tm"], b["f_tm_ucsc"], b["f_gc"], b["f_end_gc"],
                       a["name"] + "R", b["r_seq"], b["r_len"], b["r_tm"], b["r_tm_ucsc"], b["r_gc"], b["r_end_gc"],
                       b["d_tm"], b["product"], b["product_tm"], b["penalty"], b["placement"],
                       f'{b["target_cov"]}/{a["target_len"]}', b["chrom"], b["start"], b["end"],
                       a["tier"] + (" (relaxed)" if a["relaxed"] else ""),
                       snp_text(b["f_snps"]), snp_text(b["r_snps"]),
                       ", ".join(s["rs"] for s in b["amp_snps"]),
                       f'{b["f_repeat"]}/{b["r_repeat"]}' if (b["f_repeat"] or b["r_repeat"]) else "",
                       "; ".join(f'{s["identity"]}% {s["other"]}' for s in a["segdups"][:3]),
                       ispcr_text(b["ispcr"]), "\n".join(a["warnings"]), "\n".join(a["notes"]),
                       ucsc_pcr_url(b["f_seq"], b["r_seq"], assembly)])
        for f in g["failures"]:
            ws.append([g["gene"], target_label(g, f), "FAIL", g["transcript"], g["refseq"]]
                      + [""] * (len(PRIMER_COLS) - 8) + ["No primer pair found. " + f.get("reason", ""), "", ""])
    style_rows(ws, mono_cols=SEQ_COLS["Primers"], link_cols=(len(PRIMER_COLS),), status_col=3)
    ws.auto_filter.ref = ws.dimensions

    ws2 = wb.create_sheet("Alternatives")
    header(ws2, ["Gene / region", "Exon", "Rank", "Forward primer (5'→3')", "F Tm (°C)", "F Tm UCSC (°C)",
                 "Reverse primer (5'→3')", "R Tm (°C)", "R Tm UCSC (°C)", "ΔTm F-R (°C)",
                 "Product (bp)", "Product Tm (°C)", "Penalty (Primer Express)", "Placement", "SNPs F", "SNPs R",
                 "Repeat bases F/R", "In-silico PCR", "Amplicon start", "Amplicon end"],
           [12, 7, 6, 30, 9, 11, 30, 9, 11, 9, 9, 9, 9, 10, 24, 24, 9, 28, 12, 12])
    for g in results:
        for a in g["amplicons"]:
            for b in a["pairs"]:
                ws2.append([g["gene"], target_label(g, a), b["rank"], b["f_seq"], b["f_tm"], b["f_tm_ucsc"],
                            b["r_seq"], b["r_tm"], b["r_tm_ucsc"], b["d_tm"],
                            b["product"], b["product_tm"], b["penalty"], b["placement"],
                            snp_text(b["f_snps"]), snp_text(b["r_snps"]),
                            f'{b["f_repeat"]}/{b["r_repeat"]}', ispcr_text(b["ispcr"]), b["start"], b["end"]])
    style_rows(ws2, mono_cols=SEQ_COLS["Alternatives"])
    ws2.auto_filter.ref = ws2.dimensions

    ws3 = wb.create_sheet("Amplicons")
    header(ws3, ["Gene / region", "Exon", "Product (bp)", "Product Tm (°C)",
                 "Amplicon sequence (exon UPPER, intron lower)"], [12, 7, 9, 9, 120])
    for g in results:
        for a in g["amplicons"]:
            b = a["pairs"][0]
            ws3.append([g["gene"], target_label(g, a), b["product"], b["product_tm"], b["amplicon"]])
    style_rows(ws3, mono_cols=SEQ_COLS["Amplicons"])

    ws4 = wb.create_sheet("Genes")
    header(ws4, ["Gene / region", "Kind", "Transcript", "RefSeq (MANE)", "Selection", "Chr", "Strand",
                 "Exons (total)", "Targets", "OK", "To check", "Failed", "Common SNPs in region"],
           [12, 7, 19, 15, 16, 5, 7, 9, 9, 7, 8, 7, 10])
    for g in results:
        st = [a["status"] for a in g["amplicons"]]
        ws4.append([g["gene"], g["kind"], g["transcript"], g["refseq"], g["selection"], g["chrom"],
                    "+" if g["strand"] == 1 else "-", g["exon_count"],
                    len(g["amplicons"]) + len(g["failures"]), st.count("OK"), st.count("CHECK"),
                    len(g["failures"]), g.get("n_snps", "") if g["annotated"].get("snps") else "not checked"])
    style_rows(ws4)

    ws5 = wb.create_sheet("Problems")
    header(ws5, ["Gene / region", "Exon", "Target (bp)", "Problem"], [12, 8, 9, 110])
    for g in results:
        for f in g["failures"]:
            ws5.append([g["gene"], target_label(g, f), f["target_len"], "No primer pair found. " + f.get("reason", "")])
    for e in errors:
        ws5.append([e["gene"], "", "", e["error"]])
    style_rows(ws5)

    ws6 = wb.create_sheet("Parameters")
    header(ws6, ["Setting", "Value", "Meaning"], [26, 18, 95], freeze="A2")
    for row in parameter_rows(p, assembly, ispcr_ran):
        ws6.append(row)
    style_rows(ws6)

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)


MEANING = {
    "target": "coding = the coding part of each exon is the target; full = the whole exon",
    "min_overlap": "The amplicon must cover at least this many target bases (or the whole target if it is shorter)",
    "flank": "Primers may reach this far into the flanking intron (bp); the lowest penalty wins, wherever it sits",
    "product_min": "Minimal amplicon length (bp) - Primer Express: 50",
    "product_max": "Maximal amplicon length (bp) - Primer Express: 150",
    "product_opt": "Amplicon length that costs no penalty (bp) - Primer Express: 50",
    "wt_product": "Penalty per bp above product_opt - Primer Express: 5",
    "wt_size": "Penalty per nt a primer is away from size_opt - Primer Express: 1",
    "penalty_max": "A best pair at or above this penalty is flagged - the cut-off used with Primer Express",
    "tm_min": "Minimal primer Tm (°C) - Primer Express: 58",
    "tm_opt": "Optimal primer Tm (°C) - Primer Express: 59",
    "tm_max": "Maximal primer Tm (°C) - Primer Express: 60",
    "size_opt": "Optimal primer length (nt) - Primer Express: 20",
    "gc_min": "Minimal primer GC % - Primer Express: 30",
    "gc_max": "Maximal primer GC % - Primer Express: 80",
    "max_end_gc": "Most G/C allowed among the five 3'-terminal bases - Primer Express: 2",
    "max_poly_g": "Longest run of G - Primer Express: fewer than four consecutive G",
    "max_poly_x": "Longest run of any base - Primer Express accepts five C or A in a row",
    "tm_formula": "primer-express = Breslauer 1986 + Schildkraut salt term at 50 mM Na+, 5 nM oligo, as "
                  "measured on Primer Express output; santalucia = SantaLucia 1998 with the salt values below",
    "gc_clamp": "Primer3 GC clamp (0: Primer Express has none, and it would fight max_end_gc)",
    "max_tm_diff": "Largest |Tm(F) - Tm(R)| accepted for a pair (°C)",
    "salt_mono": "Monovalent cations for the Tm calculation (mM)",
    "salt_div": "Divalent cations, Mg2+, for the Tm calculation (mM)",
    "dntp": "Total dNTP for the Tm calculation (mM)",
    "dna_conc": "Annealing oligo concentration for the Tm calculation (nM)",
    "num_return": "Pairs kept per target (best + alternatives)",
    "relax": "If no pair is found, retry once with relaxed settings (Tm ±1 °C, 3' G/C ≤ 3, runs one longer, product ≤ 200 bp, "
             "primers up to flank + 50 bp into the intron) - flagged",
    "snp_maf": f"SNPs from UCSC {SNP_TRACK} at or above this MAF are avoided under primers",
    "snp_3prime": "A SNP is never accepted within this many bases of a primer's 3' end",
    "fallback": "Allow tiers repeat-5p and snp-5p when no clean pair exists (always flagged)",
    "ispcr_max_size": "In-silico PCR: largest product considered (bp)",
}


def parameter_rows(p: Params, assembly: str, ispcr_ran: bool) -> list[list]:
    rows = [["Assembly", assembly, "Genome build used for coordinates"],
            ["Transcript choice", "MANE Select", "Falls back to Ensembl canonical; overridable in the gene list"],
            ["Repeats", "masked", f"Ensembl soft-masking (RepeatMasker, TRF, Dust) + UCSC {RMSK_TRACK}"],
            ["Homology", "flagged", f"UCSC {SEGDUP_TRACK} overlap; in-silico PCR " +
             ("run on all candidate pairs (isPcr)" if ispcr_ran else "NOT run - use the UCSC links")]]
    rows += [[k, str(v), MEANING.get(k, "Primer3 setting")] for k, v in vars(p).items()]
    rows += [[f"Tier: {k}", "", v] for k, v in TIERS.items()]
    rows += [[f"Placement: {k}", "", v] for k, v in PLACEMENTS.items()]
    rows += [["Status OK", "", "Clean tier, normal settings, no segmental duplication, no off-target product"],
             ["Penalty", "", "Primer Express's penalty, reproduced: 5 x (amplicon - 50) + |F length - 20| + "
                             "|R length - 20|. Lower is better; pairs are ranked on it after the SNP, repeat "
                             "and in-silico PCR checks"],
             ["Generated", time.strftime("%Y-%m-%d %H:%M"),
              f"qpcrprimer_batch.py v{__version__}, primer3-py {primer3.__version__}"]]
    return rows


def write_json(results, errors, p: Params, out_dir: str, assembly: str, ispcr_ran: bool):
    d = Path(out_dir)
    (d / "genes").mkdir(parents=True, exist_ok=True)
    index = []
    for g in results:
        slim = dict(g, amplicons=[{k: v for k, v in a.items() if k != "candidates"} for a in g["amplicons"]],
                    failures=[{k: v for k, v in f.items() if k != "candidates"} for f in g["failures"]])
        (d / "genes" / f'{g["gene"]}.json').write_text(json.dumps(slim, separators=(",", ":")))
        st = [a["status"] for a in g["amplicons"]]
        index.append({"gene": g["gene"], "kind": g["kind"], "transcript": g["transcript"], "refseq": g["refseq"],
                      "chrom": g["chrom"], "start": min(e["start"] for e in g["exons"]),
                      "end": max(e["end"] for e in g["exons"]),
                      "amplicons": len(st), "ok": st.count("OK"),
                      "check": st.count("CHECK"), "failures": len(g["failures"])})
    meta = {"assembly": assembly, "ucsc_db": UCSC_DB[assembly], "generated": time.strftime("%Y-%m-%d"),
            "version": __version__, "params": vars(p), "ispcr": ispcr_ran, "tiers": TIERS,
            "placements": PLACEMENTS, "genes": index, "errors": errors}
    (d / "index.json").write_text(json.dumps(meta, indent=1))


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def salt_params(salt, formula: str) -> dict:
    if salt is None:
        salt = [50.0, 0.0, 0.0, 5.0] if formula == "primer-express" else [50.0, 1.5, 0.6, 50.0]
    return dict(salt_mono=salt[0], salt_div=salt[1], dntp=salt[2], dna_conc=salt[3])


def main(argv=None, ensembl=None, ucsc=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes", nargs="?", help="text file with one gene symbol per line (optional transcript ID after it)")
    ap.add_argument("--regions", help="text file with one region per line: NAME chr7:117559590-117559800")
    ap.add_argument("-o", "--out", default="qpcr_primers.xlsx", help="Excel output file")
    ap.add_argument("--json-dir", help="also write per-gene JSON for the website here")
    ap.add_argument("--assembly", choices=list(SERVERS), default="GRCh38")
    ap.add_argument("--cache", default=".cache", help="cache directory for web responses ('' to disable)")
    d = Params()
    g1 = ap.add_argument_group("Placement")
    g1.add_argument("--target", choices=["coding", "full"], default=d.target)
    g1.add_argument("--min-overlap", type=int, default=d.min_overlap,
                    help="target bases the amplicon must cover (default 40, or all if fewer)")
    g1.add_argument("--flank", type=int, default=d.flank,
                    help="how far primers may reach into the intron; the lowest penalty wins (default 100)")
    g2 = ap.add_argument_group("Primer Express rules (Primer3)")
    g2.add_argument("--product", type=int, nargs=3, metavar=("MIN", "OPT", "MAX"),
                    default=[d.product_min, d.product_opt, d.product_max],
                    help="amplicon bp; OPT is where the penalty starts (default 50 50 150)")
    g2.add_argument("--tm", type=float, nargs=3, metavar=("MIN", "OPT", "MAX"), default=[d.tm_min, d.tm_opt, d.tm_max])
    g2.add_argument("--size", type=int, nargs=3, metavar=("MIN", "OPT", "MAX"), default=[d.size_min, d.size_opt, d.size_max])
    g2.add_argument("--gc", type=float, nargs=2, metavar=("MIN", "MAX"), default=[d.gc_min, d.gc_max])
    g2.add_argument("--max-end-gc", type=int, default=d.max_end_gc, help="G/C among the last 5 bases (default 2)")
    g2.add_argument("--max-poly-g", type=int, default=d.max_poly_g, help="longest run of G (default 3)")
    g2.add_argument("--max-poly-x", type=int, default=d.max_poly_x, help="longest run of any base (default 5)")
    g2.add_argument("--penalty-max", type=float, default=d.penalty_max,
                    help="flag a best pair at or above this Primer Express penalty (default 100)")
    g2.add_argument("--max-tm-diff", type=float, default=d.max_tm_diff)
    g2.add_argument("--alternatives", type=int, default=d.num_return, help="pairs to keep per target")
    g2.add_argument("--no-relax", action="store_true")
    g2.add_argument("--tm-formula", choices=["primer-express", "santalucia"], default=d.tm_formula,
                    help="primer-express (default) reproduces Primer Express's Tm; santalucia = SantaLucia 1998 "
                         "under --salt (default for it: 50 1.5 0.6 50, a PCR buffer)")
    g2.add_argument("--salt", type=float, nargs=4, metavar=("MONO_MM", "DIV_MM", "DNTP_MM", "DNA_NM"),
                    help="Tm conditions: monovalent mM, divalent mM, dNTP mM, oligo nM "
                         "(default 50 0 0 5 for primer-express)")
    g3 = ap.add_argument_group("SNPs, repeats, uniqueness")
    g3.add_argument("--snp-maf", type=float, default=d.snp_maf, help="avoid SNPs with MAF >= this (default 0.01)")
    g3.add_argument("--snp-3prime", type=int, default=d.snp_3prime, help="no SNP allowed in this many 3' bases")
    g3.add_argument("--strict", action="store_true", help="only accept primers with no SNP and no repeat at all")
    g3.add_argument("--no-ucsc", action="store_true", help="skip UCSC SNP/repeat/segdup annotation")
    g3.add_argument("--ucsc-delay", type=float, default=1.0, help="seconds between UCSC API requests")
    g3.add_argument("--no-mane-track", action="store_true", help="never use the UCSC MANE track (slow)")
    g3.add_argument("--seq-batch", type=int, default=50, help="regions per Ensembl POST sequence request")
    g3.add_argument("--ispcr", help="path to UCSC isPcr executable (genome-wide uniqueness check)")
    g3.add_argument("--genome-2bit", help="genome .2bit for isPcr (hg38.2bit / hg19.2bit)")
    g3.add_argument("--ispcr-max-size", type=int, default=d.ispcr_max_size)
    a = ap.parse_args(argv)
    if not a.genes and not a.regions:
        ap.error("give a gene list, --regions, or both")

    try:
        p = Params(target=a.target, min_overlap=a.min_overlap, flank=a.flank,
                   product_min=a.product[0], product_opt=a.product[1], product_max=a.product[2],
                   tm_min=a.tm[0], tm_opt=a.tm[1], tm_max=a.tm[2],
                   size_min=a.size[0], size_opt=a.size[1], size_max=a.size[2], gc_min=a.gc[0], gc_max=a.gc[1],
                   max_end_gc=a.max_end_gc, max_poly_g=a.max_poly_g, max_poly_x=a.max_poly_x,
                   max_tm_diff=a.max_tm_diff, penalty_max=a.penalty_max, tm_formula=a.tm_formula,
                   num_return=a.alternatives, relax=not a.no_relax, snp_maf=a.snp_maf,
                   snp_3prime=a.snp_3prime, fallback=not a.strict, ispcr_max_size=a.ispcr_max_size,
                   **salt_params(a.salt, a.tm_formula))
        regions = read_region_list(a.regions) if a.regions else []
    except ValueError as ex:
        ap.error(str(ex))
    if bool(a.ispcr) != bool(a.genome_2bit):
        ap.error("--ispcr and --genome-2bit go together")

    def log(msg):
        print(msg, file=sys.stderr, flush=True)

    ens = ensembl or Ensembl(a.assembly, a.cache or None)
    if ucsc is None and not a.no_ucsc:
        ucsc = UCSC(a.assembly, a.cache or None, a.ucsc_delay)
    pad = p.flank_for(True) + 20
    genes = read_gene_list(a.genes) if a.genes else []
    results, errors = [], []

    # ---- phase 1: structures (MANE track first, Ensembl for the rest) ----
    models = []
    if genes:
        mane = ManeTable(None if (a.no_ucsc or a.no_mane_track) else ucsc, log)
        from_mane = 0
        for i, (sym, tx) in enumerate(genes, 1):
            try:
                m = None if tx else mane.model(sym, pad)
                if m is not None:
                    from_mane += 1
                else:
                    m = build_model(ens, sym, tx, pad, with_seq=False)
                models.append((sym, m))
            except Exception as ex:
                errors.append({"gene": sym, "error": f"{type(ex).__name__}: {ex}"})
                log(f"[{i}/{len(genes)}] {sym:<10} ERROR {ex}")
        log(f"Structures: {len(models)}/{len(genes)} ({from_mane} from the MANE track, "
            f"{len(models) - from_mane} from Ensembl)")
    taken = {sym.upper() for sym, _ in genes}
    for name, chrom, s, e in regions:
        if name.upper() in taken:
            errors.append({"gene": name, "error": "ValueError: region name clashes with a gene symbol"})
            continue
        models.append((name, region_model(name, chrom, s, e, pad)))

    # ---- phase 2: sequences, in batches ---------------------------------
    log(f"Fetching soft-masked sequence for {len(models)} genes/regions ...")
    fetch_sequences(ens, [m for _, m in models], a.seq_batch, log)
    for sym, m in models:
        if not m.seq:
            errors.append({"gene": sym, "error": "LookupError: no sequence returned by Ensembl"})
            log(f"{sym:<10} ERROR no sequence returned")
    models = [(sym, m) for sym, m in models if m.seq]

    # ---- phase 3: annotate and design -----------------------------------
    for i, (sym, m) in enumerate(models, 1):
        try:
            annotate(m, None if a.no_ucsc else ucsc, p, log)
            g = design_gene(m, p, a.assembly)
            g["n_snps"] = len(m.snps)
            results.append(g)
            st = [x["status"] for x in g["amplicons"]]
            log(f"[{i}/{len(models)}] {g['gene']:<10} {g['transcript'] or g['kind']:<20} {len(st)} targets "
                f"({st.count('OK')} OK, {st.count('CHECK')} to check), {len(g['failures'])} failed")
        except Exception as ex:
            errors.append({"gene": sym, "error": f"{type(ex).__name__}: {ex}"})
            log(f"[{i}/{len(models)}] {sym:<10} ERROR {ex}")

    ispcr_ran = False
    if a.ispcr:
        ispcr_ran = run_ispcr(results, a.ispcr, a.genome_2bit, p, log)
        if ispcr_ran:
            for g in results:
                finalize_gene(g, p)

    write_excel(results, errors, p, a.out, a.assembly, ispcr_ran)
    if a.json_dir:
        write_json(results, errors, p, a.json_dir, a.assembly, ispcr_ran)
    amps = [x for g in results for x in g["amplicons"]]
    log(f"Done: {len(results)} genes/regions, {len(amps)} targets with primers "
        f"({sum(x['status'] == 'OK' for x in amps)} OK), "
        f"{sum(len(g['failures']) for g in results)} failed targets, {len(errors)} errors -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
