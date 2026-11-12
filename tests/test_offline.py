"""Offline tests: synthetic genome, fake Ensembl / UCSC responses, fake isPcr.

Checks
  * primers match the genome at the reported coordinates on both strands;
  * every pair meets the Primer Express rules (Tm 58-60, GC 30-80 %, <= 2 G/C
    in the last five bases, no run of 4, amplicon 50-150 bp);
  * the amplicon covers the target: inside a long exon, across a short one;
  * common SNPs are avoided, or tolerated only outside the 3'-terminal bases
    and reported; the MAF cut-off and --strict behave;
  * repeats are avoided, or tolerated only in the 5' part and reported;
  * segmental duplications are flagged;
  * an in-silico PCR off-target product demotes a pair;
  * custom regions parse and design like a one-exon gene.

Run:  python tests/test_offline.py   (or python -m pytest tests/)
"""
import json
import random
import stat
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import qpcrprimer_batch as qpb  # noqa: E402

rc = qpb.rc
random.seed(7)
CHROM = "".join(random.choice("ACGT") for _ in range(60000))
CHROM = CHROM[:9000] + CHROM[9000:9400].lower() + CHROM[9400:]  # a soft-masked repeat

GENES = {
    # exon 2 is 31 bp: too short for a 50-bp amplicon, so the pair must span it
    "PLUSGENE": dict(strand=1, exons=[(5000, 5300), (8000, 8030), (12000, 13400), (20000, 20600)],
                     tl=(5150, 20250)),
    "MINUSGENE": dict(strand=-1, exons=[(30000, 30500), (33000, 33090), (36000, 37200), (40000, 40250)],
                      tl=(30300, 40100)),
}


class FakeEnsembl:
    def __init__(self, *a, **k):
        pass

    def gene(self, sym):
        if sym not in GENES:
            raise LookupError(f"No valid lookup found for symbol {sym}")
        g = GENES[sym]
        return {"id": f"ENSG_{sym}", "display_name": sym, "Transcript": [{
            "id": f"ENST_{sym}", "version": 1, "strand": g["strand"], "seq_region_name": "1",
            "is_canonical": 1, "MANE": [{"type": "MANE_Select", "refseq_match": f"NM_{sym}.1"}],
            "Exon": [{"id": f"ENSE{i}", "start": s, "end": e} for i, (s, e) in enumerate(g["exons"])],
            "Translation": {"start": g["tl"][0], "end": g["tl"][1], "length": 999}}]}

    def sequence(self, chrom, start, end, strand):
        s = CHROM[start - 1:end]
        return s if strand == 1 else rc(s)

    def sequences(self, regions):
        self.posts = getattr(self, "posts", 0) + 1
        out = {}
        for q in regions:
            chrom, rest = q.split(":", 1)
            span, strand = rest.rsplit(":", 1)
            a, b = (int(x) for x in span.split(".."))
            out[q] = self.sequence(chrom, a, b, int(strand))
        return out


class FakeUCSC:
    def __init__(self, snps=(), repeats=(), segdups=()):
        self._s, self._r, self._d = list(snps), list(repeats), list(segdups)

    def snps(self, chrom, start, end, min_maf):
        return [s for s in self._s if s["end"] >= start and s["start"] <= end and (s["maf"] or 1) >= min_maf]

    def repeats(self, chrom, start, end):
        return self._r

    def segdups(self, chrom, start, end):
        return self._d

    def mane(self, chrom):
        """The UCSC MANE track as the API serves it: 0-based half-open."""
        if chrom != "chr1":
            return []
        rows = []
        for sym, g in GENES.items():
            ex = sorted(g["exons"])
            st = ex[0][0] - 1
            rows.append({
                "chrom": "chr1", "chromStart": st, "chromEnd": ex[-1][1],
                "name": f"ENST_{sym}.1", "geneName": f"ENSG_{sym}", "geneName2": sym,
                "ncbiId": f"NM_{sym}.1", "maneStat": "MANE Select",
                "strand": "+" if g["strand"] == 1 else "-",
                "thickStart": g["tl"][0] - 1, "thickEnd": g["tl"][1],
                "blockSizes": ",".join(str(e - s + 1) for s, e in ex) + ",",
                "chromStarts": ",".join(str(s - 1 - st) for s, _ in ex) + ",",
            })
        return rows


def snp(pos, maf=0.2, rs=None):
    return {"rs": rs or f"rs{pos}", "start": pos, "end": pos, "maf": maf, "class": "snv"}


def quiet(*_):
    pass


def run(sym, p=None, ucsc=None):
    p = p or qpb.Params()
    m = qpb.build_model(FakeEnsembl(), sym, None, p.flank_for(True) + 20)
    qpb.annotate(m, ucsc if ucsc is not None else FakeUCSC(), p, log=quiet)
    return m, qpb.design_gene(m, p, "GRCh38")


def run_region(name, chrom, a, b, p=None, ucsc=None):
    p = p or qpb.Params()
    m = qpb.region_model(name, chrom, a, b, p.flank_for(True) + 20)
    m.seq = FakeEnsembl().sequence(chrom, m.region_start, m.region_end, 1)
    qpb.annotate(m, ucsc if ucsc is not None else FakeUCSC(), p, log=quiet)
    return m, qpb.design_gene(m, p, "GRCh38")


def amp(g, label):
    return next(a for a in g["amplicons"] if a["label"] == label)


def genomic(a, b):
    return CHROM[a - 1:b].upper()


def exon_of(g, a):
    return next(e for e in g["exons"] if e["n"] == a["exons"][0])


def check_pair(m, g, a, b, p):
    seq = genomic(b["start"], b["end"])
    if m.strand == -1:
        seq = rc(seq)
    assert len(seq) == b["product"]
    assert seq.startswith(b["f_seq"]) and seq.endswith(rc(b["r_seq"])), (g["gene"], a["label"])
    assert b["amplicon"].upper() == seq
    pmin, pmax = p.product_range(a["relaxed"])
    assert pmin <= b["product"] <= pmax, (a["label"], b["product"])
    # coverage of the target
    cov = min(b["end"], a["target_end"]) - max(b["start"], a["target_start"]) + 1
    assert cov == b["target_cov"] >= a["min_cov"], (a["label"], cov, a["min_cov"])
    # placement as declared
    e = exon_of(g, a)
    inside = b["start"] >= e["start"] and b["end"] <= e["end"]
    assert b["placement"] == ("exonic" if inside else "spanning")
    fl = p.flank_for(a["relaxed"])
    assert b["start"] >= e["start"] - fl and b["end"] <= e["end"] + fl
    if a["relaxed"]:
        return
    # the Primer Express rules
    for s, tm, gc in ((b["f_seq"], b["f_tm"], b["f_gc"]), (b["r_seq"], b["r_tm"], b["r_gc"])):
        assert p.tm_min <= tm <= p.tm_max, (s, tm)
        assert p.gc_min <= gc <= p.gc_max, (s, gc)
        assert qpb.end_gc(s) <= p.max_end_gc, s
        assert qpb.max_run(s) <= p.max_poly_x, s
        assert p.size_min <= len(s) <= p.size_max
    assert b["d_tm"] <= p.max_tm_diff + 1e-9


def check_gene(sym, p=None, ucsc=None):
    p = p or qpb.Params()
    m, g = run(sym, p, ucsc)
    assert g["refseq"] == f"NM_{sym}.1" and g["selection"] == "MANE Select" and g["kind"] == "gene"
    assert not g["failures"], g["failures"]
    coding = [e["n"] for e in g["exons"] if e["coding"]]
    assert [a["exons"][0] for a in g["amplicons"]] == coding, "one target per coding exon"
    for a in g["amplicons"]:
        assert a["name"] == f'{sym}_qE{a["label"]}'
        for b in a["candidates"]:
            check_pair(m, g, a, b, p)
    return m, g


# ---------------------------------------------------------------- rules
def test_defaults_are_the_primer_express_rules():
    p = qpb.Params()
    assert (p.tm_min, p.tm_opt, p.tm_max) == (58.0, 59.0, 60.0)
    assert (p.gc_min, p.gc_max) == (30.0, 80.0)
    assert (p.product_min, p.product_max) == (50, 150)
    assert p.size_opt == 20 and p.max_end_gc == 2 and p.max_poly_x == 3 and p.max_tm_diff == 2.0
    g = p.primer3_global()
    assert g["PRIMER_MAX_END_GC"] == 2 and g["PRIMER_MAX_POLY_X"] == 3 and g["PRIMER_GC_CLAMP"] == 0
    assert (g["PRIMER_MIN_TM"], g["PRIMER_OPT_TM"], g["PRIMER_MAX_TM"]) == (58.0, 59.0, 60.0)


def test_helpers():
    assert qpb.end_gc("AAAAAAAATTTGC") == 2 and qpb.end_gc("AAAAGCGCG") == 5
    assert qpb.max_run("ACGGGGT") == 4 and qpb.max_run("ACGT") == 1


def test_plus_strand():
    _, g = check_gene("PLUSGENE")
    assert all(a["status"] == "OK" and a["tier"] == "clean" for a in g["amplicons"])
    # long exons: the amplicon stays inside the exon
    for lab in ("1", "3", "4"):
        assert amp(g, lab)["placement"] == "exonic", lab
    # 31-bp exon: the pair spans it and covers all of it
    a = amp(g, "2")
    b = a["pairs"][0]
    assert a["placement"] == "spanning" and a["min_cov"] == 31 and b["target_cov"] == 31
    assert b["start"] < 8000 and b["end"] > 8030
    assert any("flanking" in n for n in a["notes"])


def test_minus_strand():
    _, g = check_gene("MINUSGENE")
    assert g["exons"][0]["start"] == 40000
    assert [a["label"] for a in g["amplicons"]] == ["1", "2", "3", "4"]
    assert all(a["status"] == "OK" for a in g["amplicons"])


def test_full_exon_target():
    _, g = check_gene("PLUSGENE", qpb.Params(target="full"))
    a = amp(g, "1")
    assert (a["target_start"], a["target_end"]) == (5000, 5300)


def test_coding_target_is_covered_even_with_a_long_utr():
    """Exon 4 of PLUSGENE is 601 bp, of which only 251 are coding: the
    amplicon must cover coding bases, not sit in the 3' UTR."""
    _, g = check_gene("PLUSGENE")
    a = amp(g, "4")
    assert (a["target_start"], a["target_end"]) == (20000, 20250)
    for b in a["candidates"]:
        assert b["target_cov"] >= 40 and b["start"] <= 20250 - 39


def test_impossible_settings_rejected():
    for kw in (dict(product_min=200, product_max=150), dict(product_min=30),
               dict(min_overlap=0), dict(min_overlap=80), dict(target="intron"),
               dict(tm_min=61)):
        try:
            qpb.Params(**kw)
        except ValueError:
            continue
        raise AssertionError(kw)


def test_alternatives_are_distinct_pairs():
    _, g = run("PLUSGENE")
    a = amp(g, "3")
    pairs = {(b["f_seq"], b["r_seq"]) for b in a["pairs"]}
    assert len(a["pairs"]) == 3 and len(pairs) == 3
    assert [b["rank"] for b in a["pairs"]] == [1, 2, 3]
    pens = [b["penalty"] for b in a["candidates"]]
    assert pens == sorted(pens), "candidates come in Primer3 penalty order"


# ---------------------------------------------------------------- SNPs
def test_snp_under_best_primer_is_avoided():
    for sym in GENES:
        m0, g0 = run(sym)
        a0 = amp(g0, "3")
        b0 = a0["pairs"][0]
        f_3 = m0.to_genomic(m0.to_local(b0["start"] if m0.strand == 1 else b0["end"]) + b0["f_len"] - 1)
        r_mid = (b0["end"] - b0["r_len"] // 2) if m0.strand == 1 else (b0["start"] + b0["r_len"] // 2)
        uc = FakeUCSC(snps=[snp(f_3), snp(r_mid)])
        _, g = check_gene(sym, ucsc=uc)
        a = amp(g, "3")
        b = a["pairs"][0]
        assert a["tier"] == "clean" and a["status"] == "OK"
        assert not b["f_snps"] and not b["r_snps"]
        assert (b["f_seq"], b["r_seq"]) != (b0["f_seq"], b0["r_seq"])


def dense_snps():
    # a SNP every 12 bp around the 31-bp exon 2: no primer can be SNP-free
    return FakeUCSC(snps=[snp(x, 0.05) for x in range(7880, 8160, 12)])


def test_dense_snps_fall_back_to_5prime_only():
    _, g = check_gene("PLUSGENE", ucsc=dense_snps())
    a = amp(g, "2")
    assert a["tier"] == "snp-5p" and a["status"] == "CHECK"
    for b in a["candidates"]:
        assert b["f_snps"] and all(s["pos3"] > 5 for s in b["f_snps"] + b["r_snps"])
    assert any("SNP rs" in w and "heterozygous deletion" in w for w in a["warnings"])
    assert amp(g, "1")["status"] == "OK"


def test_strict_mode_refuses_snps():
    _, g = run("PLUSGENE", qpb.Params(fallback=False), dense_snps())
    assert [f["label"] for f in g["failures"]] == ["2"]


def test_maf_threshold():
    rare = FakeUCSC(snps=[snp(x, 0.02) for x in range(7880, 8160, 12)])
    _, g = run("PLUSGENE", qpb.Params(snp_maf=0.05), rare)
    assert amp(g, "2")["tier"] == "clean"


def test_snps_between_primers_are_noted_not_flagged():
    m0, g0 = run("PLUSGENE")
    b0 = amp(g0, "3")["pairs"][0]
    mid = (b0["start"] + b0["end"]) // 2
    _, g = run("PLUSGENE", ucsc=FakeUCSC(snps=[snp(mid, rs="rs999")]))
    a = amp(g, "3")
    assert a["status"] == "OK" and a["pairs"][0]["f_seq"] == b0["f_seq"]
    assert a["pairs"][0]["amp_snps"] == [{"rs": "rs999", "maf": 0.2}]
    assert any("rs999" in n for n in a["notes"])


# ---------------------------------------------------------------- repeats
def test_repeat_under_best_primer_is_avoided():
    _, g0 = run("PLUSGENE")
    b0 = amp(g0, "3")["pairs"][0]
    rep = {"start": b0["start"] - 5, "end": b0["start"] + 30, "name": "AluSx", "class": "SINE"}
    _, g = check_gene("PLUSGENE", ucsc=FakeUCSC(repeats=[rep]))
    a = amp(g, "3")
    assert a["tier"] == "clean" and a["pairs"][0]["f_repeat"] == 0 and a["pairs"][0]["r_repeat"] == 0


def test_repeat_only_5prime_allowed_as_fallback():
    # Repeats on both sides of exon 2 leave 40 unique bases around it: less
    # than the shortest amplicon, so a primer must put its 5' part in a repeat.
    reps = [{"start": 7860, "end": 7995, "name": "L1MC", "class": "LINE"},
            {"start": 8036, "end": 8170, "name": "MIRb", "class": "SINE"}]
    _, g = check_gene("PLUSGENE", ucsc=FakeUCSC(repeats=reps))
    a = amp(g, "2")
    assert a["tier"] == "repeat-5p" and a["status"] == "CHECK"
    for b in a["candidates"]:
        f_end = b["start"] + b["f_len"] - 1
        r_3 = b["end"] - b["r_len"] + 1
        assert b["f_repeat"] + b["r_repeat"] > 0 and f_end > 7995 and r_3 < 8036
    assert any("L1MC" in w for w in a["warnings"])


# ---------------------------------------------------------------- homology
def test_segdup_flag():
    sd = {"start": 19000, "end": 21000, "other": "chr9:1000-3000", "identity": 98.7}
    _, g = run("PLUSGENE", ucsc=FakeUCSC(segdups=[sd]))
    a = amp(g, "4")
    assert a["status"] == "CHECK" and any("98.7%" in w for w in a["warnings"])
    assert amp(g, "1")["status"] == "OK"


def fake_ispcr(td, lines):
    mp = Path(td) / "map.json"
    mp.write_text(json.dumps(lines))
    fake = Path(td) / "isPcr"
    fake.write_text(f"""#!{sys.executable}
import json, sys
m = json.load(open({str(mp)!r}))
q, out = sys.argv[2], sys.argv[3]
with open(out, "w") as fh:
    for line in open(q):
        for l in m.get(line.split()[0], []):
            fh.write(l + "\\n")
""")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    return str(fake)


def test_ispcr_off_target_demotes_pair():
    p = qpb.Params()
    _, g = run("PLUSGENE", p)
    a = amp(g, "3")
    first = a["candidates"][0]
    lines = {}
    for ai, x in enumerate(g["amplicons"]):
        for pi, b in enumerate(x["candidates"]):
            pid = f"q0_{ai}_{pi}"
            lines[pid] = [f"chr1\t{b['start'] - 1}\t{b['end']}\t{pid}\t1000\t+"]
            if b is first:
                lines[pid] += [f"chr9\t5000\t5090\t{pid}\t950\t+", f"chr1_KI270706v1_alt\t1\t90\t{pid}\t900\t+"]
    with tempfile.TemporaryDirectory() as td:
        assert qpb.run_ispcr([g], fake_ispcr(td, lines), "genome.2bit", p, log=quiet)
    qpb.finalize_gene(g, p)
    a = amp(g, "3")
    assert (a["pairs"][0]["f_seq"], a["pairs"][0]["r_seq"]) != (first["f_seq"], first["r_seq"])
    assert a["pairs"][0]["ispcr"] == {"checked": True, "on_target": True, "off_target": []}
    demoted = next(b for b in a["candidates"] if b is first)
    assert demoted["ispcr"]["off_target"] == ["chr9:5001-5090 (90 bp)"]  # alt contig ignored
    assert a["status"] == "OK"
    assert qpb.ispcr_text(a["pairs"][0]["ispcr"]) == "specific (1 product)"


def test_ispcr_without_any_specific_pair_is_flagged():
    p = qpb.Params()
    _, g = run("PLUSGENE", p)
    lines = {}
    for ai, x in enumerate(g["amplicons"]):
        for pi, b in enumerate(x["candidates"]):
            pid = f"q0_{ai}_{pi}"
            lines[pid] = [f"chr1\t{b['start'] - 1}\t{b['end']}\t{pid}\t1000\t+"]
            if x["label"] == "1":
                lines[pid].append(f"chr4\t100\t220\t{pid}\t990\t-")
    with tempfile.TemporaryDirectory() as td:
        qpb.run_ispcr([g], fake_ispcr(td, lines), "genome.2bit", p, log=quiet)
    qpb.finalize_gene(g, p)
    a = amp(g, "1")
    assert a["status"] == "CHECK" and any("off-target" in w for w in a["warnings"])
    assert qpb.ispcr_text(a["pairs"][0]["ispcr"]).startswith("2 products")


def test_ispcr_failure_is_reported_not_fatal():
    p = qpb.Params()
    _, g = run("PLUSGENE", p)
    assert qpb.run_ispcr([g], "/nonexistent/isPcr", "x.2bit", p, log=quiet) is False


# ---------------------------------------------------------------- regions
def test_parse_region():
    assert qpb.parse_region("chr7:117,559,590-117,559,800") == ("7", 117559590, 117559800)
    assert qpb.parse_region("X:200-100") == ("X", 100, 200)
    assert qpb.parse_region("chrM:10-20") == ("MT", 10, 20)
    for bad in ("chr7:abc", "7-100", "chr99:1-2", "chr1:1-60000"):
        try:
            qpb.parse_region(bad)
        except ValueError:
            continue
        raise AssertionError(bad)


def test_read_region_list():
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "r.txt"
        f.write_text("# comment\nDEL_1 chr1:5000-5300\nchr1:8000-8030  tiny\nchr1:12000-12400\n\n")
        assert qpb.read_region_list(str(f)) == [("DEL_1", "1", 5000, 5300), ("tiny", "1", 8000, 8030),
                                                 ("chr1_12000_12400", "1", 12000, 12400)]
        f.write_text("A chr1:1-100\na chr1:200-300\n")
        try:
            qpb.read_region_list(str(f))
        except ValueError as ex:
            assert "twice" in str(ex)
        else:
            raise AssertionError("duplicate names must be refused")


def test_region_design():
    p = qpb.Params()
    m, g = run_region("MYREGION", "1", 12100, 12400, p)
    assert g["kind"] == "region" and g["selection"] == "custom region" and not g["failures"]
    a = g["amplicons"][0]
    assert a["name"] == "MYREGION_q" and a["label"] == "" and a["placement"] == "exonic"
    for b in a["candidates"]:
        check_pair(m, g, a, b, p)
        assert 12100 <= b["start"] and b["end"] <= 12400
    # a 31-bp region is spanned, like a short exon
    m, g = run_region("SMALL", "1", 8000, 8030, p)
    a = g["amplicons"][0]
    assert a["placement"] == "spanning" and a["pairs"][0]["target_cov"] == 31


# ---------------------------------------------------------------- outputs
def test_outputs(out=Path("test_output")):
    out.mkdir(exist_ok=True)
    p = qpb.Params()
    uc = FakeUCSC(snps=[snp(x, 0.05) for x in range(7880, 8160, 12)],
                  segdups=[{"start": 36000, "end": 37300, "other": "chr5:10-1310", "identity": 99.1}])
    res = [run(s, p, uc)[1] for s in GENES] + [run_region("REG1", "1", 12100, 12400, p, uc)[1]]
    errs = [{"gene": "NOTAGENE", "error": "LookupError: No valid lookup found for symbol NOTAGENE"}]
    qpb.write_excel(res, errs, p, str(out / "test_primers.xlsx"), "GRCh38", False)
    qpb.write_json(res, errs, p, str(out / "data"), "GRCh38", False)
    from openpyxl import load_workbook
    wb = load_workbook(out / "test_primers.xlsx")
    assert wb.sheetnames == ["Primers", "Alternatives", "Amplicons", "Genes", "Problems", "Parameters"]
    ws = wb["Primers"]
    assert ws.max_row - 1 == sum(len(g["amplicons"]) + len(g["failures"]) for g in res)
    head = [c.value for c in ws[1]]
    assert head == qpb.PRIMER_COLS
    rows = [[c.value for c in r] for r in ws.iter_rows(min_row=2)]
    reg = next(r for r in rows if r[0] == "REG1")
    assert reg[1] == "region" and reg[5] == "REG1_qF" and reg[12] == "REG1_qR"
    ex2 = next(r for r in rows if r[0] == "PLUSGENE" and r[1] == "2")
    assert ex2[2] == "CHECK" and ex2[5] == "PLUSGENE_qE2F"
    j = json.loads((out / "data" / "genes" / "PLUSGENE.json").read_text())
    assert "candidates" not in j["amplicons"][0]
    meta = json.loads((out / "data" / "index.json").read_text())
    assert {x["gene"]: x["kind"] for x in meta["genes"]} == {"PLUSGENE": "gene", "MINUSGENE": "gene", "REG1": "region"}
    plus = next(x for x in meta["genes"] if x["gene"] == "PLUSGENE")
    assert (plus["start"], plus["end"]) == (5000, 20600)


def test_cli_end_to_end():
    with tempfile.TemporaryDirectory() as td:
        gl = Path(td) / "genes.txt"
        gl.write_text("PLUSGENE\nMINUSGENE\nNOTAGENE\n# comment\n")
        rg = Path(td) / "regions.txt"
        rg.write_text("REG1 chr1:12100-12400\nPLUSGENE chr1:1-200\n")
        rc_ = qpb.main([str(gl), "--regions", str(rg), "-o", f"{td}/site/x.xlsx",
                        "--json-dir", f"{td}/site/data", "--cache", ""],
                       ensembl=FakeEnsembl(), ucsc=FakeUCSC())
        assert rc_ == 0 and Path(f"{td}/site/x.xlsx").exists()
        meta = json.loads(Path(f"{td}/site/data/index.json").read_text())
        assert [x["gene"] for x in meta["genes"]] == ["PLUSGENE", "MINUSGENE", "REG1"]
        assert sorted(e["gene"] for e in meta["errors"]) == ["NOTAGENE", "PLUSGENE"]  # name clash refused
        # regions alone are enough
        rc_ = qpb.main(["--regions", str(rg), "-o", f"{td}/r.xlsx", "--cache", ""],
                       ensembl=FakeEnsembl(), ucsc=FakeUCSC())
        assert rc_ == 0


def test_mane_track_matches_the_ensembl_model():
    p = qpb.Params()
    pad = p.flank + 20
    mane = qpb.ManeTable(FakeUCSC(), log=quiet)
    assert set(mane.by_symbol) == {"PLUSGENE", "MINUSGENE"}
    for sym in GENES:
        a = mane.model(sym, pad)
        b = qpb.build_model(FakeEnsembl(), sym, None, pad, with_seq=False)
        assert a.strand == b.strand and a.chrom == b.chrom, sym
        assert a.transcript_version == b.transcript_version and a.refseq == b.refseq, sym
        assert a.region_start == b.region_start and a.region_end == b.region_end, sym
        assert [(e.number, e.start, e.end, e.cds_start, e.cds_end) for e in a.exons] == \
               [(e.number, e.start, e.end, e.cds_start, e.cds_end) for e in b.exons], sym
    assert mane.model("NOTAGENE", pad) is None


def test_sequences_are_fetched_in_one_batch():
    p = qpb.Params()
    ens = FakeEnsembl()
    mane = qpb.ManeTable(FakeUCSC(), log=quiet)
    models = [mane.model(sym, p.flank + 20) for sym in GENES] + [qpb.region_model("R", "1", 100, 400, 50)]
    logged = []
    qpb.fetch_sequences(ens, models, batch=50, log=logged.append)
    assert ens.posts == 1 and not [x for x in logged if "warning" in x], logged
    for m in models:
        assert m.seq == ens.sequence(m.chrom, m.region_start, m.region_end, m.strand)


def test_batch_failure_falls_back_to_single_requests():
    class Broken(FakeEnsembl):
        def sequences(self, regions):
            raise RuntimeError("batch endpoint down")
    mane = qpb.ManeTable(FakeUCSC(), log=quiet)
    models = [mane.model(sym, 120) for sym in GENES]
    qpb.fetch_sequences(Broken(), models, batch=50, log=quiet)
    assert all(m.seq for m in models)


def test_ucsc_tm_matches_what_hgpcr_prints():
    """Values UCSC in-silico PCR printed for a POLG pair (pinned in
    exonprimer-batch too)."""
    assert qpb.tm_ucsc("CCACACCTGACCTGTAGCTG") == 59.3
    assert qpb.tm_ucsc("GCAAAGGCAACAGAGACACG") == 62.0


def test_salt_conditions_reach_primer3():
    p = qpb.Params(salt_div=0.0, dntp=0.0)
    g = p.primer3_global()
    assert g["PRIMER_SALT_DIVALENT"] == 0.0 and g["PRIMER_DNTP_CONC"] == 0.0
    _, gg = run("PLUSGENE", p)
    import primer3
    for b in amp(gg, "3")["pairs"]:
        assert b["f_tm"] == round(primer3.calc_tm(b["f_seq"], mv_conc=50, dv_conc=0, dntp_conc=0,
                                                  dna_conc=50), 1)


def test_noncoding_transcript_is_refused_clearly():
    class NonCoding(FakeEnsembl):
        def gene(self, sym):
            g = super().gene(sym)
            for t in g["Transcript"]:
                t.pop("Translation")
            return g
    try:
        qpb.build_model(NonCoding(), "PLUSGENE", None, 100, with_seq=False)
    except LookupError as ex:
        assert "protein-coding" in str(ex)
    else:
        raise AssertionError("a gene with no coding transcript must be refused")


def test_merge_parts():
    """CI designs genes in chunks and the custom regions in a job of their
    own; merge_parts must stitch them into one site and one workbook."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import merge_parts
    assert merge_parts.SEQ_COLS == qpb.SEQ_COLS, "sequence columns must match the single-run writer"
    p = qpb.Params()
    with tempfile.TemporaryDirectory() as td:
        parts = Path(td) / "parts"
        for name, res in (("000", [run("PLUSGENE", p)[1]]), ("001", [run("MINUSGENE", p)[1]]),
                          ("regions", [run_region("REG1", "1", 12100, 12400, p)[1]])):
            qpb.write_excel(res, [], p, str(parts / name / "part.xlsx"), "GRCh38", True)
            qpb.write_json(res, [], p, str(parts / name / "data"), "GRCh38", True)
        site = Path(td) / "site"
        old = sys.argv
        sys.argv = ["merge_parts.py", str(parts), str(site)]
        try:
            merge_parts.main()
        finally:
            sys.argv = old
        meta = json.loads((site / "data" / "index.json").read_text())
        assert [g["gene"] for g in meta["genes"]] == ["MINUSGENE", "PLUSGENE", "REG1"], "genes first, then regions"
        assert meta["ispcr"] is True and all("start" in g for g in meta["genes"])
        assert sorted(f.name for f in (site / "data" / "genes").iterdir()) == ["MINUSGENE.json", "PLUSGENE.json", "REG1.json"]
        from openpyxl import load_workbook
        wb = load_workbook(site / "qpcr_primers.xlsx", read_only=True)
        assert wb.sheetnames == ["Primers", "Alternatives", "Amplicons", "Genes", "Problems", "Parameters"]
        assert sum(1 for _ in wb["Primers"].iter_rows()) == 1 + 4 + 4 + 1


def test_compare_primer_express(capsys=None):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import compare_primer_express as cpe
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "pe.tsv"
        f.write_text("Name\tSequence\tTm\tPenalty\n"
                     "A_F\tCCACACCTGACCTGTAGCTG\t59,0\t12\n"
                     "A_R\tGCAAAGGCAACAGAGACACG\t59.5\n"
                     "B_F\tGGGGCCCCGGGGCCCCGGGC\n")
        rows = cpe.read_rows(str(f))
        assert [r["name"] for r in rows] == ["A_F", "A_R", "B_F"]
        assert rows[0]["pe_tm"] == 59.0 and rows[0]["pe_penalty"] == 12 and rows[2]["pe_tm"] is None
        c = cpe.check("GGGGCCCCGGGGCCCCGGGC", qpb.Params())
        assert not c["rules"]["3'G/C"] and not c["rules"]["run"] and not c["rules"]["GC"]
        assert cpe.main([str(f)]) == 0


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
