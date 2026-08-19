#!/usr/bin/env python3
"""
test_coverage_logic.py -- regression tests for the two coordinate-logic
issues that the integration test cannot catch.

The integration test asserts that files exist and contain data. It does NOT
assert genomic coordinates or coverage percentages, so it would happily pass
while the numbers were wrong. These tests assert exact expected values.

Run:  python3 test/test_coverage_logic.py
"""

import os
import subprocess
import sys
import tempfile

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
RANK = os.path.join(BIN, "rank_dna_regions.py")

FAILURES = []


def check(name, got, want):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: got {got!r}, want {want!r}")
    if not ok:
        FAILURES.append(name)


def paf_line(qname, qlen, qs, qe, tname, tlen, ts, te, strand="+", pid=95):
    aln = qe - qs
    return "\t".join(map(str, [
        qname, qlen, qs, qe, strand, tname, tlen, ts, te,
        int(aln * pid / 100), aln, 60]))


def run_rank(paf_text, qtl_start, qtl_end, extra=None):
    d = tempfile.mkdtemp()
    paf = os.path.join(d, "t.paf")
    out = os.path.join(d, "r.tsv")
    open(paf, "w").write(paf_text + "\n")
    cmd = [sys.executable, RANK, "--paf", paf, "--qtl-id", "Q",
           "--src-chrom", "chr1", "--src-start", str(qtl_start),
           "--src-end", str(qtl_end), "--out", out, "--min-aligned-bp", "100"]
    if extra:
        cmd += extra
    subprocess.run(cmd, check=True, capture_output=True)
    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    hdr = {c: i for i, c in enumerate(rows[0])}
    return hdr, rows[1:]


def test_qtl_vs_flank_coverage():
    """The reviewer's scenario, verbatim.

    source chr1, QTL = 2.0-3.0 Mb, source_flank 100 kb
    -> query slice = 1.9-3.1 Mb, header 'chr1:1900001-3100000'

      scafA maps source 2.0-2.8 Mb  => QTL coverage MUST be 80%
      scafB maps source 1.9-2.0 Mb  => flank only, QTL coverage MUST be 0%

    Before the fix, coverage was measured against the whole 1.2 Mb slice,
    so scafB scored non-zero and could influence ranking.
    """
    print("test_qtl_vs_flank_coverage")
    q = "chr1:1900001-3100000"
    # local coords: slice starts at real 1,900,000 (0-based)
    paf = "\n".join([
        paf_line(q, 1200000, 100000, 900000, "scafA", 999999, 1000, 801000),
        paf_line(q, 1200000, 20000, 100000, "scafB", 999999, 1000, 81000),
    ])
    hdr, rows = run_rank(paf, 2000000, 3000000)
    by_seq = {r[hdr["tgt_seqid"]]: r for r in rows}

    check("scafA qtl_cov_pct", by_seq["scafA"][hdr["qtl_cov_pct"]], "80.00")
    check("scafA qtl_cov_bp", by_seq["scafA"][hdr["qtl_cov_bp"]], "800000")
    check("scafB qtl_cov_pct (flank only)", by_seq["scafB"][hdr["qtl_cov_pct"]], "0.00")
    check("scafB qtl_cov_bp (flank only)", by_seq["scafB"][hdr["qtl_cov_bp"]], "0")
    check("scafB left_flank_cov_pct", by_seq["scafB"][hdr["left_flank_cov_pct"]], "80.00")
    check("scafA ranks first", by_seq["scafA"][hdr["rank"]], "1")
    check("scafB ranks second", by_seq["scafB"][hdr["rank"]], "2")


def test_flank_alignment_cannot_outrank_qtl_alignment():
    """A scaffold with MORE total alignment, all of it in the flank, must
    still rank below one with less alignment that actually touches the QTL."""
    print("test_flank_alignment_cannot_outrank_qtl_alignment")
    q = "chr1:1900001-3100000"
    paf = "\n".join([
        # 50 kb, inside the QTL
        paf_line(q, 1200000, 100000, 150000, "smallQTL", 999999, 1000, 51000),
        # 100 kb, entirely flank
        paf_line(q, 1200000, 0, 100000, "bigFlank", 999999, 1000, 101000),
    ])
    hdr, rows = run_rank(paf, 2000000, 3000000)
    by_seq = {r[hdr["tgt_seqid"]]: r for r in rows}
    check("smallQTL ranks first despite less alignment",
          by_seq["smallQTL"][hdr["rank"]], "1")
    check("bigFlank aligned_bp is larger", by_seq["bigFlank"][hdr["aligned_bp"]], "100000")
    check("bigFlank qtl_cov_bp is zero", by_seq["bigFlank"][hdr["qtl_cov_bp"]], "0")


def test_target_side_clustering():
    """One target sequence with a genuine locus and a distant repeat hit must
    yield TWO regions, not one giant span covering everything between."""
    print("test_target_side_clustering")
    q = "chr1:1900001-3100000"
    paf = "\n".join([
        paf_line(q, 1200000, 100000, 300000, "chr2", 500000000, 20000000, 22000000),
        paf_line(q, 1200000, 50000, 70000, "chr2", 500000000, 480000000, 480100000),
    ])

    # with clustering (default gap)
    hdr, rows = run_rank(paf, 2000000, 3000000)
    check("clustering yields 2 regions", len(rows), 2)
    spans = sorted(int(r[hdr["tgt_span"]]) for r in rows)
    check("no region spans the 460 Mb gap", max(spans) < 100_000_000, True)

    # without clustering, for contrast: one giant region
    hdr2, rows2 = run_rank(paf, 2000000, 3000000,
                           extra=["--cluster-max-gap", "999999999"])
    check("disabling clustering gives 1 region", len(rows2), 1)
    check("that region is absurdly large (>400 Mb)",
          int(rows2[0][hdr2["tgt_span"]]) > 400_000_000, True)


def test_clustering_is_noop_on_fragmented_target():
    """On scaffolds shorter than the gap threshold, clustering must change
    nothing - one scaffold, one region."""
    print("test_clustering_is_noop_on_fragmented_target")
    q = "chr1:1900001-3100000"
    paf = "\n".join([
        paf_line(q, 1200000, 100000, 200000, "scaf1", 57000, 1000, 21000),
        paf_line(q, 1200000, 200000, 300000, "scaf1", 57000, 25000, 45000),
    ])
    hdr, rows = run_rank(paf, 2000000, 3000000)
    check("one scaffold -> one region", len(rows), 1)
    check("both blocks retained", rows[0][hdr["n_aln_blocks"]], "2")


if __name__ == "__main__":
    test_qtl_vs_flank_coverage()
    test_flank_alignment_cannot_outrank_qtl_alignment()
    test_target_side_clustering()
    test_clustering_is_noop_on_fragmented_target()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {', '.join(FAILURES)}")
        sys.exit(1)
    print("all coverage-logic tests passed")
