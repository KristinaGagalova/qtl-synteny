#!/usr/bin/env python3
"""
rank_dna_regions.py -- STEP 1/2: define the candidate syntenic intervals for
one QTL from DNA alignment alone, independent of any protein evidence.

Selection is deliberately simple: report EVERY target sequence that aligned
(after basic quality filters), ranked by how much of the QTL it covers,
capped only by --top-n (0 = unlimited). This is an inspection tool for QTLs
on assemblies of any quality, from a handful of chromosomes to hundreds of
thousands of scaffolds -- it should show everything by default and let the
user narrow down, not decide on their behalf.

QTL vs flank, and why the distinction matters
----------------------------------------------
The alignment QUERY is the QTL plus --source_flank padding on each side (see
extract_source_region.py), because the flanks help establish synteny and
orientation even when the QTL's own sequence aligns poorly. But padding the
QUERY does not mean the flanks should count towards "how much of the QTL
does this scaffold explain" -- a scaffold that aligns strongly to 100 kb of
flank and barely touches the actual QTL should not outrank one that covers
the QTL itself. So every alignment block is converted back to real
source-genome coordinates (the query FASTA header encodes its own offset,
`chrom:start-end`, exactly like build_synteny_html.py's _unoffset()) and then
CLIPPED to [--src-start, --src-end) before it contributes to qtl_cov_bp,
qtl_cov_pct, or ranking. Flank-only alignment is reported separately
(left_flank_cov_bp/pct, right_flank_cov_bp/pct) as context, never as
evidence for coverage or rank.

Paralogy grouping: in a polyploid genome the same source region can have
several genuine copies (homoeologs) among the candidates. Any target
sequences whose QTL-clipped covered interval substantially overlaps are
grouped together (region_group) as possible copies of one another. Which
one is the best match is NOT decided here -- that needs protein evidence
(RBH, miniprot), which does not exist yet at this stage. annotate_regions.py
picks the primary within each group once that evidence is available; here
every row just gets its group id and a provisional flag (largest QTL
coverage in the group) that annotate_regions.py may override.
"""

import argparse
import sys
from collections import defaultdict


def merge_len(intervals):
    """Total bp covered by [start, end) intervals, overlap counted once."""
    if not intervals:
        return 0, []
    iv = sorted(intervals)
    merged = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(e - s for s, e in merged), merged


def clip(intervals, lo, hi):
    """Restrict every interval to [lo, hi), dropping empty results."""
    out = []
    for s, e in intervals:
        cs, ce = max(s, lo), min(e, hi)
        if ce > cs:
            out.append((cs, ce))
    return out


def unoffset(name):
    """`seqid:start-end` (1-based, from extract_source_region.py) -> (seqid,
    offset). Same convention as build_synteny_html.py's _unoffset() -- the
    query FASTA header states its own real genome coordinates, so PAF query
    positions can be converted back without needing the flank size passed
    in separately (robust to the slice being clamped at a chromosome start)."""
    if ":" in name and "-" in name.rsplit(":", 1)[-1]:
        base, span = name.rsplit(":", 1)
        try:
            start = int(span.split("-")[0])
            return base, start - 1
        except ValueError:
            pass
    return name, 0


def overlap_frac(a, b):
    """Reciprocal overlap fraction between two interval sets (0..1)."""
    _, ma = merge_len(a)
    _, mb = merge_len(b)
    inter = 0
    i = j = 0
    while i < len(ma) and j < len(mb):
        lo = max(ma[i][0], mb[j][0])
        hi = min(ma[i][1], mb[j][1])
        if hi > lo:
            inter += hi - lo
        if ma[i][1] < mb[j][1]:
            i += 1
        else:
            j += 1
    len_a = sum(e - s for s, e in ma)
    len_b = sum(e - s for s, e in mb)
    if len_a == 0 or len_b == 0:
        return 0.0
    return min(inter / len_a, inter / len_b)


def group_by_overlap(seqid_iv, min_frac):
    """Union-find clustering of scaffolds whose QTL-clipped intervals overlap
    by at least `min_frac` (reciprocal). Returns {seqid: group_id}, 1-based.

    NOTE this is single-linkage (transitive): if A-B overlap and B-C overlap
    but A-C do not, all three still end up in one group. That is a known,
    accepted trade-off (documented in the README) rather than a bug -- it
    keeps grouping simple and cheap. If it causes oversized groups in a
    given dataset, lower --group-overlap-frac.
    """
    ids = list(seqid_iv)
    parent = {s: s for s in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if overlap_frac(seqid_iv[ids[i]], seqid_iv[ids[j]]) >= min_frac:
                union(ids[i], ids[j])

    roots = {}
    groups = {}
    for s in ids:
        r = find(s)
        if r not in roots:
            roots[r] = len(roots) + 1
        groups[s] = roots[r]
    return groups


def cluster_blocks(blocks, max_gap):
    """Split one target sequence's alignment blocks into locally-clustered
    regions, so a "region" means a syntenic locus rather than a whole
    chromosome.

    Blocks are sorted by target position and cut wherever the gap to the
    next block exceeds `max_gap`. Without this, a chromosome-scale target
    with a real hit at 20 Mb and a repeat hit at 480 Mb produces a single
    460 Mb "region": a huge lane full of irrelevant genes, and protein
    evidence 460 Mb away counted as support for the locus. On a fragmented
    target (scaffolds shorter than max_gap) this is a no-op, every scaffold
    yields exactly one cluster.

    Returns a list of block-lists.
    """
    if not blocks:
        return []
    ordered = sorted(blocks, key=lambda b: (b[0], b[1]))
    clusters = [[ordered[0]]]
    for b in ordered[1:]:
        prev_end = max(x[1] for x in clusters[-1])
        if b[0] - prev_end > max_gap:
            clusters.append([b])
        else:
            clusters[-1].append(b)
    return clusters


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True, help="source QTL slice vs whole target")
    ap.add_argument("--qtl-id", required=True)
    ap.add_argument("--src-chrom", required=True)
    ap.add_argument("--src-start", type=int, required=True,
                    help="QTL start, real source-genome coordinates (0-based)")
    ap.add_argument("--src-end", type=int, required=True,
                    help="QTL end, real source-genome coordinates")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=0,
                    help="max regions to keep, ranked by QTL coverage; "
                         "0 = unlimited (default: report everything found)")
    ap.add_argument("--min-aligned-bp", type=int, default=500,
                    help="basic noise filter, not a selection strategy")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--min-pct-id", type=float, default=0.0,
                    help="drop targets whose alignment identity is below "
                         "this before ranking (a homoeolog from another "
                         "subgenome often covers well but at lower identity)")
    ap.add_argument("--group-overlap-frac", type=float, default=0.5,
                    help="scaffolds are grouped as possible genome copies "
                         "when their QTL-clipped covered intervals "
                         "reciprocally overlap by at least this fraction")
    ap.add_argument("--cluster-max-gap", type=int, default=500000,
                    help="split one target sequence into separate regions "
                         "when consecutive alignment blocks are further "
                         "apart than this. Prevents a chromosome-scale "
                         "target collapsing a real locus and a distant "
                         "repeat hit into one giant region. No effect on "
                         "scaffolds shorter than this value.")
    ap.add_argument("--flank", type=int, default=20000,
                    help="bp padding around each displayed TARGET window "
                         "(unrelated to --source_flank, which pads the "
                         "source QUERY upstream of this script)")
    args = ap.parse_args()

    qtl_start, qtl_end = args.src_start, args.src_end
    qtl_len = max(qtl_end - qtl_start, 1)

    by_t = defaultdict(list)
    slice_offset = None
    slice_qlen = 0
    with open(args.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                if int(f[11]) < args.min_mapq:
                    continue
                qname, qoff = unoffset(f[0])
                if slice_offset is None:
                    slice_offset = qoff
                qlen = int(f[1])
                slice_qlen = max(slice_qlen, qlen)
                # local PAF query coords -> real source-genome coordinates
                real_qs = qoff + int(f[2])
                real_qe = qoff + int(f[3])
                by_t[f[5]].append((int(f[7]), int(f[8]), int(f[10]), int(f[9]),
                                   f[4], real_qs, real_qe))
            except ValueError:
                continue

    if slice_offset is None:
        slice_offset = 0
    slice_start = slice_offset
    slice_end = slice_offset + slice_qlen

    rows = []
    qtl_iv_of = {}       # region key -> QTL-clipped source intervals
    left_len = max(qtl_start - slice_start, 1)
    right_len = max(slice_end - qtl_end, 1)

    # The unit of analysis is (target sequence, local alignment cluster),
    # not the whole target sequence - see cluster_blocks().
    for seqid, all_blocks in by_t.items():
        for blocks in cluster_blocks(all_blocks, args.cluster_max_gap):
            aligned_bp = sum(b[2] for b in blocks)
            if aligned_bp < args.min_aligned_bp:
                continue
            matched = sum(b[3] for b in blocks)
            pct_id = 100.0 * matched / aligned_bp if aligned_bp else 0.0
            if pct_id < args.min_pct_id:
                continue
            strands = defaultdict(int)
            for b in blocks:
                strands[b[4]] += b[2]

            full_iv = [(b[5], b[6]) for b in blocks]
            qtl_iv = clip(full_iv, qtl_start, qtl_end)
            left_iv = clip(full_iv, slice_start, qtl_start)
            right_iv = clip(full_iv, qtl_end, slice_end)

            qtl_cov, _ = merge_len(qtl_iv)
            left_cov, _ = merge_len(left_iv)
            right_cov, _ = merge_len(right_iv)
            full_cov, _ = merge_len(full_iv)
            tgt_cov, _ = merge_len([(b[0], b[1]) for b in blocks])

            tgt_start = max(min(b[0] for b in blocks) - args.flank, 0)
            tgt_end = max(b[1] for b in blocks) + args.flank
            # unique key per region, since one seqid may now yield several
            key = (seqid, tgt_start, tgt_end)
            qtl_iv_of[key] = qtl_iv

            rows.append(dict(
                key=key,
                tgt_seqid=seqid, aligned_bp=aligned_bp, n_aln_blocks=len(blocks),
                pct_id=pct_id,
                qtl_cov_bp=qtl_cov, qtl_cov_pct=100.0 * qtl_cov / qtl_len,
                left_flank_cov_bp=left_cov, left_flank_cov_pct=100.0 * left_cov / left_len,
                right_flank_cov_bp=right_cov, right_flank_cov_pct=100.0 * right_cov / right_len,
                full_query_cov_bp=full_cov,
                tgt_cov_bp=tgt_cov,
                tgt_start=tgt_start, tgt_end=tgt_end,
                strand=max(strands, key=strands.get)))

    # Rank by QTL coverage - NOT full-query coverage. A scaffold that only
    # aligns to the flanks must not outrank one that explains the QTL.
    rows.sort(key=lambda r: (-r["qtl_cov_bp"], -r["pct_id"]))

    n_before_cap = len(rows)
    if args.top_n > 0:
        rows = rows[: args.top_n]

    # Paralogy grouping uses QTL-clipped intervals - two scaffolds should
    # not be called "possible genome copies" just because they share flank
    # homology; what matters is whether they explain the same part of the
    # QTL itself.
    # keyed by region (seqid, start, end), not seqid - one scaffold can now
    # contribute several independent regions
    kept_iv = {r["key"]: qtl_iv_of[r["key"]] for r in rows}
    groups = group_by_overlap(kept_iv, args.group_overlap_frac) if kept_iv else {}
    group_size = defaultdict(int)
    for g in groups.values():
        group_size[g] += 1
    best_in_group = {}
    for r in rows:
        g = groups[r["key"]]
        if g not in best_in_group or r["qtl_cov_bp"] > best_in_group[g][1]:
            best_in_group[g] = (r["tgt_seqid"], r["qtl_cov_bp"])

    union_all, _ = merge_len([iv for ivs in qtl_iv_of.values() for iv in ivs])
    union_kept, _ = merge_len([iv for ivs in kept_iv.values() for iv in ivs])

    with open(args.out, "w") as out:
        out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\tqtl_len\t"
                  "slice_start\tslice_end\trank\t"
                  "tgt_seqid\taligned_bp\tn_aln_blocks\tpct_id\tstrand\t"
                  "qtl_cov_bp\tqtl_cov_pct\t"
                  "left_flank_cov_bp\tleft_flank_cov_pct\t"
                  "right_flank_cov_bp\tright_flank_cov_pct\t"
                  "full_query_cov_bp\ttgt_cov_bp\t"
                  "tgt_start\ttgt_end\ttgt_span\t"
                  "region_group\tgroup_size\tprovisional_group_best\t"
                  "union_qtl_cov_bp\tunion_qtl_cov_pct\n")
        for rank, r in enumerate(rows, start=1):
            g = groups[r["key"]]
            prov_best = int(best_in_group[g][0] == r["tgt_seqid"]
                            and best_in_group[g][1] == r["qtl_cov_bp"])
            out.write(f"{args.qtl_id}\t{args.src_chrom}\t{qtl_start}\t{qtl_end}\t"
                      f"{qtl_len}\t{slice_start}\t{slice_end}\t{rank}\t"
                      f"{r['tgt_seqid']}\t{r['aligned_bp']}\t{r['n_aln_blocks']}\t"
                      f"{r['pct_id']:.1f}\t{r['strand']}\t"
                      f"{r['qtl_cov_bp']}\t{r['qtl_cov_pct']:.2f}\t"
                      f"{r['left_flank_cov_bp']}\t{r['left_flank_cov_pct']:.2f}\t"
                      f"{r['right_flank_cov_bp']}\t{r['right_flank_cov_pct']:.2f}\t"
                      f"{r['full_query_cov_bp']}\t{r['tgt_cov_bp']}\t"
                      f"{r['tgt_start']}\t{r['tgt_end']}\t{r['tgt_end']-r['tgt_start']}\t"
                      f"{g}\t{group_size[g]}\t{prov_best}\t"
                      f"{union_kept}\t{100.0*union_kept/qtl_len:.2f}\n")

    n_groups = len(set(groups.values()))
    sys.stderr.write(
        f"[rank_dna_regions] {args.qtl_id}: {len(by_t)} target seqs hit, "
        f"{len(rows)}/{n_before_cap} kept"
        f"{f' (capped by --top-n {args.top_n})' if args.top_n > 0 and n_before_cap > len(rows) else ''}, "
        f"{n_groups} paralogy group(s); QTL covered "
        f"{union_kept:,}/{qtl_len:,} bp ({100.0*union_kept/qtl_len:.1f}%) by kept "
        f"targets, {union_all:,} bp by all hits (flank-only alignment excluded "
        f"from these numbers)\n")


if __name__ == "__main__":
    main()
