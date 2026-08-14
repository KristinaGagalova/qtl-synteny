#!/usr/bin/env python3
"""
rank_dna_regions.py -- STEP 1/2: define the candidate syntenic intervals for
one QTL from DNA alignment alone, independent of any protein evidence.

Selection is deliberately simple: report EVERY target sequence that aligned
(after basic quality filters), ranked by how much of the source interval it
covers, capped only by --top-n (0 = unlimited). This is an inspection tool
for QTLs on assemblies of any quality, from a handful of chromosomes to
hundreds of thousands of scaffolds -- it should show everything by default
and let the user narrow down, not decide on their behalf.

Paralogy grouping: in a polyploid genome the same source region can have
several genuine copies (homoeologs) among the candidates. Any target
sequences whose covered SOURCE interval substantially overlaps are grouped
together (region_group) as possible copies of one another. Which one is the
best match is NOT decided here -- that needs protein evidence (RBH,
miniprot), which does not exist yet at this stage. annotate_regions.py picks
the primary within each group once that evidence is available; here every
row just gets its group id and a provisional flag (largest DNA coverage in
the group) that annotate_regions.py may override.
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
    """Union-find clustering of scaffolds whose source intervals overlap by
    at least `min_frac` (reciprocal). Returns {seqid: group_id}, 1-based."""
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True, help="source QTL slice vs whole target")
    ap.add_argument("--qtl-id", required=True)
    ap.add_argument("--src-chrom", required=True)
    ap.add_argument("--src-start", type=int, required=True)
    ap.add_argument("--src-end", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=0,
                    help="max regions to keep, ranked by source coverage; "
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
                         "when their covered source intervals reciprocally "
                         "overlap by at least this fraction")
    ap.add_argument("--flank", type=int, default=20000)
    args = ap.parse_args()

    by_t = defaultdict(list)
    qlen_seen = set()
    with open(args.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                if int(f[11]) < args.min_mapq:
                    continue
                by_t[f[5]].append((int(f[7]), int(f[8]), int(f[10]), int(f[9]),
                                   f[4], int(f[2]), int(f[3])))
                qlen_seen.add(int(f[1]))
            except ValueError:
                continue

    src_slice_len = max(qlen_seen) if qlen_seen else (args.src_end - args.src_start)

    rows = []
    src_iv_of = {}
    for seqid, blocks in by_t.items():
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

        src_iv = [(b[5], b[6]) for b in blocks]
        src_cov, _ = merge_len(src_iv)
        tgt_cov, _ = merge_len([(b[0], b[1]) for b in blocks])
        src_iv_of[seqid] = src_iv

        rows.append(dict(
            tgt_seqid=seqid, aligned_bp=aligned_bp, n_aln_blocks=len(blocks),
            pct_id=pct_id, src_cov_bp=src_cov,
            src_cov_pct=100.0 * src_cov / src_slice_len if src_slice_len else 0.0,
            tgt_cov_bp=tgt_cov,
            tgt_start=max(min(b[0] for b in blocks) - args.flank, 0),
            tgt_end=max(b[1] for b in blocks) + args.flank,
            strand=max(strands, key=strands.get)))

    rows.sort(key=lambda r: (-r["src_cov_bp"], -r["pct_id"]))

    n_before_cap = len(rows)
    if args.top_n > 0:
        rows = rows[: args.top_n]

    kept_iv = {r["tgt_seqid"]: src_iv_of[r["tgt_seqid"]] for r in rows}
    groups = group_by_overlap(kept_iv, args.group_overlap_frac) if kept_iv else {}
    group_size = defaultdict(int)
    for g in groups.values():
        group_size[g] += 1
    best_in_group = {}
    for r in rows:
        g = groups[r["tgt_seqid"]]
        if g not in best_in_group or r["src_cov_bp"] > best_in_group[g][1]:
            best_in_group[g] = (r["tgt_seqid"], r["src_cov_bp"])

    union_all, _ = merge_len([iv for ivs in src_iv_of.values() for iv in ivs])
    union_kept, _ = merge_len([iv for ivs in kept_iv.values() for iv in ivs])

    with open(args.out, "w") as out:
        out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\tsrc_slice_len\trank\t"
                  "tgt_seqid\taligned_bp\tn_aln_blocks\tpct_id\tstrand\t"
                  "src_cov_bp\tsrc_cov_pct\ttgt_cov_bp\t"
                  "tgt_start\ttgt_end\ttgt_span\t"
                  "region_group\tgroup_size\tprovisional_group_best\t"
                  "union_src_cov_bp\tunion_src_cov_pct\n")
        for rank, r in enumerate(rows, start=1):
            g = groups[r["tgt_seqid"]]
            prov_best = int(best_in_group[g][0] == r["tgt_seqid"])
            out.write(f"{args.qtl_id}\t{args.src_chrom}\t{args.src_start}\t"
                      f"{args.src_end}\t{src_slice_len}\t{rank}\t{r['tgt_seqid']}\t"
                      f"{r['aligned_bp']}\t{r['n_aln_blocks']}\t{r['pct_id']:.1f}\t"
                      f"{r['strand']}\t{r['src_cov_bp']}\t{r['src_cov_pct']:.2f}\t"
                      f"{r['tgt_cov_bp']}\t{r['tgt_start']}\t{r['tgt_end']}\t"
                      f"{r['tgt_end']-r['tgt_start']}\t{g}\t{group_size[g]}\t"
                      f"{prov_best}\t{union_kept}\t"
                      f"{100.0*union_kept/src_slice_len if src_slice_len else 0:.2f}\n")

    n_groups = len(set(groups.values()))
    sys.stderr.write(
        f"[rank_dna_regions] {args.qtl_id}: {len(by_t)} target seqs hit, "
        f"{len(rows)}/{n_before_cap} kept"
        f"{f' (capped by --top-n {args.top_n})' if args.top_n > 0 and n_before_cap > len(rows) else ''}, "
        f"{n_groups} paralogy group(s); source covered "
        f"{union_kept:,}/{src_slice_len:,} bp "
        f"({100.0*union_kept/src_slice_len if src_slice_len else 0:.1f}%) by kept targets, "
        f"{union_all:,} bp by all hits\n")


if __name__ == "__main__":
    main()
