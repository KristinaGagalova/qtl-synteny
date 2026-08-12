#!/usr/bin/env python3
"""
rank_dna_regions.py -- STEP 1/2: define the candidate syntenic intervals for
one QTL from DNA alignment alone.

The source QTL interval is aligned against the whole target assembly; this
ranks the target sequences it hit by total aligned bp and keeps the top N.
Deliberately independent of any protein/ortholog evidence, so regions are
still found when RBH comes back empty (which it does on a fragmented target
with strong homoeologous competition).
"""
import argparse
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True, help="source QTL slice vs whole target")
    ap.add_argument("--qtl-id", required=True)
    ap.add_argument("--src-chrom", required=True)
    ap.add_argument("--src-start", type=int, required=True)
    ap.add_argument("--src-end", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=5)
    ap.add_argument("--min-aligned-bp", type=int, default=2000)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--flank", type=int, default=20000)
    args = ap.parse_args()

    by_t = defaultdict(list)
    with open(args.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                if int(f[11]) < args.min_mapq:
                    continue
                by_t[f[5]].append((int(f[7]), int(f[8]), int(f[10]), int(f[9]), f[4]))
            except ValueError:
                continue

    rows = []
    for seqid, blocks in by_t.items():
        aligned_bp = sum(b[2] for b in blocks)
        if aligned_bp < args.min_aligned_bp:
            continue
        matched = sum(b[3] for b in blocks)
        strands = defaultdict(int)
        for b in blocks:
            strands[b[4]] += b[2]
        rows.append(dict(
            tgt_seqid=seqid, aligned_bp=aligned_bp, n_aln_blocks=len(blocks),
            pct_id=100.0 * matched / aligned_bp if aligned_bp else 0.0,
            tgt_start=max(min(b[0] for b in blocks) - args.flank, 0),
            tgt_end=max(b[1] for b in blocks) + args.flank,
            strand=max(strands, key=strands.get)))

    rows.sort(key=lambda r: -r["aligned_bp"])
    if args.top_n > 0:
        rows = rows[:args.top_n]

    with open(args.out, "w") as out:
        out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\trank\ttgt_seqid\t"
                  "aligned_bp\tn_aln_blocks\tpct_id\tstrand\t"
                  "tgt_start\ttgt_end\ttgt_span\n")
        for rank, r in enumerate(rows, start=1):
            out.write(f"{args.qtl_id}\t{args.src_chrom}\t{args.src_start}\t"
                      f"{args.src_end}\t{rank}\t{r['tgt_seqid']}\t{r['aligned_bp']}\t"
                      f"{r['n_aln_blocks']}\t{r['pct_id']:.1f}\t{r['strand']}\t"
                      f"{r['tgt_start']}\t{r['tgt_end']}\t{r['tgt_end']-r['tgt_start']}\n")

    sys.stderr.write(f"[rank_dna_regions] {args.qtl_id}: {len(by_t)} target seqs hit, "
                     f"{len(rows)} kept\n")


if __name__ == "__main__":
    main()
