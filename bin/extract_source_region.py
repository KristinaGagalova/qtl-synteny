#!/usr/bin/env python3
"""
extract_source_region.py -- pull just the source QTL interval (plus flank)
out of the source assembly, for one QTL.

Decoupled from gene/ortholog evidence on purpose: this is the query for a
DNA-vs-whole-target-index alignment, which must be able to run even when no
protein links exist between the two varieties for this QTL.
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-fasta", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--flank", type=int, default=100000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    s2 = max(args.start - args.flank, 0)
    written = False
    with open(args.source_fasta) as fh, open(args.out, "w") as fo:
        name, buf, in_seq = None, [], False

        def flush():
            nonlocal written
            if name != args.chrom or not buf:
                return
            seq = "".join(buf)
            e2 = min(args.end + args.flank, len(seq))
            if e2 <= s2:
                return
            sub = seq[s2:e2]
            fo.write(f">{name}:{s2+1}-{e2}\n")
            for k in range(0, len(sub), 60):
                fo.write(sub[k:k+60] + "\n")
            written = True

        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                buf = []
                in_seq = (name == args.chrom)
            elif in_seq:
                buf.append(line.strip())
        flush()

    if not written:
        sys.stderr.write(f"[extract_source_region] WARNING: chrom '{args.chrom}' "
                         f"not found in {args.source_fasta}\n")
    else:
        sys.stderr.write(f"[extract_source_region] wrote {args.chrom}:{s2+1}-"
                         f"{min(args.end+args.flank, s2+10**12)}\n")


if __name__ == "__main__":
    main()
