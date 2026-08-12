#!/usr/bin/env python3
"""
gff_to_gene_bed.py -- flatten a GFF3/GTF into a gene table used by the
synteny viewer.

One representative transcript per gene (longest CDS), so the display is not
cluttered with isoforms that all sit on the same locus.

Output BED-like TSV: seqid  start  end  gene_id  transcript_id  strand
"""

import argparse
import re
import sys
from collections import defaultdict

ATTR = re.compile(r"([^;=\s]+)=([^;]*)")


def attrs(s):
    d = dict(ATTR.findall(s))
    if not d:  # GTF style
        d = dict(re.findall(r'(\S+) "([^"]*)"', s))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--feature", default="mRNA", help="transcript feature type")
    ap.add_argument("--id-attr", default="ID")
    ap.add_argument("--parent-attr", default="Parent")
    ap.add_argument("--strip-prefix", action="store_true",
                    help="drop a leading 'type:' from IDs (Ensembl-style GFF3)")
    args = ap.parse_args()

    def clean(x):
        return x.split(":", 1)[1] if (args.strip_prefix and ":" in x) else x

    tx = {}
    cdslen = defaultdict(int)
    with open(args.gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] == args.feature:
                tid = clean(a.get(args.id_attr, ""))
                if not tid:
                    continue
                tx[tid] = (f[0], int(f[3]) - 1, int(f[4]), f[6],
                           clean(a.get(args.parent_attr, tid)))
            elif f[2] == "CDS":
                for par in a.get(args.parent_attr, "").split(","):
                    par = clean(par)
                    if par:
                        cdslen[par] += int(f[4]) - int(f[3]) + 1

    best = {}
    for tid, (c, s, e, strand, gene) in tx.items():
        L = cdslen.get(tid, e - s)
        if gene not in best or L > best[gene][0]:
            best[gene] = (L, tid, c, s, e, strand)

    with open(args.out, "w") as out:
        out.write("seqid\tstart\tend\tgene_id\ttranscript_id\tstrand\n")
        for gene, (L, tid, c, s, e, strand) in sorted(
                best.items(), key=lambda x: (x[1][2], x[1][3])):
            out.write(f"{c}\t{s}\t{e}\t{gene}\t{tid}\t{strand}\n")

    sys.stderr.write(f"[gff_to_gene_bed] {len(best)} genes from {len(tx)} transcripts\n")


if __name__ == "__main__":
    main()
