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

# Ensembl/EI-style GFF3 prefixes IDs with the feature type:
#   ID=gene:TraesCAD_..._01G000100   Parent=gene:TraesCAD_..._01G000100
#   ID=transcript:TraesCAD_..._01G000100.1
# Expression tables and protein FASTAs use the bare accession, so these
# prefixes must come off or nothing downstream will match. Only these known
# type words are stripped, so a real accession containing ':' is left alone.
TYPE_PREFIX = re.compile(
    r"^(?:gene|transcript|mRNA|CDS|protein|exon|rna|ncRNA):", re.IGNORECASE)


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
        if not x:
            return x
        # always remove a known feature-type prefix (safe, targeted)
        x = TYPE_PREFIX.sub("", x, count=1)
        # --strip-prefix additionally removes ANY leading 'word:' for
        # non-standard annotations
        if args.strip_prefix and ":" in x:
            x = x.split(":", 1)[1]
        return x

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
                # Prefer an explicit gene_id attribute when the annotation
                # provides one - it is the accession the counts table uses.
                gene = a.get("gene_id") or a.get(args.parent_attr, tid)
                tx[tid] = (f[0], int(f[3]) - 1, int(f[4]), f[6], clean(gene))
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
    if best:
        k = next(iter(best))
        _l, _t, _c, _s, _e, _st = best[k]
        sys.stderr.write(f"[gff_to_gene_bed] example gene_id={k!r} "
                         f"transcript_id={_t!r}\n")


if __name__ == "__main__":
    main()
