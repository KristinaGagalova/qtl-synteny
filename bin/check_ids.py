#!/usr/bin/env python3
"""
check_ids.py -- diagnose why miniprot / RBH links are empty.

Nothing here reruns the pipeline: it reads files you already have and
reports whether the identifiers actually join.

The pipeline joins protein evidence to genes by EXACT STRING MATCH at four
points, and every one of them has to agree:

  A  gene table  gene_id       <-> counts table column 1
  B  gene table  transcript_id <-> source protein FASTA header (miniprot Target=)
  C  rbh.m8 column 1           <-> source gene table (gene or transcript id)
  D  rbh.m8 column 2           <-> target gene table (gene or transcript id)

Usage
-----
  check_ids.py \
      --source-genes results/annotation/Norin.source.genes.tsv \
      --target-genes results/annotation/Cadenza.target.genes.tsv \
      --miniprot     results/homology/Norin.miniprot.gff \
      --rbh          results/orthology/Norin.rbh.m8 \
      --source-pep   Triticum_aestivum_norin61...pep.fa \
      --target-pep   Triticum_aestivum_cadenza...pep.fa \
      --source-expr  RawCountsStar-Norin.tsv \
      --target-expr  RawCountsStar-Cadenza.tsv

Every argument is optional; whatever you give it, it checks.
"""

import argparse
import gzip
import re
import sys

ATTR = re.compile(r"([A-Za-z_]+)=([^;]+)")
TYPE_PREFIX = re.compile(
    r"^(?:gene|transcript|mRNA|CDS|protein|exon|rna|ncRNA):", re.IGNORECASE)


def op(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def sample(s, n=3):
    return ", ".join(repr(x) for x in list(s)[:n]) if s else "(none)"


def gene_ids(path):
    genes, tx = set(), set()
    with op(path) as fh:
        fh.readline()
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                genes.add(f[3])
                tx.add(f[4])
    return genes, tx


def pep_ids(path):
    ids = set()
    with op(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def miniprot_targets(path):
    ids = set()
    with op(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            a = dict(ATTR.findall(f[8]))
            t = a.get("Target", "").split()
            if t:
                ids.add(t[0])
    return ids


def m8_cols(path):
    q, t = set(), set()
    with op(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                q.add(f[0]); t.add(f[1])
    return q, t


def expr_ids(path):
    ids = set()
    with op(path) as fh:
        first = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            f = line.split("\t") if "\t" in line else line.split()
            ids.add(f[0])
    return ids


def variants(x):
    out = {x}
    b = TYPE_PREFIX.sub("", x, count=1)
    out.add(b)
    out.add(b.rsplit(".", 1)[0])
    out.add(b + ".1")
    return out


def report(name, a, a_lbl, b, b_lbl):
    if not a or not b:
        print(f"  {name}: SKIPPED (missing input)")
        return
    direct = len(a & b)
    loose = sum(1 for x in a if variants(x) & b)
    verdict = "OK" if direct else ("PREFIX/SUFFIX MISMATCH" if loose else "NO OVERLAP")
    print(f"  {name}: {verdict}")
    print(f"      {a_lbl}: {len(a):>7} ids   e.g. {sample(a)}")
    print(f"      {b_lbl}: {len(b):>7} ids   e.g. {sample(b)}")
    print(f"      exact matches: {direct}    after normalising: {loose}")
    if not direct and loose:
        ex = next((x for x in a if variants(x) & b), None)
        if ex:
            hit = next(iter(variants(ex) & b))
            print(f"      -> {ex!r} only matches as {hit!r}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("source-genes", "target-genes", "miniprot", "rbh",
              "source-pep", "target-pep", "source-expr", "target-expr"):
        ap.add_argument(f"--{a}")
    args = ap.parse_args()

    sg = tg = sgt = tgt_t = set()
    if args.source_genes:
        sg, sgt = gene_ids(args.source_genes)
        print(f"source gene table : {len(sg)} genes, {len(sgt)} transcripts")
        print(f"   gene_id   e.g. {sample(sg)}")
        print(f"   transcript e.g. {sample(sgt)}")
    if args.target_genes:
        tg, tgt_t = gene_ids(args.target_genes)
        print(f"target gene table : {len(tg)} genes, {len(tgt_t)} transcripts")
        print(f"   gene_id   e.g. {sample(tg)}")
        print(f"   transcript e.g. {sample(tgt_t)}")
    print()

    if args.source_pep:
        sp = pep_ids(args.source_pep)
        print("B. source proteins -> source gene table")
        report("   protein FASTA vs transcript_id", sp, "pep", sgt | sg, "gene table")
        print()

    if args.miniprot:
        mt = miniprot_targets(args.miniprot)
        print("B'. miniprot Target= -> source gene table")
        report("   miniprot vs gene table", mt, "miniprot", sgt | sg, "gene table")
        print()

    if args.rbh:
        q, t = m8_cols(args.rbh)
        print("C. rbh.m8 query -> source gene table")
        report("   rbh col1 vs source", q, "rbh q", sgt | sg, "source genes")
        print()
        print("D. rbh.m8 target -> target gene table")
        report("   rbh col2 vs target", t, "rbh t", tgt_t | tg, "target genes")
        print()

    if args.source_expr:
        e = expr_ids(args.source_expr)
        print("A. source counts -> source gene table")
        report("   counts col1 vs gene_id", e, "counts", sg, "gene_id")
        print()
    if args.target_expr:
        e = expr_ids(args.target_expr)
        print("A. target counts -> target gene table")
        report("   counts col1 vs gene_id", e, "counts", tg, "gene_id")
        print()


if __name__ == "__main__":
    main()
