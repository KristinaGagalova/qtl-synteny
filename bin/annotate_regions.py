#!/usr/bin/env python3
"""
annotate_regions.py -- STEPS 3/4: layer protein evidence onto the syntenic
intervals that DNA alignment already defined.

Inputs are the DNA-ranked regions for one QTL plus two independent,
genome-wide protein evidence sets. Both are FILTERED to the regions - the
question is not "where does this protein go in the whole genome" but "does
it land inside the syntenic interval we already identified".

  step 3  miniprot : source proteins spliced onto the target GENOME. A hit
                     is kept when it falls inside a candidate region. This
                     works even where the target annotation is missing or
                     wrong, because it does not use the target GFF at all.
  step 4  RBH      : reciprocal best hits between the two ANNOTATIONS
                     (MMseqs2 easy-rbh). Kept when the target gene lies in a
                     candidate region. Stronger evidence, but requires both
                     annotations to be good.

Outputs
-------
--out-links    one row per protein-level link, with `track` = miniprot | rbh,
               so the viewer can draw them as separate tracks or overlaid.
--out-genes    every gene to display: all source genes inside the QTL, plus
               all target genes inside the candidate regions. This is the set
               the expression panel uses, and it is deliberately NOT limited
               to genes that have a link - a target gene with no orthologue
               call is still worth seeing expression for.
--out-regions  the input regions with protein-evidence counts attached.
"""

import argparse
import re
import sys
from collections import defaultdict

ATTR = re.compile(r"([A-Za-z_]+)=([^;]+)")


def read_genes(path):
    """-> (by_id, by_seq) where by_id maps gene AND transcript id."""
    by_id, by_seq = {}, defaultdict(list)
    with open(path) as fh:
        fh.readline()
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            seqid, s, e, gid, tid, strand = f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]
            rec = (seqid, s, e, strand, gid)
            by_id[gid] = rec
            if tid and tid != gid:
                by_id[tid] = rec
            by_seq[seqid].append((s, e, gid, strand, tid))
    for k in by_seq:
        by_seq[k].sort()
    return by_id, by_seq


def read_regions(path, qtl_id):
    regs = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[i["qtl_id"]] != qtl_id:
                continue
            regs.append(dict(
                rank=int(f[i["rank"]]), tgt_seqid=f[i["tgt_seqid"]],
                tgt_start=int(f[i["tgt_start"]]), tgt_end=int(f[i["tgt_end"]]),
                aligned_bp=int(f[i["aligned_bp"]]),
                n_aln_blocks=int(f[i["n_aln_blocks"]]),
                pct_id=float(f[i["pct_id"]]), strand=f[i["strand"]],
                raw=f, hdr=hdr))
    regs.sort(key=lambda r: r["rank"])
    return regs


def in_regions(regs, seqid, start, end):
    for r in regs:
        if r["tgt_seqid"] == seqid and start < r["tgt_end"] and end > r["tgt_start"]:
            return r
    return None


def parse_miniprot(path, keep_ids=None):
    """mRNA lines of miniprot --gff -> list of hits."""
    hits = []
    if not path:
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            a = dict(ATTR.findall(f[8]))
            tgt = a.get("Target", "").split()
            if not tgt:
                continue
            pid = tgt[0]
            if keep_ids is not None and pid not in keep_ids:
                continue
            hits.append(dict(protein=pid, seqid=f[0], start=int(f[3]) - 1,
                             end=int(f[4]), strand=f[6],
                             ident=float(a.get("Identity", 0)),
                             rank=int(a.get("Rank", 1))))
    return hits


def read_m8(path):
    """query -> list of (target, pident, bits) best first."""
    hits = defaultdict(list)
    if not path:
        return hits
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                hits[f[0]].append((f[1], float(f[2]), float(f[11])))
            except ValueError:
                continue
    for q in hits:
        hits[q].sort(key=lambda x: -x[2])
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--qtl-id", required=True)
    ap.add_argument("--qtl-chrom", required=True)
    ap.add_argument("--qtl-start", type=int, required=True)
    ap.add_argument("--qtl-end", type=int, required=True)
    ap.add_argument("--source-genes", required=True)
    ap.add_argument("--target-genes", required=True)
    ap.add_argument("--miniprot")
    ap.add_argument("--rbh")
    ap.add_argument("--out-links", required=True)
    ap.add_argument("--out-genes", required=True)
    ap.add_argument("--out-regions", required=True)
    ap.add_argument("--miniprot-min-ident", type=float, default=0.5)
    args = ap.parse_args()

    src_by_id, src_by_seq = read_genes(args.source_genes)
    tgt_by_id, tgt_by_seq = read_genes(args.target_genes)
    regs = read_regions(args.regions, args.qtl_id)

    # ---- genes inside the QTL on the source side
    src_genes = [(s, e, gid, strand, tid)
                 for (s, e, gid, strand, tid) in src_by_seq.get(args.qtl_chrom, [])
                 if s < args.qtl_end and e > args.qtl_start]
    src_ids = {gid for (_s, _e, gid, _st, _t) in src_genes}
    src_tids = {tid for (_s, _e, _g, _st, tid) in src_genes if tid}
    qtl_prot_ids = src_ids | src_tids

    # ---- genes inside the candidate regions on the target side
    tgt_genes = []
    for r in regs:
        for (s, e, gid, strand, tid) in tgt_by_seq.get(r["tgt_seqid"], []):
            if s < r["tgt_end"] and e > r["tgt_start"]:
                tgt_genes.append((r["tgt_seqid"], s, e, gid, strand, tid, r["rank"]))
    seen = set()
    tgt_genes = [g for g in tgt_genes if not (g[3] in seen or seen.add(g[3]))]

    links = []
    per_region = defaultdict(lambda: {"miniprot": 0, "rbh": 0})

    # ---- STEP 3: miniprot placements of the QTL's source proteins
    mp = parse_miniprot(args.miniprot, keep_ids=qtl_prot_ids)
    n_mp_total = len(mp)
    for h in mp:
        if h["ident"] < args.miniprot_min_ident:
            continue
        r = in_regions(regs, h["seqid"], h["start"], h["end"])
        if r is None:
            continue
        srec = src_by_id.get(h["protein"])
        if srec is None:
            continue
        sseq, ss, se, sstrand, sgid = srec
        # is there an annotated target gene at this placement?
        overlapping = [g for g in tgt_by_seq.get(h["seqid"], [])
                       if g[0] < h["end"] and g[1] > h["start"]]
        tgene = overlapping[0][2] if overlapping else ""
        links.append(dict(track="miniprot", src_gene=sgid, src_start=ss, src_end=se,
                          src_strand=sstrand, tgt_gene=tgene, tgt_seqid=h["seqid"],
                          tgt_start=h["start"], tgt_end=h["end"],
                          tgt_strand=h["strand"], pident=100.0 * h["ident"],
                          bits=0, rank=r["rank"]))
        per_region[r["tgt_seqid"]]["miniprot"] += 1

    # ---- STEP 4: reciprocal best hits between the two annotations
    rbh = read_m8(args.rbh)
    n_rbh_total = sum(len(v) for v in rbh.values())
    for q, lst in rbh.items():
        srec = src_by_id.get(q)
        if srec is None:
            continue
        sseq, ss, se, sstrand, sgid = srec
        if sgid not in src_ids:
            continue          # only genes inside this QTL
        for (t, pid, bits) in lst:
            trec = tgt_by_id.get(t)
            if trec is None:
                continue
            tseq, ts, te, tstrand, tgid = trec
            r = in_regions(regs, tseq, ts, te)
            if r is None:
                continue      # RBH exists but outside the syntenic interval
            links.append(dict(track="rbh", src_gene=sgid, src_start=ss, src_end=se,
                              src_strand=sstrand, tgt_gene=tgid, tgt_seqid=tseq,
                              tgt_start=ts, tgt_end=te, tgt_strand=tstrand,
                              pident=pid, bits=int(bits), rank=r["rank"]))
            per_region[tseq]["rbh"] += 1

    with open(args.out_links, "w") as out:
        out.write("qtl_id\ttrack\tsrc_gene\tsrc_seqid\tsrc_start\tsrc_end\tsrc_strand\t"
                  "tgt_gene\ttgt_seqid\ttgt_start\ttgt_end\ttgt_strand\t"
                  "pident\tbits\tregion_rank\n")
        for l in links:
            out.write(f"{args.qtl_id}\t{l['track']}\t{l['src_gene']}\t{args.qtl_chrom}\t"
                      f"{l['src_start']}\t{l['src_end']}\t{l['src_strand']}\t"
                      f"{l['tgt_gene']}\t{l['tgt_seqid']}\t{l['tgt_start']}\t"
                      f"{l['tgt_end']}\t{l['tgt_strand']}\t{l['pident']:.1f}\t"
                      f"{l['bits']}\t{l['rank']}\n")

    with open(args.out_genes, "w") as out:
        out.write("qtl_id\tside\tgene_id\ttranscript_id\tseqid\tstart\tend\tstrand\tregion_rank\n")
        for (s, e, gid, strand, tid) in src_genes:
            out.write(f"{args.qtl_id}\tsource\t{gid}\t{tid}\t{args.qtl_chrom}\t"
                      f"{s}\t{e}\t{strand}\t0\n")
        for (seqid, s, e, gid, strand, tid, rank) in tgt_genes:
            out.write(f"{args.qtl_id}\ttarget\t{gid}\t{tid}\t{seqid}\t{s}\t{e}\t"
                      f"{strand}\t{rank}\n")

    with open(args.out_regions, "w") as out:
        out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\trank\ttgt_seqid\t"
                  "aligned_bp\tn_aln_blocks\tpct_id\tstrand\ttgt_start\ttgt_end\t"
                  "tgt_span\tn_tgt_genes\tn_miniprot\tn_rbh\n")
        for r in regs:
            ng = sum(1 for g in tgt_genes if g[0] == r["tgt_seqid"])
            c = per_region[r["tgt_seqid"]]
            out.write(f"{args.qtl_id}\t{args.qtl_chrom}\t{args.qtl_start}\t{args.qtl_end}\t"
                      f"{r['rank']}\t{r['tgt_seqid']}\t{r['aligned_bp']}\t"
                      f"{r['n_aln_blocks']}\t{r['pct_id']:.1f}\t{r['strand']}\t"
                      f"{r['tgt_start']}\t{r['tgt_end']}\t{r['tgt_end']-r['tgt_start']}\t"
                      f"{ng}\t{c['miniprot']}\t{c['rbh']}\n")

    n_mp = sum(1 for l in links if l["track"] == "miniprot")
    n_rb = sum(1 for l in links if l["track"] == "rbh")
    sys.stderr.write(
        f"[annotate_regions] {args.qtl_id}: {len(regs)} regions, "
        f"{len(src_genes)} source genes in QTL, {len(tgt_genes)} target genes in regions, "
        f"{n_mp}/{n_mp_total} miniprot hits in-region, {n_rb}/{n_rbh_total} RBH in-region\n")


if __name__ == "__main__":
    main()
