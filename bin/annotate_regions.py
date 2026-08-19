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
from collections import Counter, defaultdict

ATTR = re.compile(r"([A-Za-z_]+)=([^;]+)")

# Ensembl/EI GFF3 prefixes IDs with the feature type. The gene table has
# these stripped, but protein FASTA headers (and therefore miniprot Target=
# and the MMseqs m8 columns) may or may not - it depends entirely on how the
# proteome file was produced. Rather than assume, index every gene under all
# plausible spellings so the join works either way.
TYPE_PREFIX = re.compile(
    r"^(?:gene|transcript|mRNA|CDS|protein|exon|rna|ncRNA):", re.IGNORECASE)


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


def _unoffset(name):
    """Same convention as rank_dna_regions.py / build_synteny_html.py:
    `seqid:start-end` (1-based) -> (seqid, offset)."""
    if ":" in name and "-" in name.rsplit(":", 1)[-1]:
        base, span = name.rsplit(":", 1)
        try:
            return base, int(span.split("-")[0]) - 1
        except ValueError:
            pass
    return name, 0


def read_paf_src_intervals(path, qtl_start, qtl_end):
    """target_seqid -> list of QTL-CLIPPED (query_start, query_end) blocks,
    in real source-genome coordinates.

    Needed here to compute the two-way source coverage split (scaffolds
    WITH vs WITHOUT protein evidence). This clips to the QTL exactly like
    rank_dna_regions.py does, for the same reason: the alignment query
    includes source_flank padding, and a scaffold that only aligns to the
    flank must not count towards "coverage of the QTL"."""
    by_t = defaultdict(list)
    if not path:
        return by_t
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                qseq, qoff = _unoffset(f[0])
                s, e = qoff + int(f[2]), qoff + int(f[3])
                cs, ce = max(s, qtl_start), min(e, qtl_end)
                if ce > cs:
                    by_t[f[5]].append((cs, ce))
            except ValueError:
                continue
    return by_t


TX_SUFFIX = re.compile(r"\.\d+$")


def norm_id(x):
    """Canonical form of an identifier, so both sides of a join collapse to
    the same string regardless of how the file spells it:

        'transcript:TraesCAD_..._01G000100.1' -> 'TraesCAD_..._01G000100'
        'TraesCAD_..._01G000100.1'            -> 'TraesCAD_..._01G000100'
        'TraesCAD_..._01G000100'              -> 'TraesCAD_..._01G000100'

    Normalising BOTH sides is what makes this work; expanding one side into
    variants does not, because a bare id can never generate a prefixed one.
    """
    if not x:
        return x
    x = TYPE_PREFIX.sub("", x, count=1)
    return TX_SUFFIX.sub("", x)


def read_genes(path):
    """-> (by_id, by_seq).

    by_id is deliberately permissive: every gene is indexed under its gene
    id, its transcript id, and the prefix/suffix variants of both, because
    the proteome file that miniprot and MMseqs saw may spell identifiers
    differently from the GFF (an Ensembl 'transcript:' prefix, a trailing
    '.1'). Exact-match-only joins here silently produced zero links.
    """
    by_id, by_seq = {}, defaultdict(list)
    with open(path) as fh:
        fh.readline()
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            seqid, s, e, gid, tid, strand = f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]
            rec = (seqid, s, e, strand, gid)
            for key in {gid, tid, norm_id(gid), norm_id(tid)}:
                if key:
                    by_id.setdefault(key, rec)
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
            def g(col, default="0"):
                return f[i[col]] if col in i else default
            regs.append(dict(
                rank=int(f[i["rank"]]), tgt_seqid=f[i["tgt_seqid"]],
                tgt_start=int(f[i["tgt_start"]]), tgt_end=int(f[i["tgt_end"]]),
                aligned_bp=int(f[i["aligned_bp"]]),
                n_aln_blocks=int(f[i["n_aln_blocks"]]),
                pct_id=float(f[i["pct_id"]]), strand=f[i["strand"]],
                qtl_len=int(g("qtl_len")),
                qtl_cov_bp=int(g("qtl_cov_bp")),
                qtl_cov_pct=float(g("qtl_cov_pct")),
                left_flank_cov_bp=int(g("left_flank_cov_bp")),
                left_flank_cov_pct=float(g("left_flank_cov_pct")),
                right_flank_cov_bp=int(g("right_flank_cov_bp")),
                right_flank_cov_pct=float(g("right_flank_cov_pct")),
                full_query_cov_bp=int(g("full_query_cov_bp")),
                tgt_cov_bp=int(g("tgt_cov_bp")),
                region_group=int(g("region_group", "0")),
                group_size=int(g("group_size", "1")),
                provisional_group_best=int(g("provisional_group_best", "1")),
                union_qtl_cov_bp=int(g("union_qtl_cov_bp")),
                union_qtl_cov_pct=float(g("union_qtl_cov_pct")),
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
    ap.add_argument("--out-coverage", required=True)
    ap.add_argument("--paf", help="source QTL slice vs whole target (for "
                    "the with/without-gene-evidence coverage split)")
    ap.add_argument("--miniprot-min-ident", type=float, default=0.5)
    ap.add_argument("--homoeolog-min-frac", type=float, default=0.6,
                    help="fraction of a scaffold's RBH partners that must sit "
                         "on the QTL's own chromosome to call it orthologous")
    ap.add_argument("--homoeolog-evidence", choices=["both", "rbh", "miniprot"],
                    default="both",
                    help="which evidence decides the homoeolog call. "
                         "'miniprot' needs no target annotation at all.")
    ap.add_argument("--drop-homoeologs", action="store_true",
                    help="omit regions called HOMOEOLOG_LIKELY entirely "
                         "(default: keep and flag them)")
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
    # membership is tested on the normalised form, so the proteome may spell
    # ids however it likes
    qtl_prot_ids = {norm_id(x) for x in (src_ids | src_tids)}

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
    mp_all = parse_miniprot(args.miniprot)
    mp_all_unfiltered = mp_all
    mp_seen = len(mp_all)
    mp_example = mp_all[0]["protein"] if mp_all else None
    mp = [h for h in mp_all if norm_id(h["protein"]) in qtl_prot_ids]
    n_mp_total = len(mp)
    for h in mp:
        if h["ident"] < args.miniprot_min_ident:
            continue
        r = in_regions(regs, h["seqid"], h["start"], h["end"])
        if r is None:
            continue
        srec = src_by_id.get(h["protein"]) or src_by_id.get(norm_id(h["protein"]))
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
    rbh_example = next(iter(rbh), None)
    rbh_resolved = sum(1 for q in rbh
                       if src_by_id.get(q) or src_by_id.get(norm_id(q)))

    # ---- homoeolog discrimination -------------------------------------
    # Bread wheat is hexaploid: every locus has A, B and D copies at ~92-96%
    # identity, while two cultivars of the SAME subgenome are ~99%+. A
    # homoeologous scaffold therefore aligns well and can outrank the true
    # orthologue on coverage alone.
    #
    # Signal used: take EVERY gene on a candidate scaffold, look up its
    # reciprocal best hit anywhere in the source genome, and see which source
    # chromosome those partners sit on. The orthologous scaffold's genes match
    # the QTL's own chromosome; a homoeolog's genes match a different member
    # of the same homoeologous group.
    # --- evidence 1: reciprocal best hits (needs the target annotation)
    tgt_gene_src_chrom = {}
    for q, lst in rbh.items():
        srec = src_by_id.get(q) or src_by_id.get(norm_id(q))
        if srec is None:
            continue
        s_seq = srec[0]
        for (t, _pid, _bits) in lst:
            trec = tgt_by_id.get(t) or tgt_by_id.get(norm_id(t))
            if trec is None:
                continue
            tgt_gene_src_chrom.setdefault(trec[4], s_seq)

    # --- evidence 2: miniprot placements (needs NO target annotation, so it
    # still works on scaffolds that are unannotated or badly annotated)
    #
    # Only a protein's BEST placement (Rank=1) is counted. Every member of a
    # homoeologous triad will align acceptably to all three copies, so
    # counting every hit would dilute the signal to noise; the best placement
    # is the one that actually discriminates.
    mp_by_seq = defaultdict(list)
    for h in mp_all_unfiltered:
        if h["rank"] != 1:
            continue
        srec = src_by_id.get(h["protein"]) or src_by_id.get(norm_id(h["protein"]))
        if srec is None:
            continue
        mp_by_seq[h["seqid"]].append((h["start"], h["end"], srec[0]))

    def chrom_tally(seqid, region):
        """(rbh Counter, miniprot Counter) of source chromosomes for one
        region -- scoped to THIS region row specifically (by rank), not to
        every gene that happens to share the scaffold name. Once a scaffold
        can produce more than one region (see rank_dna_regions.py's target-
        side clustering), a gene belonging to a distant, unrelated cluster
        on the same scaffold must not contribute to this region's tally."""
        c_rbh = Counter()
        for (sq, _s, _e, gid, _st, _tid, rk) in tgt_genes:
            if rk != region["rank"]:
                continue
            c = tgt_gene_src_chrom.get(gid)
            if c:
                c_rbh[c] += 1
        c_mp = Counter()
        for (hs, he, schrom) in mp_by_seq.get(seqid, []):
            if hs < region["tgt_end"] and he > region["tgt_start"]:
                c_mp[schrom] += 1
        return c_rbh, c_mp

    def homoeolog_call(region):
        seqid = region["tgt_seqid"]
        c_rbh, c_mp = chrom_tally(seqid, region)
        use = {"both": c_rbh + c_mp, "rbh": c_rbh, "miniprot": c_mp}[args.homoeolog_evidence]
        n_rbh, n_mp = sum(c_rbh.values()), sum(c_mp.values())
        if not use:
            return ("NA", 0, n_rbh, n_mp, 0.0, "UNKNOWN")
        top, n_top = use.most_common(1)[0]
        total = sum(use.values())
        frac = use.get(args.qtl_chrom, 0) / total
        if frac >= args.homoeolog_min_frac:
            call = "ORTHOLOG_LIKELY"
        elif top != args.qtl_chrom and n_top / total >= args.homoeolog_min_frac:
            call = f"HOMOEOLOG_LIKELY({top})"
        else:
            call = "UNCLEAR"
        # flag when the two independent lines of evidence disagree
        if n_rbh and n_mp:
            t_rbh = c_rbh.most_common(1)[0][0]
            t_mp = c_mp.most_common(1)[0][0]
            if t_rbh != t_mp:
                call += ";EVIDENCE_CONFLICT"
        return (top, total, n_rbh, n_mp, frac, call)


    for q, lst in rbh.items():
        srec = src_by_id.get(q) or src_by_id.get(norm_id(q))
        if srec is None:
            continue
        sseq, ss, se, sstrand, sgid = srec
        if sgid not in src_ids:
            continue          # only genes inside this QTL
        for (t, pid, bits) in lst:
            trec = tgt_by_id.get(t) or tgt_by_id.get(norm_id(t))
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

    # ---- primary-in-group, using PROTEIN evidence (not DNA coverage) -----
    # rank_dna_regions.py already grouped scaffolds whose covered source
    # interval overlaps substantially (possible genome copies). Which one is
    # the true match is decided HERE, now that gene evidence exists: the
    # member of each group with the most combined RBH + miniprot support
    # wins. Ties fall back to the DNA-coverage provisional flag.
    def evidence_count(seqid):
        c = per_region[seqid]
        return c["rbh"] + c["miniprot"]

    by_group = defaultdict(list)
    for r in regs:
        by_group[r["region_group"]].append(r)
    is_primary = {}
    for g, members in by_group.items():
        members.sort(key=lambda r: (-evidence_count(r["tgt_seqid"]),
                                    -r["provisional_group_best"],
                                    -r["qtl_cov_bp"]))
        for k, r in enumerate(members):
            is_primary[r["tgt_seqid"]] = (k == 0)

    # ---- coverage of the source split by gene evidence (two "logos") ------
    # Union coverage from ALL DNA hits (not just kept/displayed regions),
    # split into scaffolds that carry >=1 RBH or miniprot link vs those that
    # do not. Independent of what --top-n or the viewer ends up showing.
    paf_iv = read_paf_src_intervals(args.paf, args.qtl_start, args.qtl_end)
    with_gene_iv, without_gene_iv = [], []
    for seqid, blocks in paf_iv.items():
        bucket = with_gene_iv if evidence_count(seqid) > 0 else without_gene_iv
        bucket.extend(blocks)
    cov_with_bp, _ = merge_len(with_gene_iv)
    cov_without_bp, _ = merge_len(without_gene_iv)
    cov_total_bp, _ = merge_len(with_gene_iv + without_gene_iv)
    qtl_len = (args.qtl_end - args.qtl_start) or 1
    pct = lambda bp: 100.0 * bp / qtl_len

    with open(args.out_coverage, "w") as cov:
        cov.write("qtl_id\tqtl_len\tcov_with_genes_bp\tcov_with_genes_pct\t"
                  "cov_without_genes_bp\tcov_without_genes_pct\t"
                  "cov_total_bp\tcov_total_pct\n")
        cov.write(f"{args.qtl_id}\t{qtl_len}\t{cov_with_bp}\t{pct(cov_with_bp):.2f}\t"
                  f"{cov_without_bp}\t{pct(cov_without_bp):.2f}\t"
                  f"{cov_total_bp}\t{pct(cov_total_bp):.2f}\n")

    with open(args.out_regions, "w") as out:
        out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\tqtl_len\trank\t"
                  "tgt_seqid\taligned_bp\tn_aln_blocks\tpct_id\tstrand\t"
                  "qtl_cov_bp\tqtl_cov_pct\t"
                  "left_flank_cov_bp\tleft_flank_cov_pct\t"
                  "right_flank_cov_bp\tright_flank_cov_pct\t"
                  "full_query_cov_bp\t"
                  "tgt_cov_bp\tunion_qtl_cov_bp\tunion_qtl_cov_pct\t"
                  "tgt_start\ttgt_end\ttgt_span\tn_tgt_genes\tn_miniprot\tn_rbh\t"
                  "region_group\tgroup_size\tis_primary_in_group\tgroup_best_scaffold\t"
                  "consensus_src_chrom\tn_chrom_assigned\tn_chrom_rbh\tn_chrom_miniprot\t"
                  "frac_on_qtl_chrom\thomoeolog_call\n")
        for r in regs:
            ng = sum(1 for g in tgt_genes if g[0] == r["tgt_seqid"])
            c = per_region[r["tgt_seqid"]]
            hchrom, hn, h_nrbh, h_nmp, hfrac, hcall = homoeolog_call(r)
            if args.drop_homoeologs and hcall.startswith("HOMOEOLOG"):
                continue
            grp_members = by_group[r["region_group"]]
            grp_best = next(m["tgt_seqid"] for m in grp_members if is_primary[m["tgt_seqid"]])
            out.write(f"{args.qtl_id}\t{args.qtl_chrom}\t{args.qtl_start}\t{args.qtl_end}\t"
                      f"{r['qtl_len']}\t{r['rank']}\t{r['tgt_seqid']}\t"
                      f"{r['aligned_bp']}\t{r['n_aln_blocks']}\t{r['pct_id']:.1f}\t"
                      f"{r['strand']}\t{r['qtl_cov_bp']}\t{r['qtl_cov_pct']:.2f}\t"
                      f"{r['left_flank_cov_bp']}\t{r['left_flank_cov_pct']:.2f}\t"
                      f"{r['right_flank_cov_bp']}\t{r['right_flank_cov_pct']:.2f}\t"
                      f"{r['full_query_cov_bp']}\t"
                      f"{r['tgt_cov_bp']}\t{r['union_qtl_cov_bp']}\t"
                      f"{r['union_qtl_cov_pct']:.2f}\t"
                      f"{r['tgt_start']}\t{r['tgt_end']}\t{r['tgt_end']-r['tgt_start']}\t"
                      f"{ng}\t{c['miniprot']}\t{c['rbh']}\t"
                      f"{r['region_group']}\t{r['group_size']}\t"
                      f"{int(is_primary[r['tgt_seqid']])}\t{grp_best}\t"
                      f"{hchrom}\t{hn}\t{h_nrbh}\t{h_nmp}\t{hfrac:.2f}\t{hcall}\n")

    n_mp = sum(1 for l in links if l["track"] == "miniprot")
    n_rb = sum(1 for l in links if l["track"] == "rbh")
    calls = [homoeolog_call(r)[5] for r in regs]
    n_hom = sum(1 for c in calls if c.startswith("HOMOEOLOG"))
    n_conf = sum(1 for c in calls if "EVIDENCE_CONFLICT" in c)
    n_groups_multi = sum(1 for g, m in by_group.items() if len(m) > 1)
    if n_groups_multi:
        sys.stderr.write(
            f"[annotate_regions] {args.qtl_id}: {n_groups_multi} paralogy group(s) "
            f"have more than one candidate scaffold (possible genome copies); "
            f"the best match in each was chosen by RBH+miniprot evidence, others "
            f"kept and flagged - see region_group/is_primary_in_group\n")
    sys.stderr.write(
        f"[annotate_regions] {args.qtl_id}: source coverage - "
        f"{pct(cov_with_bp):.1f}% by gene-bearing scaffolds, "
        f"{pct(cov_without_bp):.1f}% by gene-less scaffolds, "
        f"{pct(cov_total_bp):.1f}% total\n")
    if n_hom:
        sys.stderr.write(
            f"[annotate_regions] {args.qtl_id}: {n_hom}/{len(regs)} regions look "
            f"HOMOEOLOGOUS (evidence={args.homoeolog_evidence}; their proteins point at "
            f"a source chromosome other than {args.qtl_chrom}) - see homoeolog_call\n")
    if n_conf:
        sys.stderr.write(
            f"[annotate_regions] {args.qtl_id}: {n_conf} region(s) where RBH and miniprot "
            f"disagree on the source chromosome - flagged EVIDENCE_CONFLICT\n")
    sys.stderr.write(
        f"[annotate_regions] {args.qtl_id}: {len(regs)} regions, "
        f"{len(src_genes)} source genes in QTL, {len(tgt_genes)} target genes in regions, "
        f"{n_mp}/{n_mp_total} miniprot hits in-region, {n_rb}/{n_rbh_total} RBH in-region\n")

    # Distinguish "no evidence" from "ids do not join" - these need very
    # different fixes and look identical in the viewer.
    if args.miniprot and n_mp_total == 0 and mp_seen:
        sys.stderr.write(
            f"[annotate_regions] WARNING: {mp_seen} miniprot hits in the file but none "
            f"matched a gene id for this QTL. Example miniprot Target={mp_example!r}; "
            f"example source gene id={next(iter(src_ids), None)!r}. Run bin/check_ids.py.\n")
    if args.rbh and n_rbh_total and n_rb == 0 and rbh_resolved == 0:
        sys.stderr.write(
            f"[annotate_regions] WARNING: {n_rbh_total} RBH pairs in the file but none "
            f"resolved to a gene. Example rbh query={rbh_example!r}; "
            f"example source gene id={next(iter(src_ids), None)!r}. Run bin/check_ids.py.\n")


if __name__ == "__main__":
    main()
