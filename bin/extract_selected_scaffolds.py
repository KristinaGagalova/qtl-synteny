#!/usr/bin/env python3
"""
extract_selected_scaffolds.py -- write the target sequences that were
selected for a QTL to FASTA, for downstream work.

Two outputs, because both are useful for different things:
  --out-full    the WHOLE selected scaffolds, unmodified. Use this to
                re-annotate, BLAST, or assemble against.
  --out-window  only the displayed window (region start..end, i.e. the
                aligned span plus --region_flank). Smaller, matches exactly
                what the viewer shows.

Headers carry the evidence, so the provenance survives into whatever you do
next:
  >scafX rank=1 qtl=1B1 src_cov_pct=18.0 pct_id=99.2 call=ORTHOLOG_LIKELY
"""
import argparse
import sys


def read_regions(path, qtl_id, drop_homoeologs=False):
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
            call = f[i["homoeolog_call"]] if "homoeolog_call" in i else "NA"
            if drop_homoeologs and call.startswith("HOMOEOLOG"):
                continue
            regs.append(dict(
                seqid=f[i["tgt_seqid"]], rank=int(f[i["rank"]]),
                start=int(f[i["tgt_start"]]), end=int(f[i["tgt_end"]]),
                src_cov_pct=f[i["qtl_cov_pct"]] if "qtl_cov_pct" in i else "NA",
                pct_id=f[i["pct_id"]] if "pct_id" in i else "NA",
                n_tgt_genes=int(f[i["n_tgt_genes"]]) if "n_tgt_genes" in i else 0,
                region_group=f[i["region_group"]] if "region_group" in i else "NA",
                group_size=f[i["group_size"]] if "group_size" in i else "1",
                is_primary=f[i["is_primary_in_group"]] if "is_primary_in_group" in i else "1",
                call=call))
    return regs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--qtl-id", required=True)
    ap.add_argument("--target-fasta", required=True)
    ap.add_argument("--out-full", required=True)
    ap.add_argument("--out-window", required=True)
    ap.add_argument("--out-list", required=True)
    ap.add_argument("--drop-homoeologs", action="store_true")
    ap.add_argument("--gene-priority", dest="gene_priority", action="store_true",
                    default=True,
                    help="list/write scaffolds that carry annotated genes "
                         "before gene-less ones (default: on)")
    ap.add_argument("--no-gene-priority", dest="gene_priority", action="store_false")
    args = ap.parse_args()

    regs = read_regions(args.regions, args.qtl_id, args.drop_homoeologs)
    if args.gene_priority:
        # scaffolds carrying genes first, ties broken by the original rank
        # (i.e. coverage/greedy order) so priority is genes, then evidence
        regs.sort(key=lambda r: (r["n_tgt_genes"] == 0, r["rank"]))
    else:
        regs.sort(key=lambda r: r["rank"])
    wanted = {r["seqid"]: r for r in regs}

    with open(args.out_list, "w") as lo:
        lo.write("list_order\tqtl_id\trank\ttgt_seqid\ttgt_start\ttgt_end\t"
                 "n_tgt_genes\tsrc_cov_pct\tpct_id\tregion_group\tgroup_size\t"
                 "is_primary_in_group\thomoeolog_call\n")
        for order, r in enumerate(regs, start=1):
            lo.write(f"{order}\t{args.qtl_id}\t{r['rank']}\t{r['seqid']}\t{r['start']}\t"
                     f"{r['end']}\t{r['n_tgt_genes']}\t{r['src_cov_pct']}\t{r['pct_id']}\t"
                     f"{r['region_group']}\t{r['group_size']}\t{r['is_primary']}\t"
                     f"{r['call']}\n")

    if not regs:
        open(args.out_full, "w").close()
        open(args.out_window, "w").close()
        sys.stderr.write(f"[extract_selected_scaffolds] {args.qtl_id}: no regions\n")
        return

    def wrap(fo, seq):
        for k in range(0, len(seq), 60):
            fo.write(seq[k:k + 60] + "\n")

    n_full = n_win = 0

    # Stream the genome once, but BUFFER the wanted sequences and write them
    # out afterwards in the sorted `regs` order (gene-priority first) rather
    # than genome-file order - otherwise the FASTA order would silently
    # ignore --gene-priority even though the .tsv respects it.
    seqs = {}
    with open(args.target_fasta) as fh:
        name, buf, keep = None, [], False

        def flush():
            if keep and name in wanted:
                seqs[name] = "".join(buf)

        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                keep = name in wanted
                buf = []
            elif keep:
                buf.append(line.strip())
        flush()

    def wrap(fo, seq):
        for k in range(0, len(seq), 60):
            fo.write(seq[k:k + 60] + "\n")

    with open(args.out_full, "w") as ff, open(args.out_window, "w") as fw:
        for r in regs:
            seq = seqs.get(r["seqid"])
            if seq is None:
                continue
            primary_tag = "primary" if r['is_primary'] == "1" else "alt_copy"
            tag = (f"rank={r['rank']} qtl={args.qtl_id} "
                   f"n_genes={r['n_tgt_genes']} "
                   f"src_cov_pct={r['src_cov_pct']} pct_id={r['pct_id']} "
                   f"group={r['region_group']}/{r['group_size']}({primary_tag}) "
                   f"call={r['call']}")
            ff.write(f">{r['seqid']} {tag} len={len(seq)}\n")
            wrap(ff, seq)
            n_full += 1
            s2, e2 = max(r["start"], 0), min(r["end"], len(seq))
            if e2 > s2:
                fw.write(f">{r['seqid']}:{s2+1}-{e2} {tag}\n")
                wrap(fw, seq[s2:e2])
                n_win += 1

    missing = set(wanted) - set()
    sys.stderr.write(f"[extract_selected_scaffolds] {args.qtl_id}: "
                     f"{n_full} scaffolds written ({len(regs)} selected), "
                     f"{n_win} windows\n")
    if n_full < len(regs):
        sys.stderr.write(f"[extract_selected_scaffolds] WARNING: "
                         f"{len(regs)-n_full} selected scaffold(s) not found in "
                         f"{args.target_fasta}\n")


if __name__ == "__main__":
    main()
