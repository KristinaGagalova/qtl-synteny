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
                src_cov_pct=f[i["src_cov_pct"]] if "src_cov_pct" in i else "NA",
                pct_id=f[i["pct_id"]] if "pct_id" in i else "NA",
                call=call))
    regs.sort(key=lambda r: r["rank"])
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
    args = ap.parse_args()

    regs = read_regions(args.regions, args.qtl_id, args.drop_homoeologs)
    wanted = {r["seqid"]: r for r in regs}

    with open(args.out_list, "w") as lo:
        lo.write("qtl_id\trank\ttgt_seqid\ttgt_start\ttgt_end\tsrc_cov_pct\t"
                 "pct_id\thomoeolog_call\n")
        for r in regs:
            lo.write(f"{args.qtl_id}\t{r['rank']}\t{r['seqid']}\t{r['start']}\t"
                     f"{r['end']}\t{r['src_cov_pct']}\t{r['pct_id']}\t{r['call']}\n")

    if not regs:
        open(args.out_full, "w").close()
        open(args.out_window, "w").close()
        sys.stderr.write(f"[extract_selected_scaffolds] {args.qtl_id}: no regions\n")
        return

    def wrap(fo, seq):
        for k in range(0, len(seq), 60):
            fo.write(seq[k:k + 60] + "\n")

    n_full = n_win = 0
    with open(args.target_fasta) as fh, \
         open(args.out_full, "w") as ff, open(args.out_window, "w") as fw:
        name, buf, keep = None, [], False

        def flush():
            nonlocal n_full, n_win
            if not keep or name not in wanted:
                return
            r = wanted[name]
            seq = "".join(buf)
            tag = (f"rank={r['rank']} qtl={args.qtl_id} "
                   f"src_cov_pct={r['src_cov_pct']} pct_id={r['pct_id']} "
                   f"call={r['call']}")
            ff.write(f">{name} {tag} len={len(seq)}\n")
            wrap(ff, seq)
            n_full += 1
            s2, e2 = max(r["start"], 0), min(r["end"], len(seq))
            if e2 > s2:
                fw.write(f">{name}:{s2+1}-{e2} {tag}\n")
                wrap(fw, seq[s2:e2])
                n_win += 1

        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                keep = name in wanted
                buf = []
            elif keep:
                buf.append(line.strip())
        flush()

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
