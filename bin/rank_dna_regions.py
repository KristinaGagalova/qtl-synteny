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


def union_add(merged, new):
    """Merge `new` intervals into sorted, disjoint `merged`.
    Returns (merged, gained_bp) - how much NEW territory `new` added."""
    before = sum(e - s for s, e in merged)
    iv = sorted(merged + [tuple(x) for x in new])
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out, sum(e - s for s, e in out) - before


def merge_len(intervals):
    """Total bp covered by a set of [start, end) intervals, counting any
    overlap once. Summing raw block lengths double-counts repeats and
    tandem hits, which inflates coverage badly on a repeat-rich genome."""
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True, help="source QTL slice vs whole target")
    ap.add_argument("--qtl-id", required=True)
    ap.add_argument("--src-chrom", required=True)
    ap.add_argument("--src-start", type=int, required=True)
    ap.add_argument("--src-end", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-marginal-pct", type=float, default=None,
                    help="GREEDY MODE (recommended). Repeatedly add the target "
                         "that contributes the most NEW source coverage, and "
                         "stop when the best remaining one would add less than "
                         "this percent. Handles redundancy: scaffolds covering "
                         "ground already explained are never added.")
    ap.add_argument("--target-cum-cov-pct", type=float, default=None,
                    help="greedy mode: also stop once cumulative union "
                         "coverage reaches this percent")
    ap.add_argument("--min-pct-id", type=float, default=0.0,
                    help="drop targets whose alignment identity is below this "
                         "before ranking (a homoeolog from another subgenome "
                         "often covers well but at lower identity)")
    ap.add_argument("--weight-identity", action="store_true",
                    help="greedy mode: score candidates by marginal gain "
                         "scaled by identity, not marginal gain alone")
    ap.add_argument("--top-n", type=int, default=5,
                    help="fixed count; used only when neither "
                         "--min-marginal-pct nor --min-src-cov-pct is set")
    ap.add_argument("--min-src-cov-pct", type=float, default=None,
                    help="per-scaffold floor: keep EVERY target covering at "
                         "least this percent of the source, judged "
                         "independently (no redundancy handling)")
    ap.add_argument("--max-regions", type=int, default=40,
                    help="hard ceiling on regions kept, so a pathological QTL "
                         "cannot render hundreds of lanes")
    ap.add_argument("--min-aligned-bp", type=int, default=2000)
    ap.add_argument("--min-mapq", type=int, default=0)
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
                # keep the QUERY span too: coverage of the source interval is
                # measured on the query side, not the target side
                by_t[f[5]].append((int(f[7]), int(f[8]), int(f[10]), int(f[9]),
                                   f[4], int(f[2]), int(f[3])))
                qlen_seen.add(int(f[1]))
            except ValueError:
                continue

    # length of the source slice that was aligned (the query)
    src_slice_len = max(qlen_seen) if qlen_seen else (args.src_end - args.src_start)

    rows = []
    all_src_iv = []
    for seqid, blocks in by_t.items():
        aligned_bp = sum(b[2] for b in blocks)
        if aligned_bp < args.min_aligned_bp:
            continue
        matched = sum(b[3] for b in blocks)
        strands = defaultdict(int)
        for b in blocks:
            strands[b[4]] += b[2]

        # how much of the SOURCE interval this one target sequence covers,
        # overlapping blocks counted once
        src_iv = [(b[5], b[6]) for b in blocks]
        src_cov, _ = merge_len(src_iv)
        all_src_iv.extend(src_iv)
        tgt_cov, _ = merge_len([(b[0], b[1]) for b in blocks])

        rows.append(dict(
            tgt_seqid=seqid, aligned_bp=aligned_bp, n_aln_blocks=len(blocks),
            pct_id=100.0 * matched / aligned_bp if aligned_bp else 0.0,
            src_cov_bp=src_cov,
            src_cov_pct=100.0 * src_cov / src_slice_len if src_slice_len else 0.0,
            tgt_cov_bp=tgt_cov,
            tgt_start=max(min(b[0] for b in blocks) - args.flank, 0),
            tgt_end=max(b[1] for b in blocks) + args.flank,
            strand=max(strands, key=strands.get)))

    # rank by how much of the source each target actually explains
    # identity filter runs before any ranking
    if args.min_pct_id > 0:
        n0 = len(rows)
        rows = [r for r in rows if r["pct_id"] >= args.min_pct_id]
        if n0 != len(rows):
            sys.stderr.write(f"[rank_dna_regions] {args.qtl_id}: dropped "
                             f"{n0-len(rows)} targets below {args.min_pct_id}% identity\n")

    rows.sort(key=lambda r: (-r["src_cov_bp"], -r["pct_id"]))
    src_iv_of = {r["tgt_seqid"]: [(b[5], b[6]) for b in by_t[r["tgt_seqid"]]] for r in rows}

    # ---- selection ---------------------------------------------------
    marginal = {}
    if args.min_marginal_pct is not None:
        # GREEDY SET COVER. Re-scores every remaining candidate against what
        # is already covered, each round. This is the difference from a fixed
        # ranking: a target that looks strong on its own is skipped if it
        # merely repeats ground an earlier one already explained.
        thresh_bp = args.min_marginal_pct / 100.0 * src_slice_len
        target_bp = (args.target_cum_cov_pct / 100.0 * src_slice_len
                     if args.target_cum_cov_pct else None)
        remaining = list(rows)
        covered = []
        selected = []
        while remaining:
            if args.max_regions > 0 and len(selected) >= args.max_regions:
                stop = f"hit --max-regions {args.max_regions}"
                break
            best, best_gain, best_score = None, 0, -1.0
            for r in remaining:
                _m, gain = union_add(covered, src_iv_of[r["tgt_seqid"]])
                score = gain * (r["pct_id"] / 100.0) if args.weight_identity else gain
                if score > best_score or (score == best_score and best is not None
                                          and r["pct_id"] > best["pct_id"]):
                    best, best_gain, best_score = r, gain, score
            if best is None:
                stop = "no candidates left"
                break
            if selected and best_gain < thresh_bp:
                stop = (f"next best would add only "
                        f"{100.0*best_gain/src_slice_len if src_slice_len else 0:.2f}% "
                        f"(< --min-marginal-pct {args.min_marginal_pct})")
                break
            covered, gained = union_add(covered, src_iv_of[best["tgt_seqid"]])
            marginal[best["tgt_seqid"]] = gained
            selected.append(best)
            remaining.remove(best)
            if target_bp is not None and sum(e - s for s, e in covered) >= target_bp:
                stop = f"reached --target-cum-cov-pct {args.target_cum_cov_pct}"
                break
        else:
            stop = "all targets used"
        mode = f"greedy marginal >= {args.min_marginal_pct}% ({stop})"

    elif args.min_src_cov_pct is not None and args.min_src_cov_pct > 0:
        selected = [r for r in rows if r["src_cov_pct"] >= args.min_src_cov_pct]
        mode = f"src_cov_pct >= {args.min_src_cov_pct}"
        if not selected and rows:
            selected = rows[:1]
            mode += " (none passed; kept best)"
    else:
        selected = rows[: args.top_n] if args.top_n > 0 else rows
        mode = f"top {args.top_n}"

    n_before_cap = len(selected)
    if args.max_regions > 0 and len(selected) > args.max_regions:
        selected = selected[: args.max_regions]

    kept_seqids = {r["tgt_seqid"] for r in selected}
    kept_iv = []
    for seqid, blocks in by_t.items():
        if seqid in kept_seqids:
            kept_iv.extend([(b[5], b[6]) for b in blocks])
    union_all, _ = merge_len(all_src_iv)
    union_top, _ = merge_len(kept_iv)
    rows = selected

    with open(args.out, "w") as out:
        out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\tsrc_slice_len\trank\t"
                  "tgt_seqid\taligned_bp\tn_aln_blocks\tpct_id\tstrand\t"
                  "src_cov_bp\tsrc_cov_pct\tmarginal_cov_bp\tmarginal_cov_pct\t"
                  "cum_src_cov_bp\tcum_src_cov_pct\t"
                  "tgt_cov_bp\ttgt_start\ttgt_end\ttgt_span\t"
                  "union_src_cov_bp\tunion_src_cov_pct\n")
        cum_iv = []
        for rank, r in enumerate(rows, start=1):
            cum_iv.extend([(b[5], b[6]) for b in by_t[r["tgt_seqid"]]])
            cum_cov, _ = merge_len(cum_iv)
            cum_pct = 100.0 * cum_cov / src_slice_len if src_slice_len else 0.0
            out.write(f"{args.qtl_id}\t{args.src_chrom}\t{args.src_start}\t"
                      f"{args.src_end}\t{src_slice_len}\t{rank}\t{r['tgt_seqid']}\t"
                      f"{r['aligned_bp']}\t{r['n_aln_blocks']}\t{r['pct_id']:.1f}\t"
                      f"{r['strand']}\t{r['src_cov_bp']}\t{r['src_cov_pct']:.2f}\t"
                      f"{marginal.get(r['tgt_seqid'], 0)}\t"
                      f"{100.0*marginal.get(r['tgt_seqid'], 0)/src_slice_len if src_slice_len else 0:.2f}\t"
                      f"{cum_cov}\t{cum_pct:.2f}\t{r['tgt_cov_bp']}\t"
                      f"{r['tgt_start']}\t{r['tgt_end']}\t{r['tgt_end']-r['tgt_start']}\t"
                      f"{union_top}\t{100.0*union_top/src_slice_len if src_slice_len else 0:.2f}\n")

    pct_top = 100.0 * union_top / src_slice_len if src_slice_len else 0.0
    pct_all = 100.0 * union_all / src_slice_len if src_slice_len else 0.0
    capped = (f", capped from {n_before_cap} by --max-regions {args.max_regions}"
              if n_before_cap > len(rows) else "")
    sys.stderr.write(
        f"[rank_dna_regions] {args.qtl_id}: {len(by_t)} target seqs hit, "
        f"{len(rows)} kept by {mode}{capped}; source covered "
        f"{union_top:,}/{src_slice_len:,} bp ({pct_top:.1f}%) by displayed targets, "
        f"{union_all:,} bp ({pct_all:.1f}%) by all\n")
    if pct_all - pct_top > 5.0:
        sys.stderr.write(
            f"[rank_dna_regions] {args.qtl_id}: NOTE {pct_all-pct_top:.1f}% of the source "
            f"is covered only by targets that were not displayed - loosen "
            f"--min-src-cov-pct or raise --max-regions to see them\n")


if __name__ == "__main__":
    main()
