#!/usr/bin/env python3
"""
final_gene_table.py -- one flat row per gene (source AND target), for one
QTL, consolidating everything the pipeline already computed:

    QTLname  gene  source_target  chr_name  start  stop  strand
    rbh_available  protein_remapping  homeolog_call
    avg_gene_expression  <annotation columns...>

Nothing new is computed here -- this reads genes.tsv, links.tsv and
regions.tsv (from annotate_regions.py) plus the raw expression and
annotation tables you supplied, and flattens them into one table meant to
leave the pipeline: for spreadsheets, filtering, or joining elsewhere.

Column semantics
-----------------
source_target        'source' or 'target'
rbh_available         Y/N -- SOURCE rows only (RBH is defined source->target;
                      a target row has no "does this gene have an RBH"
                      question of its own, it either IS one or isn't)
protein_remapping     SOURCE rows only: every target gene this source gene's
                      protein evidence (RBH and/or miniprot) pointed at,
                      semicolon-separated. Blank for target rows.
homeolog_call         TARGET rows only: the homoeolog_call of the region
                      (scaffold) this gene sits on -- a property of the
                      target scaffold, not of an individual source gene.
avg_gene_expression   mean of that gene's own expression row (its own side's
                      table), skipping missing values. Blank if no match.
annotation columns    union of the source and target annotation files'
                      columns (whatever the user supplied), filled from
                      whichever side's table matched this gene; blank where
                      that column does not apply to this side.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

TYPE_PREFIX = re.compile(
    r"^(?:gene|transcript|mRNA|CDS|protein|exon|rna|ncRNA):", re.IGNORECASE)
TX_SUFFIX = re.compile(r"\.\d+$")


def norm_id(x):
    if not x:
        return x
    return TX_SUFFIX.sub("", TYPE_PREFIX.sub("", x, count=1))


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def read_expression(path):
    """gene_id -> list of numeric values (None for non-numeric cells)."""
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
        hdr = first.split("\t") if "\t" in first else first.split()
        n = len(hdr) - 1
        data = {}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t") if "\t" in line else line.split()
            vals = []
            for x in f[1:n + 1]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(None)
            data[f[0]] = vals
    return data


def match_by_id(table, gene_id):
    """Exact id, then normalised -- same tolerant join used throughout."""
    if gene_id in table:
        return table[gene_id]
    return table.get(norm_id(gene_id))


def avg_expr(expr_table, gene_id):
    vals = match_by_id(expr_table, gene_id)
    if not vals:
        return ""
    nums = [v for v in vals if v is not None]
    return f"{sum(nums)/len(nums):.4f}" if nums else ""


def read_annotation(path):
    if not path or not os.path.exists(path):
        return [], {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        columns = hdr[1:]
        data = {}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            vals = f[1:len(columns) + 1]
            while len(vals) < len(columns):
                vals.append("")
            data[f[0]] = vals
    return columns, data


def region_passes(r, mode, min_gene_region_bp):
    """Decide whether one target scaffold gets a lane / a row.

    'evidence'    only scaffolds with >=1 RBH or miniprot hit (default)
    'gene_length' the above, PLUS any scaffold that carries >=1 annotated
                  gene and whose aligned length clears --min-gene-region-bp
                  -- for real synteny that annotation or protein evidence
                  missed (see README: "synteny without RBH or miniprot")
    'all'         every scaffold that aligned at all, no filtering

    Both build_synteny_html.py and final_gene_table.py MUST apply this
    exact rule with the exact same arguments, or the HTML view and the
    gene table list different genes again.
    """
    if mode == "all":
        return True
    has_evidence = (int(r.get("n_miniprot", 0) or 0) > 0
                    or int(r.get("n_rbh", 0) or 0) > 0)
    if mode == "evidence":
        return has_evidence
    if mode == "gene_length":
        has_gene = int(r.get("n_tgt_genes", 0) or 0) > 0
        aligned = int(r.get("aligned_bp", 0) or 0)
        return has_evidence or (has_gene and aligned >= min_gene_region_bp)
    return has_evidence

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qtl-id", required=True)
    ap.add_argument("--genes", required=True, help="annotate_regions.py --out-genes")
    ap.add_argument("--links", required=True, help="annotate_regions.py --out-links")
    ap.add_argument("--regions", required=True, help="annotate_regions.py --out-regions")
    ap.add_argument("--source-expr")
    ap.add_argument("--target-expr")
    ap.add_argument("--source-annotation")
    ap.add_argument("--target-annotation")
    ap.add_argument("--out", required=True)
    ap.add_argument("--region-display-mode", choices=["evidence", "gene_length", "all"],
                    default="evidence",
                    help="MUST match build_synteny_html.py's setting, or "
                         "this table and the HTML views list different "
                         "genes. See build_synteny_html.py --help.")
    ap.add_argument("--min-gene-region-bp", type=int, default=5000)
    args = ap.parse_args()

    genes = [g for g in read_tsv(args.genes) if g["qtl_id"] == args.qtl_id]
    links = [l for l in read_tsv(args.links) if l["qtl_id"] == args.qtl_id]
    regions = [r for r in read_tsv(args.regions) if r["qtl_id"] == args.qtl_id]

    # source gene -> {rbh: bool, targets: set(tgt_gene)}
    remap = defaultdict(lambda: {"rbh": False, "targets": set()})
    for l in links:
        sg = l["src_gene"]
        if l["track"] == "rbh":
            remap[sg]["rbh"] = True
        if l["tgt_gene"]:
            remap[sg]["targets"].add(l["tgt_gene"])

    # target scaffold -> homoeolog_call
    homoeolog_by_seqid = {r["tgt_seqid"]: r.get("homoeolog_call", "NA") for r in regions}

    # Match the viewer's lane visibility EXACTLY (region_passes() above) so
    # this table and the HTML views always list the same genes.
    if args.region_display_mode == "all":
        shown_seqids = None  # no filtering
    else:
        shown_seqids = {r["tgt_seqid"] for r in regions
                        if region_passes(r, args.region_display_mode, args.min_gene_region_bp)}
        n_all_tgt_seqids = {r["tgt_seqid"] for r in regions}
        n_hidden_seqids = n_all_tgt_seqids - shown_seqids
        if n_hidden_seqids:
            sys.stderr.write(
                f"[final_gene_table] {args.qtl_id}: {len(n_hidden_seqids)} "
                f"scaffold(s) excluded by --region-display-mode "
                f"{args.region_display_mode!r}: "
                f"{', '.join(sorted(n_hidden_seqids)[:5])}"
                f"{' ...' if len(n_hidden_seqids) > 5 else ''}\n")

    expr_src = read_expression(args.source_expr)
    expr_tgt = read_expression(args.target_expr)
    cols_src, ann_src = read_annotation(args.source_annotation)
    cols_tgt, ann_tgt = read_annotation(args.target_annotation)
    ann_columns = list(cols_src) + [c for c in cols_tgt if c not in cols_src]

    with open(args.out, "w") as out:
        out.write("\t".join([
            "QTLname", "gene", "source_target", "chr_name", "start", "stop",
            "strand", "rbh_available", "protein_remapping", "homeolog_call",
            "avg_gene_expression"] + ann_columns) + "\n")

        n_src = n_tgt = n_skipped = 0
        for g in genes:
            side = g["side"]
            gid = g["gene_id"]

            if side == "target" and shown_seqids is not None \
                    and g["seqid"] not in shown_seqids:
                n_skipped += 1
                continue

            if side == "source":
                n_src += 1
                rbh_available = "Y" if remap[gid]["rbh"] else "N"
                protein_remapping = ";".join(sorted(remap[gid]["targets"]))
                homeolog_call = ""
                expr = avg_expr(expr_src, gid)
                ann_hit = match_by_id(ann_src, gid)
            else:
                n_tgt += 1
                rbh_available = ""
                protein_remapping = ""
                homeolog_call = homoeolog_by_seqid.get(g["seqid"], "NA")
                expr = avg_expr(expr_tgt, gid)
                ann_hit = match_by_id(ann_tgt, gid)

            side_cols = cols_src if side == "source" else cols_tgt
            ann_row = []
            for c in ann_columns:
                if ann_hit is not None and c in side_cols:
                    ann_row.append(ann_hit[side_cols.index(c)])
                else:
                    ann_row.append("")

            out.write("\t".join([
                args.qtl_id, gid, side, g["seqid"], g["start"], g["end"],
                g["strand"], rbh_available, protein_remapping, homeolog_call,
                expr] + ann_row) + "\n")

    sys.stderr.write(f"[final_gene_table] {args.qtl_id}: {n_src} source genes, "
                     f"{n_tgt} target genes ({n_skipped} skipped as gene-less-"
                     f"scaffold), {len(ann_columns)} annotation column(s) "
                     f"-> {args.out}\n")


if __name__ == "__main__":
    main()
