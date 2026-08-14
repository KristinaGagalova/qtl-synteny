#!/usr/bin/env python3
"""
build_index_html.py -- landing page linking every per-QTL synteny view,
with the candidate-region summary inline so you can triage before clicking.
"""
import argparse, html, os, sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--qtl", required=True)
    ap.add_argument("--views", nargs="*", default=[], help="per-QTL html filenames")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    qtls = []
    with open(args.qtl) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            qtls.append((f[3] if len(f) > 3 else f"{f[0]}:{f[1]}-{f[2]}",
                         f[0], int(f[1]), int(f[2])))

    reg = defaultdict(list)
    if os.path.exists(args.regions):
        with open(args.regions) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if not line.strip():
                    continue
                r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                reg[r["qtl_id"]].append(r)

    # Views are named <safe_qtl_id>.html by BUILD_SYNTENY_HTML, so match on
    # that exact stem rather than guessing from a prefix.
    def safe_id(qid):
        return "".join(c if c.isalnum() or c in "._-" else "_" for c in qid)

    views = {os.path.basename(v) for v in args.views}

    def view_for(qid):
        want = safe_id(qid) + ".html"
        return want if want in views else None

    rows = []
    for qid, chrom, s_, e_ in qtls:
        rs = sorted(reg.get(qid, []), key=lambda x: int(x["rank"]))
        rs_gene = [r for r in rs
                  if int(r.get("n_miniprot", 0) or 0) or int(r.get("n_rbh", 0) or 0)]
        v = view_for(qid)
        label = (f'<a href="{html.escape(v)}">{html.escape(qid)}</a>'
                 if v else html.escape(qid))
        n_groups_multi = len({r["region_group"] for r in rs
                              if int(r.get("group_size", 1) or 1) > 1})
        if rs_gene:
            top = rs_gene[0]
            tgt  = html.escape(top["tgt_seqid"])
            cov  = f'{float(top.get("union_src_cov_pct", 0) or 0):.1f}%'
            tcov = f'{float(top.get("src_cov_pct", 0) or 0):.1f}%'
            pid  = f'{float(top.get("pct_id", 0)):.1f}%'
            ngen = top.get("n_tgt_genes", "0")
            nmp  = top.get("n_miniprot", "0")
            nrbh = top.get("n_rbh", "0")
        else:
            tgt = cov = tcov = pid = ngen = nmp = nrbh = "&mdash;"
        rows.append(
            "<tr>"
            f"<td>{label}</td>"
            f"<td>{html.escape(chrom)}</td>"
            f'<td class="num">{s_:,}</td>'
            f'<td class="num">{e_:,}</td>'
            f'<td class="num">{(e_-s_)/1e6:.2f} Mb</td>'
            f'<td class="num">{len(rs)}</td>'
            f"<td>{tgt}</td>"
            f'<td class="num"><b>{cov}</b></td>'
            f'<td class="num">{tcov}</td>'
            f'<td class="num">{pid}</td>'
            f'<td class="num">{ngen}</td>'
            f'<td class="num">{nmp}</td>'
            f'<td class="num">{nrbh}</td>'
            "</tr>")

    body = "\n".join(rows)
    out = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>QTL synteny views</title><style>
body{{margin:0;font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;color:#1c2027}}
header{{padding:16px 20px;border-bottom:1px solid #d8dce3;background:#f7f8fa}}
h1{{margin:0 0 4px;font-size:17px}} .sub{{color:#697386;font-size:12px}}
.wrap{{padding:16px 20px;max-width:1400px}}
table{{border-collapse:collapse;font-size:12.5px;width:100%}}
th,td{{border-bottom:1px solid #d8dce3;padding:6px 10px;text-align:left;white-space:nowrap}}
th{{background:#f7f8fa;font-weight:600}} td.num{{text-align:right;font-variant-numeric:tabular-nums}}
a{{color:#2f6f9f}} tr:hover{{background:#f7f8fa}}
</style></head><body>
<header><h1>QTL synteny views</h1>
<div class="sub">{len(qtls)} QTL intervals &middot; &quot;QTL covered&quot; is the union of all displayed targets over the source interval (overlaps counted once); &quot;top covers&quot; is the best single target &middot; click a QTL to open its view</div></header>
<div class="wrap"><table>
<thead><tr><th>QTL</th><th>chrom</th><th>start</th><th>end</th><th>length</th>
<th>regions</th><th>top target</th><th>QTL covered</th><th>top covers</th><th>% id</th>
<th>genes</th><th>miniprot</th><th>RBH</th></tr></thead>
<tbody>
{body}
</tbody></table></div></body></html>"""
    with open(args.out, "w") as fh:
        fh.write(out)
    sys.stderr.write(f"[build_index_html] {args.out} ({len(qtls)} QTLs)\n")


if __name__ == "__main__":
    main()
