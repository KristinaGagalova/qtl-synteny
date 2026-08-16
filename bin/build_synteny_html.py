#!/usr/bin/env python3
"""
build_synteny_html.py -- one self-contained interactive HTML page per QTL.

No CDN, no network: all data is embedded as JSON and all JS/CSS is inline,
so the pages work on an offline HPC filesystem and can be copied anywhere.

Panels
------
1. Synteny    source QTL region on top (gene arrows to scale), each candidate
              target scaffold below it. DNA alignment blocks from minimap2 are
              drawn as ribbons; protein links between genes are drawn as
              curves, with reciprocal-best-hit pairs highlighted.
2. Expression heatmap, genes on X and samples on Y, for source and target
              genes side by side. Click any gene to pin it.
3. Table      per-gene expression values for whatever is selected.

Everything is driven by the same selection state, so clicking a gene in the
synteny panel highlights it in the heatmap and filters the table.
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict

TYPE_PREFIX = re.compile(
    r"^(?:gene|transcript|mRNA|CDS|protein|exon|rna|ncRNA):", re.IGNORECASE)
TX_SUFFIX = re.compile(r"\.\d+$")


def norm_id(x):
    """Canonical id form; see annotate_regions.py. Both sides of the join are
    normalised, which is the only way a bare id and a prefixed one meet."""
    if not x:
        return x
    return TX_SUFFIX.sub("", TYPE_PREFIX.sub("", x, count=1))


# --------------------------------------------------------------------------
def read_tsv(path, required=None):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if required:
            missing = [c for c in required if c not in hdr]
            if missing:
                sys.exit(f"[build_synteny_html] {path}: missing columns {missing}")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return rows


def read_expression(path):
    """Counts table: first column = gene id, remaining columns = samples
    (genes down the rows, samples across the header).

    Returns (samples, {gene: [values]}). Non-numeric cells become None.
    Whitespace-delimited so tab- and space-separated files both work.
    """
    if not path or not os.path.exists(path):
        return [], {}
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
        hdr = first.split("\t") if "\t" in first else first.split()
        samples = hdr[1:]
        data = {}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t") if "\t" in line else line.split()
            vals = []
            for x in f[1:len(samples) + 1]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(None)
            while len(vals) < len(samples):
                vals.append(None)
            data[f[0]] = vals
    return samples, data


MATCH_LOG = []   # accumulated match diagnostics, written to --match-report


def report_match(kind, label, genes, table, matched, missed):
    """One consistent diagnostic for every gene->table join.

    Both the counts table and the annotation table are joined to the gene
    set by identifier, and identifiers disagree between files far more often
    than they agree (gene vs transcript accession, an Ensembl 'gene:'
    prefix, a trailing '.N'). Silence here is the worst outcome, so every
    join reports what matched, what did not, and enough example ids from
    BOTH sides to tell an id-format problem apart from genuinely absent
    genes.
    """
    n, total = len(matched), len(genes)
    pct = 100.0 * n / total if total else 0.0
    head = f"[build_synteny_html] {kind} {label}: {n}/{total} genes matched ({pct:.1f}%)"
    if table:
        head += f", {len(table)} rows in file"
    MATCH_LOG.append(head)
    sys.stderr.write(head + "\n")

    if not table:
        msg = f"  no {kind} file supplied for {label} - all values will be NA"
        MATCH_LOG.append(msg)
        return

    if n == 0 and total:
        msg = (f"  WARNING: NOTHING matched. Example gene id from the GFF: "
               f"{genes[0]['gene_id']!r}; example id in the {kind} file: "
               f"{next(iter(table))!r}. These look like different id "
               f"conventions - check the first column of your {kind} file.")
        MATCH_LOG.append(msg)
        sys.stderr.write(msg + "\n")
    elif missed:
        shown = ", ".join(repr(m) for m in missed[:10])
        more = f" (+{len(missed)-10} more)" if len(missed) > 10 else ""
        msg = (f"  {len(missed)} gene(s) from the GFF had no {kind} row, "
               f"shown as NA: {shown}{more}")
        MATCH_LOG.append(msg)
        sys.stderr.write(msg + "\n")

    # ids present in the user's file that never matched any displayed gene -
    # usually just genes outside this QTL, but a very high number alongside
    # a low match rate points at an id-format mismatch instead
    gene_keys = set()
    for g in genes:
        gene_keys.add(g["gene_id"])
        gene_keys.add(norm_id(g["gene_id"]))
    unused = [k for k in table if k not in gene_keys and norm_id(k) not in gene_keys]
    if unused and n < total:
        msg = (f"  note: {len(unused)}/{len(table)} rows in the {kind} file "
               f"matched no displayed gene (expected - most will be genes "
               f"outside this QTL)")
        MATCH_LOG.append(msg)


def match_expression(expr, genes, label):
    """Join the counts table to the displayed genes."""
    if not expr:
        report_match("expression", label, genes, expr, {}, [g["gene_id"] for g in genes])
        return {}

    index = {}
    for k, v in expr.items():
        index.setdefault(k, v)
        index.setdefault(norm_id(k), v)

    out, misses = {}, []
    for g in genes:
        gid = g["gene_id"]
        hit = index.get(gid)
        if hit is None:
            hit = index.get(norm_id(gid))
        if hit is not None:
            out[gid] = hit
        else:
            misses.append(gid)

    report_match("expression", label, genes, expr, out, misses)
    return out


def read_annotation(path):
    """Generic gene-annotation TSV: first column is the gene id, every other
    column is free-form and passed through as-is, however many there are.

    Returns (columns, {gene_id: [values]}). No assumption about what the
    extra columns mean - whatever the header says is what gets displayed.
    """
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


def match_annotation(columns, table, genes, label):
    """Same tolerant join as match_expression: exact id, then normalised."""
    if not table:
        report_match("annotation", label, genes, table, {},
                     [g["gene_id"] for g in genes])
        return {}
    index = {}
    for k, v in table.items():
        index.setdefault(k, v)
        index.setdefault(norm_id(k), v)

    out, misses = {}, []
    for g in genes:
        gid = g["gene_id"]
        hit = index.get(gid) or index.get(norm_id(gid))
        if hit is not None:
            out[gid] = hit
        else:
            misses.append(gid)

    report_match("annotation", label, genes, table, out, misses)
    if columns:
        MATCH_LOG.append(f"  {len(columns)} annotation column(s): "
                         f"{', '.join(columns)}")
    return out


def _unoffset(name):
    """`seqid:start-end` (1-based, from extract_region_fasta.py) -> (seqid, offset).

    The region FASTAs are extracted subsequences, so PAF coordinates are local
    to the slice; add the offset back to get true genome coordinates. Plain
    names (no suffix) are returned with offset 0.
    """
    if ":" in name and "-" in name.rsplit(":", 1)[-1]:
        base, span = name.rsplit(":", 1)
        try:
            start = int(span.split("-")[0])
            return base, start - 1
        except ValueError:
            pass
    return name, 0


def read_paf(path, tgt_seqids):
    """Alignment blocks source(query) -> target, restricted to displayed
    target sequences, with region-slice offsets resolved."""
    blocks = []
    if not path or not os.path.exists(path):
        return blocks
    keep = set(tgt_seqids)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, qoff = _unoffset(f[0])
            tname, toff = _unoffset(f[5])
            if tname not in keep:
                continue
            try:
                blocks.append(dict(
                    q=qname, qs=int(f[2]) + qoff, qe=int(f[3]) + qoff, strand=f[4],
                    t=tname, ts=int(f[7]) + toff, te=int(f[8]) + toff,
                    nmatch=int(f[9]), alnlen=int(f[10]), mapq=int(f[11])))
            except ValueError:
                continue
    return blocks


# --------------------------------------------------------------------------
HTML = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>__TITLE__</title>
<style>
:root{
  --bg:#ffffff; --panel:#f7f8fa; --line:#d8dce3; --ink:#1c2027; --muted:#697386;
  --src:#2f6f9f; --tgt:#3f8f6f; --rbh:#c2571a; --mp:#7a5bb5; --sel:#d92b2b;
}
*{box-sizing:border-box}
body{margin:0;font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
     color:var(--ink);background:var(--bg)}
header{padding:14px 18px;border-bottom:1px solid var(--line);background:var(--panel)}
h1{margin:0 0 4px;font-size:16px;font-weight:650}
.sub{color:var(--muted);font-size:12px}
.wrap{padding:14px 18px;max-width:1600px}
.panel{border:1px solid var(--line);border-radius:6px;margin-bottom:16px;background:#fff}
.panel h2{margin:0;padding:9px 12px;font-size:13px;font-weight:600;
          border-bottom:1px solid var(--line);background:var(--panel);
          border-radius:6px 6px 0 0}
.panel .body{padding:10px 12px;overflow-x:auto}
.controls{display:flex;gap:14px;align-items:center;flex-wrap:wrap;
          padding:8px 12px;border-bottom:1px solid var(--line);font-size:12px}
.controls label{display:flex;gap:5px;align-items:center;color:var(--muted)}
button{font:inherit;padding:3px 9px;border:1px solid var(--line);background:#fff;
       border-radius:4px;cursor:pointer}
button:hover{background:var(--panel)}
button.on{background:var(--ink);color:#fff;border-color:var(--ink)}
table{border-collapse:collapse;font-size:12px;width:100%}
th,td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;
      white-space:nowrap}
th{background:var(--panel);position:sticky;top:0;font-weight:600;cursor:pointer}
td.num{text-align:right;font-variant-numeric:tabular-nums}
tr.sel{background:#fff2f2}
.tblwrap{max-height:420px;overflow:auto;border:1px solid var(--line);border-radius:4px}
.legend{display:flex;gap:16px;flex-wrap:wrap;color:var(--muted);font-size:12px;
        padding:6px 12px}
.key{display:inline-block;width:11px;height:11px;border-radius:2px;
     vertical-align:-1px;margin-right:4px}
#tip{position:fixed;pointer-events:none;background:#1c2027;color:#fff;padding:6px 8px;
     border-radius:4px;font-size:11px;opacity:0;transition:opacity .08s;z-index:99;
     max-width:320px}
svg{display:block}
.gene{cursor:pointer}
.gene:hover{opacity:.75}
text{font:11px -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
.axis{fill:var(--muted);font-size:10px}
.empty{color:var(--muted);padding:14px 12px;font-style:italic}
.badges{display:flex;gap:10px;margin-top:10px}
.badge{flex:1;border:1px solid var(--line);border-radius:6px;padding:8px 12px;background:#fff}
.badge .num{font-size:20px;font-weight:700;line-height:1.1}
.badge .lbl{font-size:11px;color:var(--muted);margin-top:2px}
.badge.with .num{color:var(--src)}
.badge.without .num{color:#8a6d1a}
.badge.total .num{color:var(--ink)}
</style></head><body>
<header>
  <h1>__TITLE__</h1>
  <div class="sub">__SUBTITLE__</div>
  <div class="badges" id="covBadges"></div>
</header>
<div class="wrap">

  <div class="panel">
    <h2>Synteny</h2>
    <div class="controls">
      <label><input type="checkbox" id="cbAln" checked> DNA synteny (minimap2)</label>
      <label><input type="checkbox" id="cbMiniprot" checked> miniprot placements</label>
      <label><input type="checkbox" id="cbRbh" checked> reciprocal best hits</label>
      <label>Min % identity <input type="range" id="rgIdent" min="0" max="100" value="0" style="width:110px"><span id="lbIdent">0</span></label>
      <button id="btnReset" class="btnClearSel">Clear selection</button>
      <button class="btnDisplaySel">Display selected in alignment</button>
      <button id="btnClearLane" style="display:none">Show all regions</button>
      <button id="btnSynSvg">Download SVG</button>
      <button id="btnSynCsv">Download regions CSV</button>
      <span id="laneHint" class="sub" style="display:none">click a target label to isolate it</span>
    </div>
    <div class="legend">
      <span><i class="key" style="background:var(--src)"></i>source gene</span>
      <span><i class="key" style="background:var(--tgt)"></i>target gene</span>
      <span><i class="key" style="background:#b9c3d0"></i>DNA block (defines the region)</span>
      <span><i class="key" style="background:var(--mp)"></i>miniprot placement</span>
      <span><i class="key" style="background:var(--rbh)"></i>reciprocal best hit</span>
      <span>gene arrows show strand &middot; click a gene (here, in the expression table, or in annotation) to select it &middot; click a target label to isolate that region</span>
    </div>
    <div class="body"><svg id="syn"></svg></div>
  </div>

  <div class="panel" id="annPanel" style="display:none">
    <h2>Gene annotation</h2>
    <div class="controls">
      <span id="annNote" class="sub"></span>
      <button class="btnClearSel">Clear selection</button>
      <button class="btnDisplaySel">Display selected in alignment</button>
      <button id="btnAnnCsv">Download CSV</button>
    </div>
    <div class="body"><div class="tblwrap"><table id="annTbl"></table></div></div>
  </div>

  <div class="panel">
    <h2>Expression</h2>
    <div class="controls">
      <label>Scale
        <select id="selScale">
          <option value="log">log2(x+1)</option>
          <option value="raw">raw</option>
          <option value="zrow">z-score per gene</option>
        </select>
      </label>
      <label>Show
        <select id="selWhich">
          <option value="both">source + target</option>
          <option value="src">source only</option>
          <option value="tgt">target only</option>
        </select>
      </label>
      <label><input type="checkbox" id="cbLinkedOnly"> only genes with links</label>
      <label><input type="checkbox" id="cbTranspose"> transpose (genes on Y)</label>
    </div>
    <div class="body"><svg id="heat"></svg></div>
  </div>

  <div class="panel">
    <h2>Expression table</h2>
    <div class="controls">
      <span id="tblNote" class="sub">click a gene above to filter &mdash; showing all</span>
      <button class="btnClearSel">Clear selection</button>
      <button class="btnDisplaySel">Display selected in alignment</button>
      <button id="btnCsv">Download CSV</button>
    </div>
    <div class="body"><div class="tblwrap"><table id="tbl"></table></div></div>
  </div>

  <div class="panel">
    <h2>Syntenic regions &mdash; ranked by how much of the source QTL each covers</h2>
    <div class="body"><div class="tblwrap"><table id="regTbl"></table></div></div>
  </div>

</div>
<div id="tip"></div>
<script>
const DATA = __DATA__;
const tip = document.getElementById('tip');
let SEL = new Set();       // selected gene ids (multi-select)
function toggleGene(id){
  if (SEL.has(id)) SEL.delete(id); else SEL.add(id);
  redraw();
}
function selectedLabel(){
  const n = SEL.size;
  if (n === 0) return '';
  if (n === 1) return esc([...SEL][0]);
  return `${n} selected genes`;
}
let SHOW_SELECTED_ONLY = false;   // filter the idiogram + heatmap to selection
function toggleDisplaySelected(){
  SHOW_SELECTED_ONLY = !SHOW_SELECTED_ONLY;
  redraw();
}
function clearSelection(){
  // A full reset, not just an empty set: if Display was left on from a
  // previous round, the very next gene picked would be instantly filtered
  // to itself (SHOW_SELECTED_ONLY still true, SEL.size back to 1) -
  // reproducing the exact "only lets me pick one" symptom this button is
  // supposed to fix. Clearing always returns to a fresh, unfiltered state.
  SEL = new Set();
  SHOW_SELECTED_ONLY = false;
  redraw();
}
function relevantGenes(){
  // the current selection, plus every gene linked to any selected gene
  // (source<->target via RBH/miniprot) - this is what "Display selected"
  // filters down to, so a selected source gene brings its target match(es)
  // into view and vice versa.
  const rel = new Set(SEL);
  DATA.links.forEach(l => {
    if (SEL.has(l.src_gene)) rel.add(l.tgt_gene);
    if (SEL.has(l.tgt_gene)) rel.add(l.src_gene);
  });
  return rel;
}
function updateDisplayButtons(){
  document.querySelectorAll('.btnDisplaySel').forEach(b => {
    b.classList.toggle('on', SHOW_SELECTED_ONLY);
    b.textContent = SHOW_SELECTED_ONLY ? 'Showing selected only' : 'Display selected in alignment';
  });
}
let SEL_LANE = null;       // selected target region rank (isolates one lane)
const state = {aln:true, miniprot:true, rbh:true, minIdent:0,
               scale:'log', which:'both', linkedOnly:false, transpose:false};

function showTip(e, html){
  tip.innerHTML = html; tip.style.opacity = 1;
  const pad = 14;
  let x = e.clientX + pad, y = e.clientY + pad;
  if (x + 330 > window.innerWidth) x = e.clientX - 330;
  tip.style.left = x + 'px'; tip.style.top = y + 'px';
}
function hideTip(){ tip.style.opacity = 0; }
const esc = s => String(s).replace(/[&<>]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));

/* ---------------------------------------------------------------- synteny */
function drawSynteny(){
  const svg = document.getElementById('syn');

  // "Display selected in alignment": when on and something is selected,
  // only regions containing a relevant gene (selected, or linked to a
  // selected gene) get a lane at all - genuinely filtered, not just dimmed.
  const filtering = SHOW_SELECTED_ONLY && SEL.size > 0;
  const rel = filtering ? relevantGenes() : null;
  const visibleRegions = filtering
    ? DATA.regions.filter(r => DATA.target_genes.some(g => g.seqid === r.tgt_seqid && rel.has(g.gene_id)))
    : DATA.regions;

  const W = Math.max(1100, (visibleRegions.length ? 1100 : 900));
  const laneH = 54, gap = 76, padL = 150, padR = 24, top = 44;
  const trackW = W - padL - padR;
  const nLanes = 1 + visibleRegions.length;
  const H = top + nLanes * laneH + (nLanes - 1) * gap + 20;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('width', W); svg.setAttribute('height', H);
  let s = '';

  const src = DATA.source;
  const srcLen = Math.max(src.end - src.start, 1);
  const sx = p => padL + (Math.min(Math.max(p, src.start), src.end) - src.start) / srcLen * trackW;
  const srcY = top;

  // lane scales, one per VISIBLE region
  const lanes = visibleRegions.map((r, i) => {
    const len = Math.max(r.tgt_end - r.tgt_start, 1);
    return {r, y: top + (i + 1) * (laneH + gap),
            x: p => padL + (Math.min(Math.max(p, r.tgt_start), r.tgt_end) - r.tgt_start) / len * trackW};
  });
  const isDim = r => SEL_LANE != null && r.rank !== SEL_LANE;
  const DIM = 0.12;  // opacity multiplier for faded (non-selected) content

  function axis(y, label, a, b, xf){
    let o = `<line x1="${padL}" y1="${y+26}" x2="${padL+trackW}" y2="${y+26}" stroke="var(--line)"/>`;
    o += `<text x="8" y="${y+16}" font-weight="600">${esc(label)}</text>`;
    const n = 5;
    for (let i = 0; i <= n; i++){
      const p = a + (b - a) * i / n, x = xf(p);
      o += `<line x1="${x}" y1="${y+26}" x2="${x}" y2="${y+31}" stroke="var(--line)"/>`;
      o += `<text class="axis" x="${x}" y="${y+42}" text-anchor="middle">${(p/1e6).toFixed(2)}Mb</text>`;
    }
    return o;
  }

  // DNA alignment ribbons (drawn first, behind genes)
  if (state.aln){
    for (const lane of lanes){
      for (const b of DATA.alignments){
        if (b.t !== lane.r.tgt_seqid) continue;
        if (b.qe < src.start || b.qs > src.end) continue;
        if (b.te < lane.r.tgt_start || b.ts > lane.r.tgt_end) continue;
        const x1 = sx(b.qs), x2 = sx(b.qe);
        const x3 = lane.x(b.strand === '+' ? b.te : b.ts);
        const x4 = lane.x(b.strand === '+' ? b.ts : b.te);
        const yA = srcY + 26, yB = lane.y + 8;
        const inv = b.strand === '-';
        const dim = isDim(lane.r);
        const fop = (dim ? DIM : 1) * 0.5;
        s += `<path d="M${x1},${yA} L${x2},${yA} L${x3},${yB} L${x4},${yB} Z"
               fill="${dim ? '#c7ccd3' : (inv ? '#e8d5c4' : '#b9c3d0')}" fill-opacity="${fop}" stroke="none"
               ${dim ? '' : `data-tip="DNA block &middot; ${b.alnlen.toLocaleString()} bp &middot; ${(100*b.nmatch/b.alnlen).toFixed(1)}% id &middot; strand ${b.strand} &middot; MAPQ ${b.mapq}"`}></path>`;
      }
    }
  }

  // protein evidence links - miniprot and RBH are independent tracks
  {
    for (const lane of lanes){
      for (const l of DATA.links){
        if (l.tgt_seqid !== lane.r.tgt_seqid) continue;
        if (l.track === 'miniprot' && !state.miniprot) continue;
        if (l.track === 'rbh' && !state.rbh) continue;
        if (l.pident < state.minIdent) continue;
        const x1 = sx((l.src_start + l.src_end)/2);
        const x2 = lane.x((l.tgt_start + l.tgt_end)/2);
        const yA = srcY + 26, yB = lane.y + 8;
        const my = (yA + yB)/2;
        const dim = isDim(lane.r);
        const isSel = !dim && (SEL.has(l.src_gene) || SEL.has(l.tgt_gene));
        const isRbh = l.track === 'rbh';
        const col = dim ? '#c7ccd3' : (isSel ? 'var(--sel)' : (isRbh ? 'var(--rbh)' : 'var(--mp)'));
        const w = isSel ? 2.2 : (isRbh ? 1.6 : 1.0);
        const op = (dim ? DIM : 1) * (isSel ? 1 : (isRbh ? .9 : .55));
        const dash = isRbh ? '' : ' stroke-dasharray="3,2"';
        const partner = l.tgt_gene ? esc(l.tgt_gene) : '(unannotated locus)';
        s += `<path d="M${x1},${yA} C${x1},${my} ${x2},${my} ${x2},${yB}"
               fill="none" stroke="${col}" stroke-width="${w}" stroke-opacity="${op}"${dash}
               ${dim ? '' : `data-tip="<b>${isRbh ? 'RBH' : 'miniprot'}</b><br>${esc(l.src_gene)} &rarr; ${partner}<br>${l.pident.toFixed(1)}% id${l.bits ? ' &middot; bits ' + l.bits : ''}"`}></path>`;
      }
    }
  }

  function geneArrow(g, x1, x2, y, fill, dim){
    const h = 11, tipw = Math.min(6, Math.max(2, x2 - x1));
    const sel = !dim && SEL.has(g.gene_id);
    const f = dim ? '#c7ccd3' : (sel ? 'var(--sel)' : fill);
    let d;
    if (g.strand === '-')
      d = `M${x2},${y} L${x1+tipw},${y} L${x1},${y+h/2} L${x1+tipw},${y+h} L${x2},${y+h} Z`;
    else
      d = `M${x1},${y} L${x2-tipw},${y} L${x2},${y+h/2} L${x2-tipw},${y+h} L${x1},${y+h} Z`;
    const clickable = dim ? '' : `class="gene" data-gene="${esc(g.gene_id)}"`;
    const tip = dim ? '' : `data-tip="<b>${esc(g.gene_id)}</b><br>${esc(g.seqid)}:${g.start.toLocaleString()}-${g.end.toLocaleString()} (${g.strand})"`;
    return `<path ${clickable} d="${d}" fill="${f}" fill-opacity="${dim ? DIM : 1}"
             stroke="${sel ? 'var(--sel)' : 'none'}" stroke-width="${sel ? 1.5 : 0}" ${tip}></path>`;
  }

  // coverage strip: union of every displayed target's blocks, so gaps in
  // the QTL that nothing explains are visible at a glance
  {
    const iv = DATA.alignments.map(b => [Math.min(b.qs,b.qe), Math.max(b.qs,b.qe)])
                              .sort((a,b) => a[0]-b[0]);
    const merged = [];
    for (const [a,b] of iv){
      if (merged.length && a <= merged[merged.length-1][1])
        merged[merged.length-1][1] = Math.max(merged[merged.length-1][1], b);
      else merged.push([a,b]);
    }
    const cy = srcY - 14;
    s += `<rect x="${padL}" y="${cy}" width="${trackW}" height="7" fill="#eceff3" rx="2"></rect>`;
    for (const [a,b] of merged){
      const x1 = sx(a), x2 = Math.max(sx(b), sx(a)+1);
      s += `<rect x="${x1}" y="${cy}" width="${x2-x1}" height="7" fill="var(--src)" fill-opacity=".55"
             data-tip="covered ${a.toLocaleString()}-${b.toLocaleString()}"></rect>`;
    }
    const cov = DATA.coverage || {};
    s += `<text x="8" y="${cy+7}" class="axis">gene-evidence coverage ${(cov.with_genes_pct||0).toFixed(1)}%</text>`;
  }

  // source lane
  s += axis(srcY, 'SOURCE', src.start, src.end, sx);
  s += `<text x="8" y="${srcY+30}" class="axis">${esc(src.chrom)}</text>`;
  const selRegion = SEL_LANE == null ? null : DATA.regions.find(r => r.rank === SEL_LANE);
  const linkedSrcGenes = selRegion
    ? new Set(DATA.links.filter(l => l.tgt_seqid === selRegion.tgt_seqid).map(l => l.src_gene))
    : null;
  for (const g of DATA.source_genes){
    if (filtering && !rel.has(g.gene_id)) continue;
    const x1 = sx(g.start), x2 = Math.max(sx(g.end), sx(g.start) + 2.5);
    const dimGene = selRegion != null && !linkedSrcGenes.has(g.gene_id);
    s += geneArrow(g, x1, x2, srcY + 10, 'var(--src)', dimGene);
  }

  // target lanes
  for (const lane of lanes){
    const r = lane.r;
    const dim = isDim(r);
    const active = SEL_LANE === r.rank;
    const lop = dim ? DIM : 1;
    if (active){
      s += `<rect x="4" y="${lane.y-38}" width="${padL+trackW-8}" height="${laneH+34}"
             fill="var(--sel)" fill-opacity=".06" stroke="var(--sel)" stroke-opacity=".4"
             stroke-width="1" rx="4"></rect>`;
    }
    s += axis(lane.y - 18, `TARGET #${r.rank}`, r.tgt_start, r.tgt_end, lane.x);
    s += `<g class="laneLabel" data-lane="${r.rank}" opacity="${lop}" style="cursor:pointer">
      <rect x="4" y="${lane.y-14}" width="140" height="70" fill="transparent"></rect>
      <text x="8" y="${lane.y+2}" class="axis" font-weight="${active?700:400}">${esc(r.tgt_seqid)} &#128269;</text>
      <text x="8" y="${lane.y+16}" class="axis" font-weight="600">covers ${r.src_cov_pct.toFixed(1)}% of QTL</text>
      <text x="8" y="${lane.y+28}" class="axis">${r.pct_id.toFixed(0)}% id &middot; ${(r.aligned_bp/1000).toFixed(0)}kb aligned</text>
      <text x="8" y="${lane.y+40}" class="axis">${r.n_tgt_genes} genes &middot; ${r.n_miniprot} mp &middot; ${r.n_rbh} rbh</text>
      ${r.group_size > 1 ? `<text x="8" y="${lane.y+52}" class="axis" fill="${r.is_primary_in_group ? '#1a7a4c' : '#8a6d1a'}" font-weight="700"
             >${r.is_primary_in_group ? `&#9733; best match (group of ${r.group_size})` : `&#8635; possible copy of ${esc(r.group_best_scaffold)}`}</text>` : ''}
      ${(r.homoeolog_call && (r.homoeolog_call.startsWith('HOMOEOLOG') || r.homoeolog_call.includes('CONFLICT'))) ?
        `<text x="8" y="${lane.y+52}" class="axis" fill="var(--rbh)" font-weight="700">&#9888; ${esc(r.homoeolog_call)}</text>` : ''}
    </g>`;
    for (const g of DATA.target_genes){
      if (g.seqid !== r.tgt_seqid) continue;
      if (g.end < r.tgt_start || g.start > r.tgt_end) continue;
      if (filtering && !rel.has(g.gene_id)) continue;
      const x1 = lane.x(g.start), x2 = Math.max(lane.x(g.end), lane.x(g.start) + 2.5);
      s += geneArrow(g, x1, x2, lane.y + 8, 'var(--tgt)', isDim(r));
    }
  }

  if (!visibleRegions.length && filtering)
    s += `<text x="${padL}" y="${top+80}" class="axis">No target region contains the selected gene(s) or their linked partners.</text>`;
  else if (!visibleRegions.length)
    s += `<text x="${padL}" y="${top+80}" class="axis">No candidate target regions passed the thresholds for this QTL.</text>`;

  svg.innerHTML = s;
  svg.querySelectorAll('[data-tip]').forEach(el => {
    el.addEventListener('mousemove', e => showTip(e, el.dataset.tip));
    el.addEventListener('mouseleave', hideTip);
  });
  svg.querySelectorAll('[data-gene]').forEach(el => {
    el.addEventListener('click', (e) => {
      e.stopPropagation();
      toggleGene(el.dataset.gene);
    });
  });
  svg.querySelectorAll('.laneLabel').forEach(el => {
    el.addEventListener('click', () => {
      const rank = +el.dataset.lane;
      SEL_LANE = (SEL_LANE === rank) ? null : rank;
      redraw();
    });
  });
  updateLaneControls();
}

function updateLaneControls(){
  const btn = document.getElementById('btnClearLane');
  const hint = document.getElementById('laneHint');
  if (SEL_LANE != null){
    btn.style.display = '';
    btn.textContent = `Show all regions (isolating #${SEL_LANE})`;
    hint.style.display = 'none';
  } else {
    btn.style.display = 'none';
    hint.style.display = DATA.regions.length > 1 ? '' : 'none';
  }
}

/* -------------------------------------------------------------- heatmap */
function laneFilter(){
  // when a region is isolated, restrict to: that region's target genes, and
  // source genes linked into it. Returns null when nothing is isolated.
  if (SEL_LANE == null) return null;
  const r = DATA.regions.find(x => x.rank === SEL_LANE);
  if (!r) return null;
  const tgtIds = new Set(DATA.target_genes.filter(g => g.seqid === r.tgt_seqid).map(g => g.gene_id));
  const srcIds = new Set(DATA.links.filter(l => l.tgt_seqid === r.tgt_seqid).map(l => l.src_gene));
  return {tgtIds, srcIds, seqid: r.tgt_seqid};
}

function heatGenes(){
  let gs = [];
  const linked = new Set();
  DATA.links.forEach(l => { linked.add(l.src_gene); linked.add(l.tgt_gene); });
  if (state.which !== 'tgt') gs = gs.concat(DATA.source_genes.map(g => ({...g, side:'src'})));
  if (state.which !== 'src') gs = gs.concat(DATA.target_genes.map(g => ({...g, side:'tgt'})));
  if (state.linkedOnly) gs = gs.filter(g => linked.has(g.gene_id));
  const lf = laneFilter();
  if (lf) gs = gs.filter(g => g.side === 'src' ? lf.srcIds.has(g.gene_id) : lf.tgtIds.has(g.gene_id));
  if (SHOW_SELECTED_ONLY && SEL.size){
    const rel = relevantGenes();
    gs = gs.filter(g => rel.has(g.gene_id));
  }
  // Genes with no expression row are KEPT - their cells render in the
  // neutral "no data" colour. Dropping them here would make the heatmap's
  // gene set disagree with the tables and the synteny panel, which is
  // exactly the mismatch this pipeline keeps having to fix.
  return gs;
}

function scaleVals(vals, mode){
  if (mode === 'raw') return vals.slice();
  if (mode === 'log') return vals.map(v => v == null ? null : Math.log2(v + 1));
  const nums = vals.filter(v => v != null);
  if (!nums.length) return vals.slice();
  const m = nums.reduce((a,b)=>a+b,0)/nums.length;
  const sd = Math.sqrt(nums.reduce((a,b)=>a+(b-m)*(b-m),0)/nums.length) || 1;
  return vals.map(v => v == null ? null : (v - m)/sd);
}

function colour(v, lo, hi, diverging){
  if (v == null) return '#eef0f3';
  const t = hi > lo ? (v - lo)/(hi - lo) : 0.5;
  if (diverging){
    const u = Math.max(0, Math.min(1, t));
    return u < .5
      ? `rgb(${Math.round(49+(247-49)*u*2)},${Math.round(84+(247-84)*u*2)},${Math.round(141+(247-141)*u*2)})`
      : `rgb(${Math.round(247-(247-178)*(u-.5)*2)},${Math.round(247-(247-24)*(u-.5)*2)},${Math.round(247-(247-43)*(u-.5)*2)})`;
  }
  const u = Math.max(0, Math.min(1, t));
  return `rgb(${Math.round(255-(255-8)*u)},${Math.round(255-(255-64)*u)},${Math.round(255-(255-129)*u)})`;
}

function drawHeat(){
  const svg = document.getElementById('heat');
  const genes = heatGenes();
  const samples = DATA.samples_src.length >= DATA.samples_tgt.length ? DATA.samples_src : DATA.samples_tgt;
  if (!genes.length || !samples.length){
    svg.setAttribute('viewBox','0 0 900 60'); svg.setAttribute('height',60);
    svg.innerHTML = `<text x="10" y="34" class="axis">No expression data for the genes in view. `
      + `Check the log: gene ids in the counts table must match the annotation.</text>`;
    return;
  }

  const mats = genes.map(g => {
    const src = (g.side === 'src' ? DATA.expr_src : DATA.expr_tgt)[g.gene_id] || [];
    return scaleVals(src, state.scale);
  });
  let lo = Infinity, hi = -Infinity;
  mats.forEach(r => r.forEach(v => { if (v != null){ if (v<lo) lo=v; if (v>hi) hi=v; } }));
  if (!isFinite(lo)){ lo = 0; hi = 1; }
  const div = state.scale === 'zrow';
  if (div){ const m = Math.max(Math.abs(lo), Math.abs(hi)); lo = -m; hi = m; }

  const T = state.transpose;
  // rows/cols of the drawn grid
  const nCol = T ? samples.length : genes.length;
  const nRow = T ? genes.length : samples.length;
  const cw = T ? Math.max(24, Math.min(46, Math.floor(900/Math.max(nCol,1))))
               : Math.max(7, Math.min(20, Math.floor(1200/Math.max(nCol,1))));
  const ch = T ? 13 : 17;
  const padL = T ? 230 : 150, padT = T ? 130 : 92, padB = 18, padR = 90;
  const W = padL + nCol*cw + padR, H = padT + nRow*ch + padB;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('width', W); svg.setAttribute('height', H);

  const val = (gi, si) => mats[gi][si];
  let s = '';

  // column headers
  for (let c = 0; c < nCol; c++){
    const x = padL + c*cw;
    const isGene = !T;
    const g = isGene ? genes[c] : null;
    const label = isGene ? g.gene_id : samples[c];
    const sel = isGene && SEL.has(g.gene_id);
    const col = isGene ? (sel ? 'var(--sel)' : (g.side==='src'?'var(--src)':'var(--tgt)')) : 'var(--muted)';
    s += `<text class="axis" transform="translate(${x+cw/2},${padT-6}) rotate(-58)"
           text-anchor="start" fill="${col}" font-weight="${sel?700:400}">${esc(label)}</text>`;
  }
  // row labels
  for (let r = 0; r < nRow; r++){
    const y = padT + r*ch + ch/2 + 3;
    const isGene = T;
    const g = isGene ? genes[r] : null;
    const label = isGene ? g.gene_id : samples[r];
    const sel = isGene && SEL.has(g.gene_id);
    const col = isGene ? (sel ? 'var(--sel)' : (g.side==='src'?'var(--src)':'var(--tgt)')) : 'var(--muted)';
    s += `<text class="axis" x="${padL-7}" y="${y}" text-anchor="end" fill="${col}"
           font-weight="${sel?700:400}">${esc(label)}</text>`;
  }
  // cells
  for (let r = 0; r < nRow; r++){
    for (let c = 0; c < nCol; c++){
      const gi = T ? r : c, si = T ? c : r;
      const g = genes[gi], sm = samples[si], v = val(gi, si);
      s += `<rect class="gene" x="${padL + c*cw}" y="${padT + r*ch}" width="${cw-1}" height="${ch-1}"
             fill="${colour(v, lo, hi, div)}" data-gene="${esc(g.gene_id)}"
             data-tip="<b>${esc(g.gene_id)}</b><br>${esc(sm)}: ${v==null?'NA':v.toFixed(3)}"></rect>`;
    }
  }
  // selection outlines (one per selected gene present in this view)
  genes.forEach((g, idx) => {
    if (!SEL.has(g.gene_id)) return;
    if (T) s += `<rect x="${padL-1}" y="${padT + idx*ch - 1}" width="${nCol*cw+1}" height="${ch+1}"
                  fill="none" stroke="var(--sel)" stroke-width="1.6"></rect>`;
    else   s += `<rect x="${padL + idx*cw - 1}" y="${padT-2}" width="${cw+1}" height="${nRow*ch+3}"
                  fill="none" stroke="var(--sel)" stroke-width="1.6"></rect>`;
  });
  // colour key
  const kx = padL + nCol*cw + 22, kh = Math.min(120, nRow*ch);
  for (let i = 0; i < 40; i++){
    const v = lo + (hi-lo)*(1 - i/39);
    s += `<rect x="${kx}" y="${padT + i*(kh/40)}" width="12" height="${kh/40+0.6}" fill="${colour(v, lo, hi, div)}"></rect>`;
  }
  s += `<text class="axis" x="${kx+16}" y="${padT+8}">${hi.toFixed(1)}</text>`;
  s += `<text class="axis" x="${kx+16}" y="${padT+kh}">${lo.toFixed(1)}</text>`;

  svg.innerHTML = s;
  svg.querySelectorAll('[data-tip]').forEach(el => {
    el.addEventListener('mousemove', e => showTip(e, el.dataset.tip));
    el.addEventListener('mouseleave', hideTip);
  });
  svg.querySelectorAll('[data-gene]').forEach(el => {
    el.addEventListener('click', () => toggleGene(el.dataset.gene));
  });
}

/* ---------------------------------------------------------------- tables */
let sortCol = null, sortDir = 1;
function tableRows(){
  const linked = new Map();
  DATA.links.forEach(l => {
    if (!linked.has(l.src_gene)) linked.set(l.src_gene, []);
    linked.get(l.src_gene).push(l);
  });
  const rows = [];
  const add = (g, side, samples, expr) => {
    const v = expr[g.gene_id];
    const ls = side === 'src' ? (linked.get(g.gene_id) || []) : DATA.links.filter(l => l.tgt_gene === g.gene_id);
    // Every gene gets a row even with no expression data - a row of nulls
    // (rendered as NA), not a skipped one, so this table's gene count
    // matches genes.tsv / the final gene table / the annotation panel.
    rows.push({gene:g.gene_id, side, seqid:g.seqid, start:g.start, end:g.end,
               strand:g.strand,
               partner: ls.map(l => (side==='src' ? l.tgt_gene : l.src_gene) || '(unannot)').join(','),
               evidence: [...new Set(ls.map(l => l.track))].join('+'),
               vals: v || samples.map(() => null), samples,
               hasExpr: !!v});
  };
  const lf = laneFilter();
  const srcGenes = lf ? DATA.source_genes.filter(g => lf.srcIds.has(g.gene_id)) : DATA.source_genes;
  const tgtGenes = lf ? DATA.target_genes.filter(g => lf.tgtIds.has(g.gene_id)) : DATA.target_genes;
  srcGenes.forEach(g => add(g, 'src', DATA.samples_src, DATA.expr_src));
  tgtGenes.forEach(g => add(g, 'tgt', DATA.samples_tgt, DATA.expr_tgt));
  // Only hide non-selected rows once "Display selected" is on. Before that,
  // selecting a gene must just highlight it - filtering on every click would
  // hide the very rows you still need to click to build a multi-selection.
  return (SHOW_SELECTED_ONLY && SEL.size) ? rows.filter(r => SEL.has(r.gene) ||
      DATA.links.some(l => (SEL.has(l.src_gene) && l.tgt_gene===r.gene) ||
                           (SEL.has(l.tgt_gene) && l.src_gene===r.gene))) : rows;
}

function drawTable(){
  const rows = tableRows();
  const maxS = Math.max(DATA.samples_src.length, DATA.samples_tgt.length);
  const sampleCols = (DATA.samples_src.length >= DATA.samples_tgt.length ? DATA.samples_src : DATA.samples_tgt);
  const cols = ['gene','side','seqid','start','end','strand','partner','evidence'];
  let h = '<thead><tr>' + cols.map((c,i)=>`<th data-c="${i}">${c}</th>`).join('')
        + sampleCols.map((s,i)=>`<th data-c="${cols.length+i}">${esc(s)}</th>`).join('') + '</tr></thead><tbody>';
  if (sortCol != null){
    rows.sort((a,b)=>{
      const get = r => sortCol < cols.length ? r[cols[sortCol]] : (r.vals[sortCol-cols.length] ?? -Infinity);
      const x = get(a), y = get(b);
      return (typeof x === 'number' && typeof y === 'number') ? (x-y)*sortDir
           : String(x).localeCompare(String(y))*sortDir;
    });
  }
  for (const r of rows){
    h += `<tr class="${SEL.has(r.gene)?'sel':''}" data-gene="${esc(r.gene)}">`
       + `<td>${esc(r.gene)}</td><td>${r.side}</td><td>${esc(r.seqid)}</td>`
       + `<td class="num">${r.start.toLocaleString()}</td><td class="num">${r.end.toLocaleString()}</td>`
       + `<td>${r.strand}</td><td>${esc(r.partner)}</td><td>${esc(r.evidence)}</td>`;
    for (let i = 0; i < sampleCols.length; i++){
      const v = r.vals[i];
      h += `<td class="num">${v==null?'NA':v.toFixed(2)}</td>`;
    }
    h += '</tr>';
  }
  h += '</tbody>';
  const t = document.getElementById('tbl');
  t.innerHTML = h;
  t.querySelectorAll('th').forEach(th => th.addEventListener('click', ()=>{
    const c = +th.dataset.c;
    if (sortCol === c) sortDir = -sortDir; else { sortCol = c; sortDir = 1; }
    drawTable();
  }));
  t.querySelectorAll('tr[data-gene]').forEach(tr => tr.addEventListener('click', ()=>{
    toggleGene(tr.dataset.gene);
  }));
  const noteParts = [];
  if (SEL_LANE != null) noteParts.push(`isolated to region #${SEL_LANE}`);
  if (SHOW_SELECTED_ONLY && SEL.size)
    noteParts.push(`filtered to ${selectedLabel()} and linked partners`);
  else if (SEL.size)
    noteParts.push(`${SEL.size === 1 ? selectedLabel() + ' selected' : selectedLabel()} — click "Display selected" to filter`);
  document.getElementById('tblNote').textContent = noteParts.length
    ? `${noteParts.join(', ')} (${rows.length} rows)`
    : `click gene(s) to select — showing all ${rows.length} rows`;
}

function annotationRows(){
  const ann = DATA.annotation || {columns:[], src:{}, tgt:{}};
  if (!ann.columns.length) return {ann, rows:[]};
  const lf = laneFilter();
  const srcGenes = lf ? DATA.source_genes.filter(g => lf.srcIds.has(g.gene_id)) : DATA.source_genes;
  const tgtGenes = lf ? DATA.target_genes.filter(g => lf.tgtIds.has(g.gene_id)) : DATA.target_genes;
  let rows = [];
  const add = (g, side) => {
    const cols = side === 'src' ? ann.cols_src : ann.cols_tgt;
    const vals = (side === 'src' ? ann.src : ann.tgt)[g.gene_id];
    // Every gene gets a row, even with no annotation match - an "NA" row,
    // not a skipped one, so this table's gene count matches genes.tsv /
    // the final gene table exactly, the same way the expression table
    // already includes every gene regardless of expression data.
    rows.push({gene: g.gene_id, side, seqid: g.seqid,
              vals: vals || cols.map(() => 'NA'), cols});
  };
  srcGenes.forEach(g => add(g, 'src'));
  tgtGenes.forEach(g => add(g, 'tgt'));
  // Same rule as tableRows(): only filter once Display is explicitly on.
  if (SHOW_SELECTED_ONLY && SEL.size) rows = rows.filter(r => SEL.has(r.gene) ||
      DATA.links.some(l => (SEL.has(l.src_gene) && l.tgt_gene===r.gene) ||
                           (SEL.has(l.tgt_gene) && l.src_gene===r.gene)));
  return {ann, rows};
}

function drawAnnotationTable(){
  const panel = document.getElementById('annPanel');
  const {ann, rows} = annotationRows();
  if (!ann.columns.length){ panel.style.display = 'none'; return; }
  panel.style.display = '';

  const cols = ['gene', 'side', 'seqid', ...ann.columns];
  let h = '<thead><tr>' + cols.map(c => `<th>${esc(c)}</th>`).join('') + '</tr></thead><tbody>';
  for (const r of rows){
    h += `<tr class="${SEL.has(r.gene)?'sel':''}" data-gene="${esc(r.gene)}">`
       + `<td>${esc(r.gene)}</td><td>${r.side}</td><td>${esc(r.seqid)}</td>`;
    for (const c of ann.columns){
      const idx = r.cols.indexOf(c);
      h += `<td>${idx >= 0 ? esc(r.vals[idx]) : ''}</td>`;
    }
    h += '</tr>';
  }
  h += '</tbody>';
  const t = document.getElementById('annTbl');
  t.innerHTML = h;
  t.querySelectorAll('tr[data-gene]').forEach(tr => tr.addEventListener('click', () => {
    toggleGene(tr.dataset.gene);
  }));

  const noteParts = [];
  if (SEL_LANE != null) noteParts.push(`isolated to region #${SEL_LANE}`);
  if (SHOW_SELECTED_ONLY && SEL.size)
    noteParts.push(`filtered to ${selectedLabel()} and linked partners`);
  else if (SEL.size)
    noteParts.push(`${SEL.size === 1 ? selectedLabel() + ' selected' : selectedLabel()} — click "Display selected" to filter`);
  document.getElementById('annNote').textContent = noteParts.length
    ? `${noteParts.join(', ')} (${rows.length} rows)`
    : `${rows.length} annotated gene(s) in view`;
}

function downloadAnnotationCsv(){
  const {ann, rows} = annotationRows();
  const cols = ['gene', 'side', 'seqid', ...ann.columns];
  const lines = [cols.join(',')];
  for (const r of rows){
    const vals = ann.columns.map(c => {
      const idx = r.cols.indexOf(c);
      const v = idx >= 0 ? r.vals[idx] : '';
      return '"' + String(v).replace(/"/g, '""') + '"';
    });
    lines.push([r.gene, r.side, r.seqid].concat(vals).join(','));
  }
  const b = new Blob([lines.join('\\n')], {type:'text/csv'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(b);
  a.download = DATA.qtl_id.replace(/[^A-Za-z0-9._-]/g,'_') + '_annotation.csv';
  a.click();
}

function drawRegions(){
  const t = document.getElementById('regTbl');
  if (!DATA.regions.length){ t.innerHTML = '<tbody><tr><td class="empty">none</td></tr></tbody>'; return; }
  const cols = ['rank','tgt_seqid','src_cov_bp','src_cov_pct','aligned_bp','n_aln_blocks','pct_id','strand','tgt_cov_bp','n_tgt_genes','n_miniprot','n_rbh','region_group','group_size','is_primary_in_group','group_best_scaffold','consensus_src_chrom','n_chrom_rbh','n_chrom_miniprot','frac_on_qtl_chrom','homoeolog_call','tgt_start','tgt_end','tgt_span'];
  let h = '<thead><tr>' + cols.map(c=>`<th>${c}</th>`).join('') + '</tr></thead><tbody>';
  for (const r of DATA.regions){
    const active = SEL_LANE === r.rank;
    h += `<tr class="${active?'sel':''}" data-lane="${r.rank}" style="cursor:pointer">`
       + cols.map(c=>`<td class="${typeof r[c]==='number'?'num':''}">${esc(r[c])}</td>`).join('') + '</tr>';
  }
  t.innerHTML = h + '</tbody>';
  t.querySelectorAll('tr[data-lane]').forEach(tr => tr.addEventListener('click', () => {
    const rank = +tr.dataset.lane;
    SEL_LANE = (SEL_LANE === rank) ? null : rank;
    redraw();
  }));
}

function downloadCsv(){
  const rows = tableRows();
  const sampleCols = (DATA.samples_src.length >= DATA.samples_tgt.length ? DATA.samples_src : DATA.samples_tgt);
  const head = ['gene','side','seqid','start','end','strand','partner','evidence'].concat(sampleCols);
  const lines = [head.join(',')];
  for (const r of rows)
    lines.push([r.gene,r.side,r.seqid,r.start,r.end,r.strand,'"'+r.partner+'"',r.evidence]
      .concat(sampleCols.map((_,i)=> r.vals[i]==null?'NA':r.vals[i])).join(','));
  const b = new Blob([lines.join('\\n')], {type:'text/csv'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(b);
  a.download = DATA.qtl_id.replace(/[^A-Za-z0-9._-]/g,'_') + '_expression.csv';
  a.click();
}

function downloadSyntenySvg(){
  // Serialise the rendered SVG with the CSS variables resolved to literal
  // colours, otherwise the file opens colourless outside this page.
  const svg = document.getElementById('syn');
  const vars = {'--src':'#2f6f9f','--tgt':'#3f8f6f','--rbh':'#c2571a',
                '--mp':'#7a5bb5','--sel':'#d92b2b','--line':'#d8dce3',
                '--muted':'#697386','--ink':'#1c2027'};
  let inner = svg.innerHTML;
  for (const [k,v] of Object.entries(vars))
    inner = inner.split(`var(${k})`).join(v);
  const out = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${svg.getAttribute('viewBox')}" `
            + `width="${svg.getAttribute('width')}" height="${svg.getAttribute('height')}">`
            + `<style>text{font:11px -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}`
            + `.axis{fill:#697386;font-size:10px}</style>${inner}</svg>`;
  const b = new Blob([out], {type:'image/svg+xml'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(b);
  a.download = DATA.qtl_id.replace(/[^A-Za-z0-9._-]/g,'_') + '_synteny.svg';
  a.click();
}

function downloadRegionsCsv(){
  // exactly the regions currently drawn, in draw order
  const filtering = SHOW_SELECTED_ONLY && SEL.size > 0;
  const rel = filtering ? relevantGenes() : null;
  const regs = filtering
    ? DATA.regions.filter(r => DATA.target_genes.some(g => g.seqid === r.tgt_seqid && rel.has(g.gene_id)))
    : DATA.regions;
  const cols = ['rank','tgt_seqid','src_cov_bp','src_cov_pct','aligned_bp','n_aln_blocks',
                'pct_id','strand','tgt_cov_bp','n_tgt_genes','n_miniprot','n_rbh',
                'region_group','group_size','is_primary_in_group','group_best_scaffold',
                'consensus_src_chrom','n_chrom_rbh','n_chrom_miniprot','frac_on_qtl_chrom',
                'homoeolog_call','tgt_start','tgt_end','tgt_span'];
  const lines = [['qtl_id'].concat(cols).join(',')];
  for (const r of regs)
    lines.push([DATA.qtl_id].concat(cols.map(c => {
      const v = r[c];
      return typeof v === 'string' ? '"' + v.replace(/"/g,'""') + '"' : v;
    })).join(','));
  const b = new Blob([lines.join('\\n')], {type:'text/csv'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(b);
  a.download = DATA.qtl_id.replace(/[^A-Za-z0-9._-]/g,'_') + '_regions.csv';
  a.click();
}

function redraw(){
  drawSynteny(); drawHeat(); drawTable(); drawAnnotationTable(); drawRegions();
  updateDisplayButtons();
}

document.getElementById('cbAln').onchange = e => { state.aln = e.target.checked; drawSynteny(); };
document.getElementById('cbMiniprot').onchange = e => { state.miniprot = e.target.checked; drawSynteny(); };
document.getElementById('cbRbh').onchange = e => { state.rbh = e.target.checked; drawSynteny(); };
document.getElementById('rgIdent').oninput = e => {
  state.minIdent = +e.target.value;
  document.getElementById('lbIdent').textContent = e.target.value;
  drawSynteny();
};
document.getElementById('selScale').onchange = e => { state.scale = e.target.value; drawHeat(); };
document.getElementById('selWhich').onchange = e => { state.which = e.target.value; drawHeat(); };
document.getElementById('cbLinkedOnly').onchange = e => { state.linkedOnly = e.target.checked; drawHeat(); };
document.getElementById('cbTranspose').onchange = e => { state.transpose = e.target.checked; drawHeat(); };
document.getElementById('btnSynSvg').onclick = downloadSyntenySvg;
document.getElementById('btnSynCsv').onclick = downloadRegionsCsv;
document.querySelectorAll('.btnClearSel').forEach(b => b.onclick = clearSelection);
document.querySelectorAll('.btnDisplaySel').forEach(b => b.onclick = toggleDisplaySelected);
document.getElementById('btnClearLane').onclick = () => { SEL_LANE = null; redraw(); };
document.getElementById('btnCsv').onclick = downloadCsv;
document.getElementById('btnAnnCsv').onclick = downloadAnnotationCsv;

function drawCoverageBadges(){
  const c = DATA.coverage || {};
  const el = document.getElementById('covBadges');
  el.innerHTML = `
    <div class="badge with">
      <div class="num">${(c.with_genes_pct||0).toFixed(1)}%</div>
      <div class="lbl">covered by scaffolds WITH gene evidence (${(c.with_genes_bp||0).toLocaleString()} bp)</div>
    </div>
    <div class="badge without">
      <div class="num">${(c.without_genes_pct||0).toFixed(1)}%</div>
      <div class="lbl">covered ONLY by gene-less scaffolds (${(c.without_genes_bp||0).toLocaleString()} bp)</div>
    </div>
    <div class="badge total">
      <div class="num">${(c.total_pct||0).toFixed(1)}%</div>
      <div class="lbl">total DNA coverage, any evidence (${(c.total_bp||0).toLocaleString()} of ${(c.slice_len||0).toLocaleString()} bp)</div>
    </div>`;
}

drawCoverageBadges(); drawRegions(); drawAnnotationTable(); redraw();
</script></body></html>
"""


# --------------------------------------------------------------------------
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
    ap.add_argument("--qtl-chrom", required=True)
    ap.add_argument("--qtl-start", type=int, required=True)
    ap.add_argument("--qtl-end", type=int, required=True)
    ap.add_argument("--regions", required=True, help="annotate_regions.py --out-regions")
    ap.add_argument("--links", required=True, help="annotate_regions.py --out-links")
    ap.add_argument("--genes", required=True, help="annotate_regions.py --out-genes")
    ap.add_argument("--coverage", help="annotate_regions.py --out-coverage")
    ap.add_argument("--paf", help="source QTL slice vs target")
    ap.add_argument("--source-expr")
    ap.add_argument("--target-expr")
    ap.add_argument("--source-annotation",
                    help="TSV: gene id in col 1, any number of annotation "
                         "columns after, with a header row")
    ap.add_argument("--target-annotation")
    ap.add_argument("--out", required=True)
    ap.add_argument("--match-report",
                    help="write the gene/expression/annotation match "
                         "diagnostics to this file as well as stderr")
    ap.add_argument("--flank", type=int, default=0)
    ap.add_argument("--region-display-mode", choices=["evidence", "gene_length", "all"],
                    default="evidence",
                    help="'evidence' (default): only scaffolds with >=1 RBH "
                         "or miniprot hit. 'gene_length': also include "
                         "scaffolds with an annotated gene whose aligned "
                         "length clears --min-gene-region-bp, even with no "
                         "protein evidence. 'all': no filtering.")
    ap.add_argument("--min-gene-region-bp", type=int, default=5000,
                    help="gene_length mode only: minimum aligned bp for a "
                         "gene-bearing, evidence-less scaffold to be shown")
    args = ap.parse_args()

    regions = [r for r in read_tsv(args.regions) if r["qtl_id"] == args.qtl_id]
    links   = [l for l in read_tsv(args.links)   if l["qtl_id"] == args.qtl_id]
    genes   = [g for g in read_tsv(args.genes)   if g["qtl_id"] == args.qtl_id]

    all_regions = []
    for r in regions:
        all_regions.append(dict(
            rank=int(r["rank"]), tgt_seqid=r["tgt_seqid"],
            aligned_bp=int(r["aligned_bp"]), n_aln_blocks=int(r["n_aln_blocks"]),
            pct_id=float(r["pct_id"]), strand=r["strand"],
            src_cov_bp=int(r.get("src_cov_bp", 0) or 0),
            src_cov_pct=float(r.get("src_cov_pct", 0) or 0),
            homoeolog_call=r.get("homoeolog_call", "NA"),
            consensus_src_chrom=r.get("consensus_src_chrom", "NA"),
            frac_on_qtl_chrom=float(r.get("frac_on_qtl_chrom", 0) or 0),
            n_chrom_rbh=int(r.get("n_chrom_rbh", 0) or 0),
            n_chrom_miniprot=int(r.get("n_chrom_miniprot", 0) or 0),
            tgt_cov_bp=int(r.get("tgt_cov_bp", 0) or 0),
            n_tgt_genes=int(r.get("n_tgt_genes", 0) or 0),
            n_miniprot=int(r.get("n_miniprot", 0) or 0),
            n_rbh=int(r.get("n_rbh", 0) or 0),
            region_group=int(r.get("region_group", 0) or 0),
            group_size=int(r.get("group_size", 1) or 1),
            is_primary_in_group=int(r.get("is_primary_in_group", 1) or 0),
            group_best_scaffold=r.get("group_best_scaffold", ""),
            tgt_start=int(r["tgt_start"]), tgt_end=int(r["tgt_end"]),
            tgt_span=int(r["tgt_span"])))
    all_regions.sort(key=lambda x: x["rank"])

    # The viewer is for INSPECTING evidence: a scaffold with zero RBH and zero
    # miniprot hits has nothing to look at (it may still be real DNA homology,
    # but there is no protein-level signal to examine here), so it is hidden
    # by default. It still exists in --regions / the FASTA output upstream -
    # this filter only affects what gets drawn.
    reg_out = [r for r in all_regions
              if region_passes(r, args.region_display_mode, args.min_gene_region_bp)]
    n_hidden_empty = len(all_regions) - len(reg_out)
    tgt_seqids = [r["tgt_seqid"] for r in reg_out]

    link_out = []
    shown = set(tgt_seqids)
    for l in links:
        if l["tgt_seqid"] not in shown:
            continue
        link_out.append(dict(
            track=l["track"], src_gene=l["src_gene"],
            src_start=int(l["src_start"]), src_end=int(l["src_end"]),
            tgt_gene=l["tgt_gene"], tgt_seqid=l["tgt_seqid"],
            tgt_start=int(l["tgt_start"]), tgt_end=int(l["tgt_end"]),
            pident=float(l["pident"]), bits=int(float(l["bits"] or 0)),
            rbh=1 if l["track"] == "rbh" else 0))

    src_genes = [dict(gene_id=g["gene_id"], seqid=g["seqid"], start=int(g["start"]),
                      end=int(g["end"]), strand=g["strand"])
                 for g in genes if g["side"] == "source"]
    tgt_genes = [dict(gene_id=g["gene_id"], seqid=g["seqid"], start=int(g["start"]),
                      end=int(g["end"]), strand=g["strand"])
                 for g in genes if g["side"] == "target" and g["seqid"] in shown]

    view_start = max(args.qtl_start - args.flank, 0)
    view_end = args.qtl_end + args.flank

    samples_src, expr_src_all = read_expression(args.source_expr)
    samples_tgt, expr_tgt_all = read_expression(args.target_expr)
    expr_src = match_expression(expr_src_all, src_genes, "source")
    expr_tgt = match_expression(expr_tgt_all, tgt_genes, "target")

    cols_src, ann_src_all = read_annotation(args.source_annotation)
    cols_tgt, ann_tgt_all = read_annotation(args.target_annotation)
    ann_src = match_annotation(cols_src, ann_src_all, src_genes, "source")
    ann_tgt = match_annotation(cols_tgt, ann_tgt_all, tgt_genes, "target")
    # union of both sides' columns, source order first - so a table with
    # different columns per side (or only one side supplied) still renders
    # sensibly, with blanks where a column does not apply to that side
    ann_columns = list(cols_src) + [c for c in cols_tgt if c not in cols_src]

    alignments = read_paf(args.paf, tgt_seqids)

    cov_rows = read_tsv(args.coverage) if args.coverage else []
    cov_row = next((c for c in cov_rows if c["qtl_id"] == args.qtl_id), None)
    coverage = dict(
        with_genes_bp=int(cov_row["cov_with_genes_bp"]) if cov_row else 0,
        with_genes_pct=float(cov_row["cov_with_genes_pct"]) if cov_row else 0.0,
        without_genes_bp=int(cov_row["cov_without_genes_bp"]) if cov_row else 0,
        without_genes_pct=float(cov_row["cov_without_genes_pct"]) if cov_row else 0.0,
        total_bp=int(cov_row["cov_total_bp"]) if cov_row else 0,
        total_pct=float(cov_row["cov_total_pct"]) if cov_row else 0.0,
        slice_len=int(cov_row["src_slice_len"]) if cov_row else 0)

    data = dict(
        qtl_id=args.qtl_id,
        source=dict(chrom=args.qtl_chrom, start=view_start, end=view_end,
                    qtl_start=args.qtl_start, qtl_end=args.qtl_end),
        regions=reg_out, links=link_out,
        source_genes=src_genes, target_genes=tgt_genes,
        samples_src=samples_src, samples_tgt=samples_tgt,
        expr_src=expr_src, expr_tgt=expr_tgt,
        alignments=alignments, coverage=coverage,
        n_hidden_empty=n_hidden_empty, n_all_regions=len(all_regions),
        annotation=dict(columns=ann_columns, cols_src=cols_src, cols_tgt=cols_tgt,
                        src=ann_src, tgt=ann_tgt))

    n_mp = sum(1 for l in link_out if l["track"] == "miniprot")
    n_rb = sum(1 for l in link_out if l["track"] == "rbh")
    n_groups_multi = len({r["region_group"] for r in reg_out
                          if r["group_size"] > 1})
    hidden_note = (f" &middot; {n_hidden_empty} gene-less scaffold(s) hidden"
                  if n_hidden_empty else "")
    sub = (f"{args.qtl_chrom}:{args.qtl_start:,}-{args.qtl_end:,} "
           f"&middot; <b>{coverage['with_genes_pct']:.1f}% covered with gene evidence</b>, "
           f"{coverage['without_genes_pct']:.1f}% gene-less, "
           f"{coverage['total_pct']:.1f}% total "
           f"&middot; {len(reg_out)} scaffolds shown{hidden_note} "
           f"&middot; {n_groups_multi} possible-genome-copy group(s) "
           f"&middot; {len(alignments)} DNA blocks "
           f"&middot; {n_mp} miniprot, {n_rb} RBH links")

    html = (HTML.replace("__TITLE__", args.qtl_id)
                .replace("__SUBTITLE__", sub)
                .replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    with open(args.out, "w") as fh:
        fh.write(html)
    if args.match_report:
        with open(args.match_report, "w") as rf:
            rf.write(f"# gene matching report for QTL: {args.qtl_id}\n")
            rf.write(f"# {args.qtl_chrom}:{args.qtl_start}-{args.qtl_end}\n")
            rf.write(f"# genes in view: {len(src_genes)} source, "
                     f"{len(tgt_genes)} target\n#\n")
            rf.write("\n".join(MATCH_LOG) + "\n")

    sys.stderr.write(f"[build_synteny_html] {args.out} ({len(src_genes)} src, "
                     f"{len(tgt_genes)} tgt genes, {len(reg_out)}/{len(all_regions)} "
                     f"scaffolds shown, {n_mp} miniprot, {n_rb} rbh, "
                     f"{len(alignments)} aln blocks)\n")


if __name__ == "__main__":
    main()
