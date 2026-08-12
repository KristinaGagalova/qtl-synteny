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
import sys
from collections import defaultdict


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
    """First column = gene id, remaining columns = samples. Returns
    (samples, {gene: [values]}). Missing/þnon-numeric cells become None."""
    if not path or not os.path.exists(path):
        return [], {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        samples = hdr[1:]
        data = {}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            vals = []
            for x in f[1:len(samples) + 1]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(None)
            data[f[0]] = vals
    return samples, data


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
</style></head><body>
<header>
  <h1>__TITLE__</h1>
  <div class="sub">__SUBTITLE__</div>
</header>
<div class="wrap">

  <div class="panel">
    <h2>Synteny</h2>
    <div class="controls">
      <label><input type="checkbox" id="cbAln" checked> DNA synteny (minimap2)</label>
      <label><input type="checkbox" id="cbMiniprot" checked> miniprot placements</label>
      <label><input type="checkbox" id="cbRbh" checked> reciprocal best hits</label>
      <label>Min % identity <input type="range" id="rgIdent" min="0" max="100" value="0" style="width:110px"><span id="lbIdent">0</span></label>
      <button id="btnReset">Clear selection</button>
    </div>
    <div class="legend">
      <span><i class="key" style="background:var(--src)"></i>source gene</span>
      <span><i class="key" style="background:var(--tgt)"></i>target gene</span>
      <span><i class="key" style="background:#b9c3d0"></i>DNA block (defines the region)</span>
      <span><i class="key" style="background:var(--mp)"></i>miniprot placement</span>
      <span><i class="key" style="background:var(--rbh)"></i>reciprocal best hit</span>
      <span>gene arrows show strand &middot; click a gene to pin it</span>
    </div>
    <div class="body"><svg id="syn"></svg></div>
  </div>

  <div class="panel">
    <h2>Expression &mdash; genes on X, samples on Y</h2>
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
    </div>
    <div class="body"><svg id="heat"></svg></div>
  </div>

  <div class="panel">
    <h2>Expression table</h2>
    <div class="controls">
      <span id="tblNote" class="sub">click a gene above to filter &mdash; showing all</span>
      <button id="btnCsv">Download CSV</button>
    </div>
    <div class="body"><div class="tblwrap"><table id="tbl"></table></div></div>
  </div>

  <div class="panel">
    <h2>Syntenic regions (ranked by DNA alignment)</h2>
    <div class="body"><div class="tblwrap"><table id="regTbl"></table></div></div>
  </div>

</div>
<div id="tip"></div>
<script>
const DATA = __DATA__;
const tip = document.getElementById('tip');
let SEL = null;            // selected gene id
const state = {aln:true, miniprot:true, rbh:true, minIdent:0,
               scale:'log', which:'both', linkedOnly:false};

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
  const W = Math.max(1100, (DATA.regions.length ? 1100 : 900));
  const laneH = 54, gap = 46, padL = 118, padR = 24, top = 26;
  const trackW = W - padL - padR;
  const nLanes = 1 + DATA.regions.length;
  const H = top + nLanes * laneH + (nLanes - 1) * gap + 20;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('width', W); svg.setAttribute('height', H);
  let s = '';

  const src = DATA.source;
  const srcLen = Math.max(src.end - src.start, 1);
  const sx = p => padL + (Math.min(Math.max(p, src.start), src.end) - src.start) / srcLen * trackW;
  const srcY = top;

  // lane scales for each target region
  const lanes = DATA.regions.map((r, i) => {
    const len = Math.max(r.tgt_end - r.tgt_start, 1);
    return {r, y: top + (i + 1) * (laneH + gap),
            x: p => padL + (Math.min(Math.max(p, r.tgt_start), r.tgt_end) - r.tgt_start) / len * trackW};
  });

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
        s += `<path d="M${x1},${yA} L${x2},${yA} L${x3},${yB} L${x4},${yB} Z"
               fill="${inv ? '#e8d5c4' : '#b9c3d0'}" fill-opacity=".5" stroke="none"
               data-tip="DNA block &middot; ${b.alnlen.toLocaleString()} bp &middot; ${(100*b.nmatch/b.alnlen).toFixed(1)}% id &middot; strand ${b.strand} &middot; MAPQ ${b.mapq}"></path>`;
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
        const isSel = SEL && (l.src_gene === SEL || l.tgt_gene === SEL);
        const isRbh = l.track === 'rbh';
        const col = isSel ? 'var(--sel)' : (isRbh ? 'var(--rbh)' : 'var(--mp)');
        const w = isSel ? 2.2 : (isRbh ? 1.6 : 1.0);
        const op = isSel ? 1 : (isRbh ? .9 : .55);
        const dash = isRbh ? '' : ' stroke-dasharray="3,2"';
        const partner = l.tgt_gene ? esc(l.tgt_gene) : '(unannotated locus)';
        s += `<path d="M${x1},${yA} C${x1},${my} ${x2},${my} ${x2},${yB}"
               fill="none" stroke="${col}" stroke-width="${w}" stroke-opacity="${op}"${dash}
               data-tip="<b>${isRbh ? 'RBH' : 'miniprot'}</b><br>${esc(l.src_gene)} &rarr; ${partner}<br>${l.pident.toFixed(1)}% id${l.bits ? ' &middot; bits ' + l.bits : ''}"></path>`;
      }
    }
  }

  function geneArrow(g, x1, x2, y, fill){
    const h = 11, tipw = Math.min(6, Math.max(2, x2 - x1));
    const sel = SEL === g.gene_id;
    const f = sel ? 'var(--sel)' : fill;
    let d;
    if (g.strand === '-')
      d = `M${x2},${y} L${x1+tipw},${y} L${x1},${y+h/2} L${x1+tipw},${y+h} L${x2},${y+h} Z`;
    else
      d = `M${x1},${y} L${x2-tipw},${y} L${x2},${y+h/2} L${x2-tipw},${y+h} L${x1},${y+h} Z`;
    return `<path class="gene" d="${d}" fill="${f}" stroke="${sel ? 'var(--sel)' : 'none'}"
             stroke-width="${sel ? 1.5 : 0}" data-gene="${esc(g.gene_id)}"
             data-tip="<b>${esc(g.gene_id)}</b><br>${esc(g.seqid)}:${g.start.toLocaleString()}-${g.end.toLocaleString()} (${g.strand})"></path>`;
  }

  // source lane
  s += axis(srcY, 'SOURCE', src.start, src.end, sx);
  s += `<text x="8" y="${srcY+30}" class="axis">${esc(src.chrom)}</text>`;
  for (const g of DATA.source_genes){
    const x1 = sx(g.start), x2 = Math.max(sx(g.end), sx(g.start) + 2.5);
    s += geneArrow(g, x1, x2, srcY + 10, 'var(--src)');
  }

  // target lanes
  for (const lane of lanes){
    const r = lane.r;
    s += axis(lane.y - 18, `TARGET #${r.rank}`, r.tgt_start, r.tgt_end, lane.x);
    s += `<text x="8" y="${lane.y+2}" class="axis">${esc(r.tgt_seqid)}</text>`;
    s += `<text x="8" y="${lane.y+16}" class="axis">${(r.aligned_bp/1000).toFixed(0)}kb aln &middot; ${r.pct_id.toFixed(0)}%</text>`;
    s += `<text x="8" y="${lane.y+28}" class="axis">${r.n_tgt_genes} genes &middot; ${r.n_miniprot} mp &middot; ${r.n_rbh} rbh</text>`;
    for (const g of DATA.target_genes){
      if (g.seqid !== r.tgt_seqid) continue;
      if (g.end < r.tgt_start || g.start > r.tgt_end) continue;
      const x1 = lane.x(g.start), x2 = Math.max(lane.x(g.end), lane.x(g.start) + 2.5);
      s += geneArrow(g, x1, x2, lane.y + 8, 'var(--tgt)');
    }
  }

  if (!DATA.regions.length)
    s += `<text x="${padL}" y="${top+80}" class="axis">No candidate target regions passed the thresholds for this QTL.</text>`;

  svg.innerHTML = s;
  svg.querySelectorAll('[data-tip]').forEach(el => {
    el.addEventListener('mousemove', e => showTip(e, el.dataset.tip));
    el.addEventListener('mouseleave', hideTip);
  });
  svg.querySelectorAll('[data-gene]').forEach(el => {
    el.addEventListener('click', () => { SEL = (SEL === el.dataset.gene) ? null : el.dataset.gene; redraw(); });
  });
}

/* -------------------------------------------------------------- heatmap */
function heatGenes(){
  let gs = [];
  const linked = new Set();
  DATA.links.forEach(l => { linked.add(l.src_gene); linked.add(l.tgt_gene); });
  if (state.which !== 'tgt') gs = gs.concat(DATA.source_genes.map(g => ({...g, side:'src'})));
  if (state.which !== 'src') gs = gs.concat(DATA.target_genes.map(g => ({...g, side:'tgt'})));
  if (state.linkedOnly) gs = gs.filter(g => linked.has(g.gene_id));
  return gs.filter(g => (g.side === 'src' ? DATA.expr_src : DATA.expr_tgt)[g.gene_id]);
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
    svg.innerHTML = `<text x="10" y="34" class="axis">No expression data for the genes in view.</text>`;
    return;
  }
  const cw = Math.max(7, Math.min(20, Math.floor(1200/genes.length)));
  const ch = 17, padL = 150, padT = 92, padB = 18, padR = 90;
  const W = padL + genes.length*cw + padR, H = padT + samples.length*ch + padB;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('width', W); svg.setAttribute('height', H);

  const mats = genes.map(g => {
    const src = (g.side === 'src' ? DATA.expr_src : DATA.expr_tgt)[g.gene_id] || [];
    return scaleVals(src, state.scale);
  });
  let lo = Infinity, hi = -Infinity;
  mats.forEach(r => r.forEach(v => { if (v != null){ if (v<lo) lo=v; if (v>hi) hi=v; } }));
  if (!isFinite(lo)){ lo = 0; hi = 1; }
  const div = state.scale === 'zrow';
  if (div){ const m = Math.max(Math.abs(lo), Math.abs(hi)); lo = -m; hi = m; }

  let s = '';
  genes.forEach((g, i) => {
    const x = padL + i*cw;
    const sel = SEL === g.gene_id;
    s += `<text class="axis" transform="translate(${x+cw/2},${padT-6}) rotate(-58)"
           text-anchor="start" fill="${sel ? 'var(--sel)' : (g.side==='src'?'var(--src)':'var(--tgt)')}"
           font-weight="${sel?700:400}">${esc(g.gene_id)}</text>`;
    samples.forEach((sm, j) => {
      const v = mats[i][j];
      s += `<rect class="gene" x="${x}" y="${padT + j*ch}" width="${cw-1}" height="${ch-1}"
             fill="${colour(v, lo, hi, div)}" data-gene="${esc(g.gene_id)}"
             data-tip="<b>${esc(g.gene_id)}</b><br>${esc(sm)}: ${v==null?'NA':v.toFixed(3)}"></rect>`;
    });
    if (sel)
      s += `<rect x="${x-1}" y="${padT-2}" width="${cw+1}" height="${samples.length*ch+3}"
             fill="none" stroke="var(--sel)" stroke-width="1.6"></rect>`;
  });
  samples.forEach((sm, j) => {
    s += `<text class="axis" x="${padL-7}" y="${padT + j*ch + ch/2 + 3}" text-anchor="end">${esc(sm)}</text>`;
  });
  // colour key
  const kx = padL + genes.length*cw + 22, kh = Math.min(120, samples.length*ch);
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
    el.addEventListener('click', () => { SEL = (SEL === el.dataset.gene) ? null : el.dataset.gene; redraw(); });
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
    if (!v) return;
    const ls = side === 'src' ? (linked.get(g.gene_id) || []) : DATA.links.filter(l => l.tgt_gene === g.gene_id);
    rows.push({gene:g.gene_id, side, seqid:g.seqid, start:g.start, end:g.end,
               strand:g.strand,
               partner: ls.map(l => (side==='src' ? l.tgt_gene : l.src_gene) || '(unannot)').join(','),
               evidence: [...new Set(ls.map(l => l.track))].join('+'),
               vals:v, samples});
  };
  DATA.source_genes.forEach(g => add(g, 'src', DATA.samples_src, DATA.expr_src));
  DATA.target_genes.forEach(g => add(g, 'tgt', DATA.samples_tgt, DATA.expr_tgt));
  return SEL ? rows.filter(r => r.gene === SEL ||
      DATA.links.some(l => (l.src_gene===SEL && l.tgt_gene===r.gene) ||
                           (l.tgt_gene===SEL && l.src_gene===r.gene))) : rows;
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
    h += `<tr class="${r.gene===SEL?'sel':''}" data-gene="${esc(r.gene)}">`
       + `<td>${esc(r.gene)}</td><td>${r.side}</td><td>${esc(r.seqid)}</td>`
       + `<td class="num">${r.start.toLocaleString()}</td><td class="num">${r.end.toLocaleString()}</td>`
       + `<td>${r.strand}</td><td>${esc(r.partner)}</td><td>${esc(r.evidence)}</td>`;
    for (let i = 0; i < sampleCols.length; i++){
      const v = r.vals[i];
      h += `<td class="num">${v==null?'':v.toFixed(2)}</td>`;
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
    SEL = (SEL === tr.dataset.gene) ? null : tr.dataset.gene; redraw();
  }));
  document.getElementById('tblNote').textContent =
    SEL ? `filtered to ${esc(SEL)} and its linked partners (${rows.length} rows)`
        : `click a gene above to filter — showing all ${rows.length} rows`;
}

function drawRegions(){
  const t = document.getElementById('regTbl');
  if (!DATA.regions.length){ t.innerHTML = '<tbody><tr><td class="empty">none</td></tr></tbody>'; return; }
  const cols = ['rank','tgt_seqid','aligned_bp','n_aln_blocks','pct_id','strand','n_tgt_genes','n_miniprot','n_rbh','tgt_start','tgt_end','tgt_span'];
  let h = '<thead><tr>' + cols.map(c=>`<th>${c}</th>`).join('') + '</tr></thead><tbody>';
  for (const r of DATA.regions)
    h += '<tr>' + cols.map(c=>`<td class="${typeof r[c]==='number'?'num':''}">${esc(r[c])}</td>`).join('') + '</tr>';
  t.innerHTML = h + '</tbody>';
}

function downloadCsv(){
  const rows = tableRows();
  const sampleCols = (DATA.samples_src.length >= DATA.samples_tgt.length ? DATA.samples_src : DATA.samples_tgt);
  const head = ['gene','side','seqid','start','end','strand','partner','evidence'].concat(sampleCols);
  const lines = [head.join(',')];
  for (const r of rows)
    lines.push([r.gene,r.side,r.seqid,r.start,r.end,r.strand,'"'+r.partner+'"',r.evidence]
      .concat(sampleCols.map((_,i)=> r.vals[i]==null?'':r.vals[i])).join(','));
  const b = new Blob([lines.join('\\n')], {type:'text/csv'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(b);
  a.download = DATA.qtl_id.replace(/[^A-Za-z0-9._-]/g,'_') + '_expression.csv';
  a.click();
}

function redraw(){ drawSynteny(); drawHeat(); drawTable(); }

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
document.getElementById('btnReset').onclick = () => { SEL = null; redraw(); };
document.getElementById('btnCsv').onclick = downloadCsv;

drawRegions(); redraw();
</script></body></html>
"""


# --------------------------------------------------------------------------
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
    ap.add_argument("--paf", help="source QTL slice vs target")
    ap.add_argument("--source-expr")
    ap.add_argument("--target-expr")
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank", type=int, default=0)
    args = ap.parse_args()

    regions = [r for r in read_tsv(args.regions) if r["qtl_id"] == args.qtl_id]
    links   = [l for l in read_tsv(args.links)   if l["qtl_id"] == args.qtl_id]
    genes   = [g for g in read_tsv(args.genes)   if g["qtl_id"] == args.qtl_id]

    reg_out, tgt_seqids = [], []
    for r in regions:
        reg_out.append(dict(
            rank=int(r["rank"]), tgt_seqid=r["tgt_seqid"],
            aligned_bp=int(r["aligned_bp"]), n_aln_blocks=int(r["n_aln_blocks"]),
            pct_id=float(r["pct_id"]), strand=r["strand"],
            n_tgt_genes=int(r.get("n_tgt_genes", 0) or 0),
            n_miniprot=int(r.get("n_miniprot", 0) or 0),
            n_rbh=int(r.get("n_rbh", 0) or 0),
            tgt_start=int(r["tgt_start"]), tgt_end=int(r["tgt_end"]),
            tgt_span=int(r["tgt_span"])))
        tgt_seqids.append(r["tgt_seqid"])
    reg_out.sort(key=lambda x: x["rank"])

    link_out = []
    for l in links:
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
                 for g in genes if g["side"] == "target"]

    view_start = max(args.qtl_start - args.flank, 0)
    view_end = args.qtl_end + args.flank

    samples_src, expr_src_all = read_expression(args.source_expr)
    samples_tgt, expr_tgt_all = read_expression(args.target_expr)
    keep_s = {g["gene_id"] for g in src_genes}
    keep_t = {g["gene_id"] for g in tgt_genes}
    expr_src = {k: v for k, v in expr_src_all.items() if k in keep_s}
    expr_tgt = {k: v for k, v in expr_tgt_all.items() if k in keep_t}

    alignments = read_paf(args.paf, tgt_seqids)

    data = dict(
        qtl_id=args.qtl_id,
        source=dict(chrom=args.qtl_chrom, start=view_start, end=view_end,
                    qtl_start=args.qtl_start, qtl_end=args.qtl_end),
        regions=reg_out, links=link_out,
        source_genes=src_genes, target_genes=tgt_genes,
        samples_src=samples_src, samples_tgt=samples_tgt,
        expr_src=expr_src, expr_tgt=expr_tgt,
        alignments=alignments)

    n_mp = sum(1 for l in link_out if l["track"] == "miniprot")
    n_rb = sum(1 for l in link_out if l["track"] == "rbh")
    sub = (f"{args.qtl_chrom}:{args.qtl_start:,}-{args.qtl_end:,} "
           f"&middot; {len(src_genes)} source genes in QTL "
           f"&middot; {len(reg_out)} syntenic regions ({len(alignments)} DNA blocks) "
           f"&middot; {len(tgt_genes)} target genes in regions "
           f"&middot; {n_mp} miniprot, {n_rb} RBH links")

    html = (HTML.replace("__TITLE__", args.qtl_id)
                .replace("__SUBTITLE__", sub)
                .replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    with open(args.out, "w") as fh:
        fh.write(html)
    sys.stderr.write(f"[build_synteny_html] {args.out} ({len(src_genes)} src, "
                     f"{len(tgt_genes)} tgt genes, {n_mp} miniprot, {n_rb} rbh, "
                     f"{len(alignments)} aln blocks)\n")


if __name__ == "__main__":
    main()
