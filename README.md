# qtl-synteny

Interactive pairwise synteny + expression views of QTL regions between two
varieties, for manual inspection.

Built for the case where the **target assembly is fragmented** (no
chromosome-scale sequences), so collapsing each QTL to a single transferred
interval is not meaningful. Instead of forcing one answer, this shows you the
source QTL against its top-ranked candidate target sequences side by side,
with genes, sequence alignment, and expression, and lets you judge.

## How it works

Four layered steps. **DNA synteny defines the regions; protein evidence is
layered on afterwards** — so regions are still found when the protein
evidence comes back empty, which happens on a fragmented target where
reciprocal-best-hit loses to homoeologous competition.

| Step | What | Tool | Role |
|---|---|---|---|
| 1 | QTL interval vs the **whole** target assembly | minimap2 | finds region hits |
| 2 | rank the target sequences hit, by aligned bp, keep top N | — | **defines the syntenic regions** |
| 3 | source proteins spliced onto the target **genome**, kept if they land inside those regions | miniprot | homology; needs no target annotation |
| 4 | reciprocal best hits between the two **annotations**, kept if the target gene is inside those regions | MMseqs2 `easy-rbh` | orthology; needs both annotations |

Steps 3 and 4 are shown as separate, independently toggleable tracks —
miniprot as dashed purple curves, RBH as solid orange. A miniprot hit with no
annotated target gene under it is labelled *(unannotated locus)*, which is
exactly the case where the target annotation is missing something real.

Genes displayed, and used by the expression panel, are **all source genes in
the QTL plus all target genes in the syntenic regions** — including those
with no protein link at all, since a target gene with no orthologue call is
still worth seeing expression for.

## What each view shows

**Synteny panel** — a **coverage strip** above the source axis shows which
parts of the QTL any displayed target explains, so unexplained gaps are
visible at a glance. Below it the source QTL with gene arrows to scale, then
each syntenic target as a lane labelled with the share of the QTL it covers,
its cumulative contribution, % identity, and gene / miniprot / RBH counts.
- **Grey ribbons** = DNA alignment blocks (orange-tinted when inverted)
- **Dashed purple curves** = miniprot placements
- **Solid orange curves** = reciprocal best hits
- Toggle each of the three independently, or raise the identity floor.

**Expression heatmap** — genes on X, samples on Y by default, with a
**transpose** toggle to put genes on Y (useful with many samples: your 24
Cadenza libraries make the default layout wide). log2(x+1), raw, or per-gene
z-score; restrict to source, target, or only genes carrying links.

**Expression table** — per-gene values with an `evidence` column showing
which tracks support each gene. Sortable, CSV-downloadable.

Clicking a gene anywhere pins it across all three panels.

**Click a target's label (or its row in the regions table) to isolate it** —
every other region fades (ribbons, links, genes, and its own label) rather
than disappearing, so you keep the full picture for context while inspecting
one candidate at a time. Only that region's source genes and target genes
feed the expression heatmap and table while isolated. "Show all regions"
clears it. Combines with gene selection: isolate a region, then click a gene
within it to pin that gene specifically.

Pages are **fully self-contained** — data embedded as JSON, all CSS/JS
inline, no CDN. They work from an offline HPC filesystem.

## Quick start

```bash
nextflow run . \
    --source_genome     norin61.fa \
    --target_genome     cadenza.fa \
    --source_gff        norin61.gff3 \
    --target_gff        cadenza.gff3 \
    --source_proteins   norin61.pep.fa \
    --target_proteins   cadenza.pep.fa \
    --qtl_bed           qtl_intervals.bed \
    --source_expression norin_expression.tsv \
    --target_expression cadenza_expression.tsv \
    --source_name Norin --target_name Cadenza \
    --outdir results \
    -profile setonix
```

Then open `results/views/index.html` and click through.

## Inputs

| Parameter | Required | Description |
|---|---|---|
| `--source_genome` | ✅ | assembly the QTLs are defined on |
| `--target_genome` | ✅ | assembly to compare against |
| `--source_gff` / `--target_gff` | ✅ | annotations for both |
| `--source_proteins` / `--target_proteins` | ✅ | proteomes for both |
| `--qtl_bed` | ✅ | `chrom  start  end  qtl_id` (BED, 0-based) |
| `--source_expression` / `--target_expression` | optional | `gene_id <TAB> sample1 <TAB> sample2 …` |
| `--source_name` / `--target_name` | optional | display names |

### Gene IDs and the counts table

Counts tables are `gene_id` in column 1, samples across the header
(genes down the rows), e.g.

```
Genes                               Cad_PN143_0h_R1  Cad_PN143_0h_R2  ...
TraesCAD_scaffold_000001_01G000100  0                0
```

Ensembl/EI GFF3 prefixes its IDs with the feature type
(`ID=gene:TraesCAD_...`, `ID=transcript:TraesCAD_....1`) while counts tables
and protein FASTAs use the bare accession. **These prefixes are now stripped
automatically** — `gene:`, `transcript:`, `mRNA:`, `CDS:`, `protein:`,
`exon:`, `rna:`, `ncRNA:` — so no flag is needed. An explicit `gene_id`
attribute on the transcript is preferred when present.

Matching then tries, in order: the gene id as-is, with any type prefix
removed, with a trailing `.1` removed, and with `.1` added. The log reports
`expression source: N/M genes matched` per QTL and warns loudly if nothing
matched, so an ID mismatch is visible rather than showing an empty heatmap.

Use `--source_gff_*` / `--target_gff_*` if the two annotations follow
different conventions (`gff_feature`, `gff_id_attr`, `gff_strip_prefix`);
each falls back to the shared `--gff_*` value.

## Region selection: report everything, cap only by count

Selection is deliberately simple: **every target sequence that aligns is
reported**, after basic identity/length filters, ranked by how much of the
source interval it covers. The only way to limit the list is
`--top_n_regions` (default `0` = unlimited). This is an inspection tool
meant to work on assemblies of any quality — a handful of chromosomes or
hundreds of thousands of scaffolds — so it shows everything by default
rather than deciding on your behalf.

```bash
nextflow run . --top_n_regions 0 ...    # default: report everything found
nextflow run . --top_n_regions 20 ...   # cap at 20 per QTL
nextflow run . --min_pct_id 97 ...      # drop obvious homoeologs before ranking
```

`--region_min_aligned` (default `500` bp) and `--region_min_mapq` are basic
noise filters, not selection strategies — they exist to drop spurious tiny
alignments, not to decide which real hits matter.

## Paralogy grouping: possible genome copies

On a polyploid genome the same source region legitimately has several real
copies among the candidates (homoeologs). Rather than pick one and discard
the rest, candidates whose covered **source** interval overlaps by at least
`--group_overlap_frac` (default `0.5`, reciprocal) are grouped as possible
copies of one another (`region_group`).

**Which one is the best match is decided by protein evidence, not DNA
coverage.** Within each group, the scaffold with the most combined RBH +
miniprot support is marked `is_primary_in_group`; DNA coverage alone only
breaks ties. The scaffold with the most DNA coverage is not always the true
orthologue — a homoeolog can align just as well or better while carrying
weaker or no gene-level support. In the viewer the primary gets a `★ best
match` label; the others get `↺ possible copy of <primary>`.

## Visualization: hide gene-less scaffolds

The synteny viewer only draws lanes for scaffolds carrying **at least one**
RBH or miniprot link. A scaffold with real DNA homology but zero protein
evidence has nothing to inspect at the gene level, and with unlimited
selection as the default there could be many of these cluttering the view.
They are not lost — `regions.tsv` and the FASTA output in
`selected_scaffolds/` still contain every mapped scaffold regardless. Use
`--show_empty_regions true` to draw them anyway.

## Two coverage badges

Every view reports source-interval coverage split by whether the covering
scaffold has any gene evidence:

- **covered by scaffolds WITH gene evidence** — union DNA coverage from
  scaffolds carrying ≥1 RBH or miniprot hit
- **covered ONLY by gene-less scaffolds** — union coverage that comes
  exclusively from scaffolds with no protein evidence
- **total DNA coverage** — union of both, any evidence at all

These are computed from **every** mapped scaffold, not just what
`--top_n_regions` kept or what the viewer draws — they answer "how much of
this QTL is explainable at all", independent of any display filtering.

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--top_n_regions` | `0` | max regions kept per QTL, ranked by source coverage; `0` = unlimited |
| `--min_pct_id` | `0` | drop targets below this alignment identity before ranking (try 97-98 for wheat) |
| `--region_min_aligned` | `500` | bp floor, noise filter only |
| `--region_min_mapq` | `0` | MAPQ floor, noise filter only |
| `--group_overlap_frac` | `0.5` | reciprocal source-overlap fraction to group scaffolds as possible genome copies |
| `--show_empty_regions` | `false` | draw lanes for scaffolds with no RBH/miniprot evidence too |
| `--homoeolog_min_frac` | `0.6` | fraction of a scaffold's RBH partners that must sit on the QTL's chromosome to call it orthologous |
| `--homoeolog_evidence` | `both` | `both`, `rbh`, or `miniprot` (the last needs no target annotation) |
| `--drop_homoeologs` | `false` | omit flagged homoeologs instead of marking them |
| `--scaffold_gene_priority` | `true` | list/write scaffolds carrying annotated genes first in `selected_scaffolds/`; ties broken by rank |
| `--region_flank` | `20000` | bp padding around each target window |
| `--source_flank` | `100000` | bp either side of the QTL, sliced and displayed |
| `--region_minimap_args` | `-x asm20` | preset for the per-QTL alignment |
| `--minimap2_index_args` | `-x asm20` | preset for the target index |
| `--miniprot_min_ident` | `0.5` | identity floor for an in-region miniprot placement |
| `--skip_miniprot` | `false` | drop step 3 |
| `--skip_rbh` | `false` | drop step 4 (then `--target_proteins` is not needed) |
| `--mmseqs_args` | `-s 5.7 -e 1e-5` | MMseqs2 `easy-rbh` options |
| `--mmseqs_also_search` | `false` | additionally run `easy-search` for non-RBH links |
| `--mmseqs_search_args` | `-s 5.7 --max-seqs 50 -e 1e-5` | options for that extra search |

Regions are ranked by **total aligned bp** from the DNA alignment, which
needs no gene or ortholog evidence at all.

## Output

```
results/
├── selected_scaffolds/
│   ├── *.scaffolds.fa            whole selected scaffolds, headers carry rank/n_genes/coverage/call
│   ├── *.windows.fa              just the displayed window of each
│   └── *.selected_scaffolds.tsv  plain list, `list_order` column reflects gene-priority
├── views/
│   ├── index.html                    landing page, one row per QTL
│   └── <qtl_id>.html                 interactive view per QTL
├── regions/
│   ├── *.dna_regions.tsv             step 2: every DNA-ranked region, with region_group
│   ├── *.regions.tsv                 + protein evidence, is_primary_in_group
│   ├── *.coverage.tsv                source coverage split: with/without gene evidence
│   ├── *.links.tsv                   protein links, `track` = miniprot | rbh
│   ├── *.genes.tsv                   genes displayed per QTL, both sides
│   └── paf/*.region.paf              per-QTL DNA alignments
├── annotation/*.genes.tsv
├── homology/*.miniprot.gff           step 3
├── orthology/*.rbh.m8                step 4
└── pipeline_info/
```

## Why this is scattered per QTL

The expensive whole-genome alignment is avoided entirely. Each QTL only
triggers a minimap2 run between its own source interval (a few Mb) and its
handful of candidate target windows (a few hundred kb each) — small jobs that
parallelise across the cluster instead of one huge memory-bound job. This is
what makes the approach viable on a fragmented multi-Gb target where
whole-genome alignment ran out of memory.

The genome-scale steps — the minimap2 target index, the miniprot placement,
and the MMseqs2 `easy-rbh` comparison — each run **once** and are reused by
every QTL.

**On a fragmented target, consider `--mmseqs_also_search true`.** Reciprocal
best hits are strict: a true orthologue that happens to lose reciprocity to a
tandem duplicate produces no link at all, and on a shattered assembly that can
leave a real candidate region invisible. The extra search adds those hits as
non-RBH links, clearly distinguished in the viewer, at the cost of one more
proteome search.

## Profiles

`docker`, `singularity`, `conda`, `slurm`, `pawsey_setonix`, `test`. Combine
them, e.g. `-profile pawsey_setonix,singularity` (matching the qtl-liftover
repo's `run_*.sbatch`).

The `pawsey_setonix` profile inherits the nf-core Pawsey config and
deliberately binds only `/scratch,/software` — **not** `/group`, which is not
mounted on every allocation and causes a fatal container-creation error when
hardcoded.

Resource labels match the qtl-liftover repo: `process_high` is 200 GB and
`process_extra_high` 500 GB, both clamped by `--max_memory`.

## Config layering gotcha

Tool arguments live in `conf/modules.config` as `ext.args`; container images
and `publishDir` live in `conf/containers.config`.

**A `withName` block inside a profile-included config replaces that
selector's entire entry**, silently wiping the `ext.args` set in
`conf/modules.config` — the process then runs with no arguments and no error.
So all profile configs set resources with `withLabel`, never `withName`.
`check_max()` is defined in `nextflow.config` (the parent) so every included
file can call it.

## Tests

```bash
./test/run_test.sh        # end-to-end on synthetic data, no aligners needed
nextflow run . -stub-run -profile test
```

The test builds a 5 Mb source chromosome and a deliberately fragmented
6-scaffold target with implanted homologous blocks, runs the whole pipeline
with shims for MMseqs2 and minimap2, and asserts the views contain embedded
data and the index links resolve. It also includes a control QTL with no
target hits, to confirm that degrades to a source-only view rather than
failing.

## Layout

```
main.nf                          thin entry point
nextflow.config                  params, profiles, check_max
nfcore_custom.config             resources by process label
conf/
  containers.config              containers + publishDir
  modules.config                 ext.args
  test.config  slurm.config  profiles/
workflows/qtl-synteny.nf         QTL_SYNTENY_FLOW
subworkflows/local/
  prepare_annotation.nf          gene tables + MMseqs2 both directions
  synteny_views.nf               region selection, per-QTL align, viewers
modules/local/                   7 processes
bin/                             python scripts (auto-added to PATH)
test/                            data generator, shims, run_test.sh
```

> If the executable bit was lost in transfer, run
> `chmod +x bin/*.py test/*.py test/*.sh test/shims/*` — Nextflow requires
> `bin/` scripts to be executable.
