#!/usr/bin/env python3
"""Synthetic two-variety dataset for the synteny viewer.

Source: one 5 Mb chromosome '1B' with genes every 50 kb.
Target: 6 small scaffolds (mimicking a fragmented assembly) carrying
        orthologues of the source QTL genes, deliberately scattered so the
        multi-scaffold display is exercised.
"""
import os, random, sys
random.seed(11)
OUT = sys.argv[1] if len(sys.argv) > 1 else "test/data"
os.makedirs(OUT, exist_ok=True)

AA = "ACDEFGHIKLMNPQRSTVWY"
def dna(n): return "".join(random.choice("ACGT") for _ in range(n))
def wrap(s, w=60): return "\n".join(s[i:i+w] for i in range(0, len(s), w))

SRC_LEN = 5_000_000
QS, QE = 2_455_000, 4_500_000
src = dna(SRC_LEN)

src_genes = []
for i, p in enumerate(range(1_000_000, SRC_LEN - 500_000, 50_000)):
    src_genes.append((f"NorGene{i:03d}", p, p + 8_000, random.choice("+-")))

scafs = [f"cad_scaffold_{n:06d}" for n in (8738, 1122, 39029, 6409, 19105, 91312)]
scaf_len = {s: 250_000 for s in scafs}
tgt_seq = {s: dna(scaf_len[s]) for s in scafs}

# orthologues: QTL genes distributed across scaffolds, most on the first
qtl_genes = [g for g in src_genes if g[1] < QE and g[2] > QS]
tgt_genes, ortho = [], {}
counts = [14, 9, 6, 4, 3, 2]
gi, k = 0, 0
for si, sc in enumerate(scafs):
    for j in range(counts[si]):
        if k >= len(qtl_genes):
            break
        p = 20_000 + j * 14_000
        gid = f"CadGene{gi:03d}"
        tgt_genes.append((gid, sc, p, p + 7_000, random.choice("+-")))
        ortho[qtl_genes[k][0]] = gid
        gi += 1
        k += 1

# implant homologous DNA so minimap2 finds real blocks
tl = {s: list(tgt_seq[s]) for s in scafs}
for sgene, tgene in ortho.items():
    s = next(g for g in src_genes if g[0] == sgene)
    t = next(g for g in tgt_genes if g[0] == tgene)
    block = src[s[1] - 4000: s[1] + 12_000]
    block = "".join(c if random.random() > 0.04 else random.choice("ACGT") for c in block)
    st = max(t[2] - 4000, 0)
    tl[t[1]][st:st + len(block)] = list(block)
tgt_seq = {s: "".join(tl[s])[:scaf_len[s]] for s in scafs}

with open(f"{OUT}/source.fa", "w") as f:
    f.write(f">1B\n{wrap(src)}\n")
with open(f"{OUT}/target.fa", "w") as f:
    for s in scafs:
        f.write(f">{s}\n{wrap(tgt_seq[s])}\n")

with open(f"{OUT}/source.gff3", "w") as g, open(f"{OUT}/source.pep.fa", "w") as p:
    g.write("##gff-version 3\n")
    for (gid, s, e, strand) in src_genes:
        g.write(f"1B\t.\tgene\t{s+1}\t{e}\t.\t{strand}\t.\tID={gid}\n")
        g.write(f"1B\t.\tmRNA\t{s+1}\t{e}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}\n")
        g.write(f"1B\t.\tCDS\t{s+1}\t{e}\t.\t{strand}\t0\tID={gid}.cds;Parent={gid}.t1\n")
        p.write(f">{gid}.t1\nM{''.join(random.choice(AA) for _ in range(199))}\n")

# target proteins: orthologues share most of the source sequence
src_pep = {}
for line in open(f"{OUT}/source.pep.fa"):
    if line.startswith(">"): cur = line[1:].strip().split(".")[0]
    else: src_pep[cur] = line.strip()
rev = {v: k for k, v in ortho.items()}
with open(f"{OUT}/target.gff3", "w") as g, open(f"{OUT}/target.pep.fa", "w") as p:
    g.write("##gff-version 3\n")
    for (gid, sc, s, e, strand) in tgt_genes:
        g.write(f"{sc}\t.\tgene\t{s+1}\t{e}\t.\t{strand}\t.\tID={gid}\n")
        g.write(f"{sc}\t.\tmRNA\t{s+1}\t{e}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}\n")
        g.write(f"{sc}\t.\tCDS\t{s+1}\t{e}\t.\t{strand}\t0\tID={gid}.cds;Parent={gid}.t1\n")
        src_of = rev.get(gid)
        if src_of and src_of in src_pep:
            base = list(src_pep[src_of])
            for _ in range(int(len(base) * 0.12)):
                base[random.randrange(len(base))] = random.choice(AA)
            seq = "".join(base)
        else:
            seq = "M" + "".join(random.choice(AA) for _ in range(199))
        p.write(f">{gid}.t1\n{seq}\n")

with open(f"{OUT}/qtl.bed", "w") as f:
    f.write(f"1B\t{QS}\t{QE}\t1B1:Norin:Snn1/SnTox1:LOD-76\n")
    f.write("1B\t600000\t900000\t1B_control:no_target\n")

def expr(path, genes, samples):
    with open(path, "w") as f:
        f.write("gene_id\t" + "\t".join(samples) + "\n")
        for g in genes:
            f.write(g + "\t" + "\t".join(f"{max(0, random.gauss(25, 18)):.2f}"
                                        for _ in samples) + "\n")
expr(f"{OUT}/expr_source.tsv", [g[0] for g in src_genes],
     ["Nor_mock_1", "Nor_mock_2", "Nor_inf_1", "Nor_inf_2", "Nor_inf_3"])
expr(f"{OUT}/expr_target.tsv", [g[0] for g in tgt_genes],
     ["Cad_mock_1", "Cad_mock_2", "Cad_inf_1", "Cad_inf_2", "Cad_inf_3"])

print(f"{OUT}: {len(src_genes)} source genes ({len(qtl_genes)} in QTL), "
      f"{len(tgt_genes)} target genes on {len(scafs)} scaffolds, "
      f"{len(ortho)} orthologue pairs")
