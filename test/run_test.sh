#!/usr/bin/env bash
# Integration test on synthetic data. mmseqs and minimap2 are replaced by
# shims, so no aligners need to be installed.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
NF="${NEXTFLOW:-nextflow}"

# coordinate-logic regression tests first - these assert exact coverage
# numbers and region boundaries, which the integration test below cannot
echo "=== coordinate-logic regression tests ==="
python3 "$HERE/test_coverage_logic.py"
echo

python3 "$HERE/make_test_data.py" "$HERE/data"
export PATH="$HERE/shims:$PATH"

cd "$ROOT"
"$NF" run . -profile test "$@"

echo; echo "=== syntenic regions (DNA-ranked, protein evidence attached) ==="
cat "$HERE"/results/regions/*.regions.tsv
echo; echo "=== views produced ==="
ls -1 "$HERE"/results/views/

python3 - "$HERE/results" <<'PY'
import glob, os, re, sys, json
d = sys.argv[1]
views = [v for v in glob.glob(f"{d}/views/*.html") if not v.endswith("index.html")]
assert views, "no per-QTL views produced"
for v in views:
    h = open(v).read()
    assert "__DATA__" not in h and "__TITLE__" not in h, f"unreplaced placeholder in {v}"
    m = re.search(r"const DATA = (\{.*?\});\n", h, re.S)
    assert m, f"no embedded DATA in {v}"
    data = json.loads(m.group(1))
    mp = sum(1 for l in data['links'] if l['track'] == 'miniprot')
    rb = sum(1 for l in data['links'] if l['track'] == 'rbh')
    print(f"{os.path.basename(v)}: {len(data['source_genes'])} src genes, "
          f"{len(data['target_genes'])} tgt genes, {len(data['regions'])} regions, "
          f"{len(data['alignments'])} DNA blocks, {mp} miniprot, {rb} rbh, "
          f"{len(data['expr_src'])}+{len(data['expr_tgt'])} expr rows")
idx = f"{d}/views/index.html"
assert os.path.exists(idx), "no index.html"
ih = open(idx).read()
assert "href=" in ih, "index has no links to views"
print("\nPASS: views built, data embedded, index links resolve")
PY
