#!/usr/bin/env bash
# Smoke test for the microbiomeutils container image.
#
# Reruns the README examples inside the image against the example data in
# data/ and compares the results to the committed reference outputs.
# The container runs with a read-only root filesystem and as an arbitrary
# (non-image) user, as it would under Apptainer on an HPC cluster.
#
# Usage: .github/scripts/smoke-test.sh IMAGE   (run from the repository root)
set -euo pipefail

image=${1:?usage: $0 IMAGE}
data="$PWD/data"
out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT

run() {
    docker run --rm --read-only --tmpfs /tmp \
        --user "$(id -u):$(id -g)" \
        -v "$data:/data:ro" -v "$out:/out" \
        "$image" "$@"
}

check() {  # check NAME ACTUAL EXPECTED
    if cmp "$2" "$3"; then
        echo "ok   $1"
    else
        echo "FAIL $1: $2 differs from $3" >&2
        exit 1
    fi
}

run -h > /dev/null

# Distance matrices must be byte-identical to the reference outputs.
run betadiv -i feature_table.tsv -m bray-curtis > "$out/bc.tsv"
check bray-curtis "$out/bc.tsv" "$data/bc_res.tsv"

run betadiv -i feature_table.tsv -m weighted-unifrac -t tree.fasttree > "$out/wuf.tsv"
check weighted-unifrac "$out/wuf.tsv" "$data/wuf_res.tsv"

run betadiv -i feature_table.tsv -m unweighted-unifrac -t tree.fasttree > "$out/uwuf.tsv"
check unweighted-unifrac "$out/uwuf.tsv" "$data/uwuf_res.tsv"

# Taxonomic summaries.
for level in 1 2 3 4 5 6 7; do
    run taxsum -i feature_table.tsv -l "$level" > "$out/taxonomy_L$level.tsv"
    check "taxsum L$level" "$out/taxonomy_L$level.tsv" "$data/taxonomy_L$level.tsv"
done
run taxsum -i feature_table.tsv -l 6 -t relative > "$out/taxonomy_L6_relative.tsv"
check "taxsum L6 relative" "$out/taxonomy_L6_relative.tsv" "$data/taxonomy_L6_relative.tsv"

# PCoA coordinates are not compared byte-for-byte: eigenvector signs and the
# last floating-point digits depend on the BLAS library. Check that pcoa runs
# and that its output feeds emperor.
run pcoa -i /out/wuf.tsv > "$out/wuf_coords.tsv"
test -s "$out/wuf_coords.tsv"
echo "ok   pcoa"

run emperor -i /out/wuf_coords.tsv -m mapping_file.tsv -o /out/plot
test -s "$out/plot/index.html"
echo "ok   emperor"

echo "All smoke tests passed."
