#!/bin/bash
# ============================================================================
# ShardPlanner CI test — assert the contract the pbs script's cost-balanced sharding relies on.
#
# longtr_joint_call.pbs cuts the catalog into shards using ShardPlanner's per-locus cost, then
# `die`s if the shard files don't contain exactly the original locus count. That guard only helps
# if the planner itself preserves the catalog, so this checks the properties the partitioner
# actually depends on:
#
#   1. exits 0 and emits a '#' header
#   2. one output line per catalog locus -- no loci gained or lost
#   3. the BED payload is passed through VERBATIM and in order, so concatenating the shards
#      reproduces the catalog byte-for-byte (which is what keeps the final VCF concat contiguous)
#   4. costs are non-negative integers, and not all zero -- an all-zero plan gives the partitioner
#      nothing to balance on and silently degenerates to one shard
#
# Reads only the BAM headers and .bai, so it is fast regardless of fixture size.
#
# Usage:  run_shard_planner_test.sh /path/to/ShardPlanner
# ============================================================================
set -euo pipefail
PLANNER="${1:?usage: run_shard_planner_test.sh /path/to/ShardPlanner}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT

"$PLANNER" --bams "$DIR/HG00097.bam,$DIR/HG00102.bam" --regions "$DIR/loci.bed" > "$OUT/plan.tsv"

# 1. header present
head -1 "$OUT/plan.tsv" | grep -q '^#ShardPlanner' \
  || { echo "PLANNER FAIL: missing '#ShardPlanner' header line" >&2; exit 1; }

# 2. one line per locus
want=$(grep -cvE '^(#|$)' "$DIR/loci.bed")
got=$(grep -vc '^#' "$OUT/plan.tsv" || true)
[ "${got:-0}" -eq "$want" ] \
  || { echo "PLANNER FAIL: expected $want locus lines, got ${got:-0}" >&2; exit 1; }

# 3. BED payload preserved verbatim and in order (cut the cost column back off)
grep -v '^#' "$OUT/plan.tsv" | cut -f2- > "$OUT/roundtrip.bed"
grep -vE '^(#|$)' "$DIR/loci.bed" > "$OUT/expect.bed"
cmp -s "$OUT/roundtrip.bed" "$OUT/expect.bed" \
  || { echo "PLANNER FAIL: BED payload not round-tripped verbatim" >&2
       diff "$OUT/expect.bed" "$OUT/roundtrip.bed" >&2 || true; exit 1; }

# 4. costs are non-negative integers and not uniformly zero
awk -F'\t' '!/^#/ {
              if ($1 !~ /^[0-9]+$/) { print "PLANNER FAIL: non-integer cost: " $1 > "/dev/stderr"; exit 1 }
              tot += $1
            }
            END { if (tot == 0) { print "PLANNER FAIL: every locus cost 0 -- nothing to balance on" > "/dev/stderr"; exit 1 } }' \
    "$OUT/plan.tsv"

echo "PLANNER OK: $got loci, header + verbatim round-trip + non-zero costs"
