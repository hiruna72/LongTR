#!/bin/bash
# ============================================================================
# LongTR CI determinism test — the same input must produce the same VCF.
#
# Guards against doc/bug_poa_random_subsample.md: HaplotypeGenerator::poa() subsamples large read
# clusters, and used a std::random_device-seeded RNG. Identical runs therefore produced different
# candidate haplotypes, different ALT alleles and different genotypes — silently, at ~0.3% of loci.
# The bug survived from 2024-02 to 2026-07 because nothing checked reproducibility.
#
# Fixture (this directory): 2 public 1000G ONT samples over ONE chr21 locus known to exercise the
# subsampling path, plus a chr21 reference with real bases only where the reads land (N elsewhere,
# full length kept so coordinates need no shifting).
#
# THREE assertions, and the third matters as much as the first two:
#   1. same seed, twice            -> byte-identical            (the regression itself)
#   2. default settings, twice     -> byte-identical            (any other nondeterminism)
#   3. two different seeds         -> MUST DIFFER               (proves the fixture still reaches the
#                                                                subsampling path; without this the
#                                                                test passes vacuously if the fixture,
#                                                                the cluster threshold or the
#                                                                admission logic ever changes)
#
# poa() is only reached when a sample has >25% of its reads unadmitted, and only subsamples when a
# cluster reaches --poa-cluster-limit (default 30). This fixture is small, so the test lowers the
# limit to 5 to reach the path; assertion 3 is what verifies that is still true.
#
# Usage:  run_determinism_test.sh /path/to/LongTR
# ============================================================================
set -euo pipefail
LONGTR="${1:?usage: run_determinism_test.sh /path/to/LongTR}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT

COMMON=(--bams "$DIR/HG00097.bam,$DIR/HG00102.bam" --bam-samps HG00097,HG00102 --bam-libs HG00097,HG00102
        --fasta "$DIR/chr21.fa.gz" --regions "$DIR/loci.bed"
        --min-reads 5 --phased-bam --output-gls --output-pls --output-phased-gls --output-filters
        --alignment-params -1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1)

run(){ "$LONGTR" "${COMMON[@]}" "${@:2}" --tr-vcf "$OUT/$1.vcf.gz" > "$OUT/$1.log" 2>&1 \
         || { echo "DETERMINISM FAIL: LongTR exited non-zero for '$1'" >&2; tail -20 "$OUT/$1.log" >&2; exit 1; }
       zcat "$OUT/$1.vcf.gz" | grep -v '^#' > "$OUT/$1.body"
       [ -s "$OUT/$1.body" ] || { echo "DETERMINISM FAIL: '$1' produced 0 VCF records" >&2; exit 1; } }

run sub_a --poa-cluster-limit 5 --seed 42
run sub_b --poa-cluster-limit 5 --seed 42
run sub_c --poa-cluster-limit 5 --seed 7
run def_a
run def_b

fail=0

if cmp -s "$OUT/sub_a.body" "$OUT/sub_b.body"; then
  echo "  ok   1/3  same seed twice -> identical ($(wc -l < "$OUT/sub_a.body") records)"
else
  echo "  FAIL 1/3  same seed twice -> DIFFERS: LongTR is not reproducible" >&2
  diff "$OUT/sub_a.body" "$OUT/sub_b.body" | head -4 | cut -c1-160 >&2
  fail=1
fi

if cmp -s "$OUT/def_a.body" "$OUT/def_b.body"; then
  echo "  ok   2/3  default settings twice -> identical"
else
  echo "  FAIL 2/3  default settings twice -> DIFFERS: nondeterminism outside the POA subsampling" >&2
  diff "$OUT/def_a.body" "$OUT/def_b.body" | head -4 | cut -c1-160 >&2
  fail=1
fi

if cmp -s "$OUT/sub_a.body" "$OUT/sub_c.body"; then
  echo "  FAIL 3/3  --seed 42 and --seed 7 agree -> the fixture no longer reaches the POA" >&2
  echo "            subsampling path, so assertions 1-2 prove nothing. Rebuild the fixture" >&2
  echo "            (see doc/bug_poa_random_subsample.md) or lower --poa-cluster-limit." >&2
  fail=1
else
  echo "  ok   3/3  different seeds -> differ (subsampling path is exercised)"
fi

[ "$fail" -eq 0 ] || { echo "DETERMINISM TEST FAILED" >&2; exit 1; }
echo "DETERMINISM OK"
