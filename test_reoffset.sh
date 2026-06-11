#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

HMM="MET-test.hmm"
FASTA="met-optimized.fa"
BUILD="build"

run() {
    local flags="$1"
    local label="$2"
    echo "=== $label ==="
    cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="$flags" > /dev/null 2>&1
    if ! make -C "$BUILD" -j$(nproc) dummer > /dev/null 2>&1; then
        echo "BUILD FAILED with flags: $flags"
        exit 1
    fi
    "$BUILD/dummer" "$HMM" "$FASTA" 2>&1 | grep "^a score="
    echo
}

run ""                                    "baseline"
run "-DDEBUG_FORCE_REOFFSET"              "force-reoffset"
run "-DDEBUG_REOFFSET -DDEBUG_FORCE_REOFFSET" "force+per-lane-offsets"

echo "--- plain FASTA smoke test ---"
if ! "$BUILD/dummer" -t10 -l200 -b0 tests/dfam-test.hmm tests/dna-test.fa > /dev/null 2>&1; then
    echo "plain FASTA: BUILD FAILED"
    exit 1
fi
echo "plain FASTA: OK"

echo
echo "all tests passed"
