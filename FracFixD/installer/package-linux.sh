#!/usr/bin/env bash
# FracFixD Linux release-packaging script.
#
# Usage:  installer/package-linux.sh
#
# What it does:
#   1. Strips the per-microarch binaries in dist/.
#   2. Computes SHA-256 for each artefact.
#   3. Writes ./FracFixD/bin/SHA256SUMS.
#
# Assumes the four binaries have already been built per the
# walk-through in docs/BUILD_LINUX.md and are in ./dist/.
#
# License: MIT
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if [[ ! -d dist ]]; then
    echo "fatal: no dist/ directory; build the binaries first" >&2
    echo "see docs/BUILD_LINUX.md" >&2
    exit 1
fi

ARTEFACTS=(
    "dist/fracfixd-linux-znver2-x86_64"
    "dist/fracfixd-linux-broadwell-x86_64"
    "dist/fracfixd-linux-generic-x86_64"
    "dist/fracfixd-cli-linux-x86_64-static"
)

for a in "${ARTEFACTS[@]}"; do
    if [[ ! -x "$a" ]]; then
        echo "fatal: missing or non-executable: $a" >&2
        exit 1
    fi
    echo "[strip] $a"
    strip --strip-unneeded "$a"
done

mkdir -p bin
cp "${ARTEFACTS[@]}" bin/

echo "[sha256] computing checksums in bin/"
( cd bin && sha256sum * > SHA256SUMS )

echo "[done] FracFixD Linux binaries staged in ./bin/"
ls -l bin/
echo
cat bin/SHA256SUMS
