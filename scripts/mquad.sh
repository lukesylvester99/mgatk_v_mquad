#!/usr/bin/env bash
set -euo pipefail

# Args:
#   1 CELL_SNP_OUT   (directory from cellsnp-lite mode 1a)
#   2 MQUAD_OUT      (output directory)
#   3 THREADS

CELL_SNP_OUT="$1"
MQUAD_OUT="$2"
THREADS="${3:-20}"

command -v mquad >/dev/null 2>&1 || { echo "ERROR: mquad not found in PATH"; exit 1; }

mkdir -p "${MQUAD_OUT}"

echo "[INFO] CELL_SNP_OUT: ${CELL_SNP_OUT}"
echo "[INFO] MQUAD_OUT   : ${MQUAD_OUT}"
echo "[INFO] THREADS     : ${THREADS}"

echo "Running MQuad ..."
mquad \
  -c "${CELL_SNP_OUT}" \
  -o "${MQUAD_OUT}" \
  -p "${THREADS}" \
  --minDP 5 \
  > "${MQUAD_OUT}/mquad.log" 2>&1

echo "MQuad has completed. Woohoo!"
