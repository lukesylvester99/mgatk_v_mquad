#!/usr/bin/env bash
set -euo pipefail

# Args:
#   1 MQUAD_OUT   (directory that contains MQuad outputs)
#   2 VIERO_PY    (path to fit_mito_clones.py)

MQUAD_OUT="$1"
VIERO_PY="$2"

[[ -d "${MQUAD_OUT}" ]] || { echo "ERROR: Missing MQUAD_OUT dir: ${MQUAD_OUT}"; exit 1; }
[[ -f "${VIERO_PY}"  ]] || { echo "ERROR: Missing vireo script: ${VIERO_PY}"; exit 1; }

VIERO_OUT="${MQUAD_OUT}/vireo"
mkdir -p "${VIERO_OUT}"

echo "[INFO] MQUAD_OUT: ${MQUAD_OUT}"
echo "[INFO] VIERO_PY : ${VIERO_PY}"
echo "[INFO] VIERO_OUT: ${VIERO_OUT}"

echo "Running vireo script in ${MQUAD_OUT} ..."
pushd "${MQUAD_OUT}" >/dev/null
python3 "${VIERO_PY}" > "${VIERO_OUT}/fit_mito_clones.log" 2>&1
popd >/dev/null

echo "fit_mito_clones.py done."
