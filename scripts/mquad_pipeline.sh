#!/usr/bin/env bash
set -euo pipefail

############
# Config #
##########

# Optional positional args (so Snakemake can call this script):
#   1 BAM
#   2 EPI_BARCODES
#   3 DEPTH_BARCODES
#   4 VARFILT_TSV
#   5 OUT_BASE
#   6 THREADS
#   7 MT_CONTIG
#
# If not provided, defaults below are used.
if [[ $# -ge 5 ]]; then
  BAM="$1"
  EPI_BARCODES="$2"
  DEPTH_BARCODES="$3"
  VARFILT_TSV="$4"
  OUT_BASE="$5"
  THREADS="${6:-20}"
  MT_CONTIG="${7:-chrM}"
else
  # Cell Ranger ATAC BAM
  BAM="$HOME/A10_trial/snake_outs/A10/outs/possorted_bam.bam"

  # base directory
  BASE_DIR="$HOME/mgatk_v_mquad_A10"

  # Barcodes (epiAneufinder-selected cells)
  EPI_BARCODES="$HOME/mgatk_v_mquad_A10/outs/A10/epi_results/epiAneufinder_results/A10_cells_for_mquad.tsv"

  # Depth-passing cells from process_mgatk
  DEPTH_BARCODES="$HOME/mgatk_v_mquad_A10/outs/A10/process_mgatk/cells_passing_depth.tsv"

  # mgatk filtered variants table (to restrict mtSNPs)
  VARFILT_TSV="$HOME/mgatk_v_mquad_A10/outs/A10/process_mgatk/variant_stats_filtered.tsv"

  # Output base directory
  OUT_BASE="${BASE_DIR}/outs"

  THREADS=20
  MT_CONTIG="chrM"
fi

CELL_SNP_OUT="${OUT_BASE}/cellsnp_per_cell_mtSNP"
MQUAD_OUT="${OUT_BASE}/mquad_mt"

mkdir -p "${CELL_SNP_OUT}" "${MQUAD_OUT}"

# sanity checks
command -v cellsnp-lite >/dev/null 2>&1 || {
  echo "ERROR: cellsnp-lite not found in PATH"; exit 1;
}
command -v mquad >/dev/null 2>&1 || {
  echo "ERROR: mquad not found in PATH"; exit 1;
}
command -v bcftools >/dev/null 2>&1 || {
  echo "ERROR: bcftools not found in PATH (needed to filter VCF)"; exit 1;
}
command -v tabix >/dev/null 2>&1 || {
  echo "ERROR: tabix not found in PATH (needed to index VCF)"; exit 1;
}

echo "[INFO] BAM      : ${BAM}"
echo "[INFO] EPI_BARCODES   : ${EPI_BARCODES}"
echo "[INFO] DEPTH_BARCODES : ${DEPTH_BARCODES}"
echo "[INFO] VARFILT_TSV    : ${VARFILT_TSV}"
echo "[INFO] OUT_BASE : ${OUT_BASE}"
echo "[INFO] MT_CONTIG: ${MT_CONTIG}"

############################################
# cellsnp-lite: pileup mtDNA in single cells
############################################
CALL_OUT="${OUT_BASE}/cellsnp_bulk_mtSNP"
mkdir -p "${CALL_OUT}"

echo "starting step 1 (mode 2b)"

cellsnp-lite \
  -s "${BAM}" \
  -O "${CALL_OUT}" \
  -p "${THREADS}" \
  --chrom "${MT_CONTIG}" \
  --cellTAG None \
  --UMItag None \
  --minMAF 0.001 \
  --minCOUNT 10 \
  --gzip \
  > "${CALL_OUT}/cellsnp-lite.call.log" 2>&1

[[ -f "${CALL_OUT}/cellSNP.base.vcf.gz" ]] || { echo "ERROR: Missing ${CALL_OUT}/cellSNP.base.vcf.gz"; exit 1; }

# NEW: ensure step-1 VCF is indexed (required downstream for -R usage patterns)
if [[ ! -f "${CALL_OUT}/cellSNP.base.vcf.gz.tbi" && ! -f "${CALL_OUT}/cellSNP.base.vcf.gz.csi" ]]; then
  echo "[INFO] Indexing step-1 VCF: ${CALL_OUT}/cellSNP.base.vcf.gz"
  tabix -f -p vcf "${CALL_OUT}/cellSNP.base.vcf.gz"
fi

echo "ending step 1 (mode 2b)"

############################################
# NEW: intersect barcodes (epi ∩ depth)
############################################
BARCODES="${OUT_BASE}/barcodes_epi_x_depth.tsv"

# EPI_BARCODES has no header; DEPTH_BARCODES has header "barcode" (from your python script)
comm -12 \
  <(sort -u "${EPI_BARCODES}") \
  <(tail -n +2 "${DEPTH_BARCODES}" | sort -u) \
  > "${BARCODES}"

echo "[INFO] Intersected barcodes written: ${BARCODES}"
echo "[INFO] n_barcodes: $(wc -l < "${BARCODES}")"

############################################
# NEW: restrict mtSNPs to mgatk filtered sites
# variant_stats_filtered.tsv has column 1 = position
############################################
VCF_IN="${CALL_OUT}/cellSNP.base.vcf.gz"
REGIONS_BED="${CALL_OUT}/mgatk_pass_sites.bed"
VCF_KEEP="${CALL_OUT}/cellSNP.base.filtered_to_mgatk.vcf.gz"

awk -v OFS="\t" -v CONTIG="${MT_CONTIG}" 'NR>1 {print CONTIG, $1-1, $1}' "${VARFILT_TSV}" > "${REGIONS_BED}"

bcftools view -R "${REGIONS_BED}" -Oz -o "${VCF_KEEP}" "${VCF_IN}"
tabix -f -p vcf "${VCF_KEEP}"

[[ -f "${VCF_KEEP}" ]] || { echo "ERROR: Missing filtered VCF: ${VCF_KEEP}"; exit 1; }

GENO_OUT="${CELL_SNP_OUT}"
mkdir -p "${GENO_OUT}"

echo "starting step 2 (mode 1a)"

cellsnp-lite \
  -s "${BAM}" \
  -b "${BARCODES}" \
  -O "${GENO_OUT}" \
  -p "${THREADS}" \
  -R "${VCF_KEEP}" \
  --UMItag None \
  --cellTAG CB \
  --chrom "${MT_CONTIG}" \
  --minMAF 0.001 \
  --minCOUNT 10 \
  --gzip \
  > "${GENO_OUT}/cellsnp-lite.geno.log" 2>&1

echo "ending step 2 (mode 1a)"

VCF_GZ="${CELL_SNP_OUT}/cellSNP.base.vcf.gz"
VCF_PLAIN="${CELL_SNP_OUT}/cellSNP.base.vcf"

if [[ ! -f "${VCF_GZ}" ]]; then
  if [[ -f "${VCF_PLAIN}" ]]; then
    echo "[INFO] Found uncompressed VCF. Compressing: ${VCF_PLAIN}"

    if command -v bgzip >/dev/null 2>&1; then
      bgzip -f "${VCF_PLAIN}"
    else
      gzip -f "${VCF_PLAIN}"
    fi

    # Optional: index if possible
    tabix -p vcf "${VCF_GZ}" || true

  else
    echo "ERROR: Neither ${VCF_GZ} nor ${VCF_PLAIN} exists in ${CELL_SNP_OUT}"
    exit 1
  fi
fi

#########
# MQuad #
#########

echo "Running MQuad ..."
mquad \
  -c "${CELL_SNP_OUT}" \
  -o "${MQUAD_OUT}" \
  -p "${THREADS}" \
  --minDP 5 \
  > "${MQUAD_OUT}/mquad.log" 2>&1

echo "MQuad has completed. Woohoo!"

#############
# vireoSNP (Python script, see below)
#############

# We now call one python script that should be placed in the base dir:
#   1) fit_mito_clones.py

# If BASE_DIR isn't defined (because we ran via args mode), infer it from OUT_BASE
if [[ -z "${BASE_DIR:-}" ]]; then
  BASE_DIR="$(cd "$(dirname "${OUT_BASE}")" && pwd)"
fi

VIERO_PY="${BASE_DIR}/fit_mito_clones.py"
VIERO_OUT="${MQUAD_OUT}/vireo"

[[ -f "${VIERO_PY}" ]] || { echo "ERROR: Missing ${VIERO_PY}"; exit 1; }

mkdir -p "${VIERO_OUT}"

echo "Running vireo script in ${MQUAD_OUT} ..."
pushd "${MQUAD_OUT}" >/dev/null

python3 "${VIERO_PY}" > "${VIERO_OUT}/fit_mito_clones.log" 2>&1

echo "fit_mito_clones.py done."

popd >/dev/null
echo "Pipeline finished."
