#!/usr/bin/env bash
set -euo pipefail

# Args:
#   1 BAM
#   2 EPI_BARCODES (no header)
#   3 DEPTH_BARCODES (has header "barcode")
#   4 VARFILT_TSV (variant_stats_filtered.tsv; col1=position with header)
#   5 OUT_BASE
#   6 THREADS
#   7 MT_CONTIG

BAM="$1"
EPI_BARCODES="$2"
DEPTH_BARCODES="$3"
VARFILT_TSV="$4"
OUT_BASE="$5"
THREADS="${6:-20}"
MT_CONTIG="${7:-chrM}"

CELL_SNP_OUT="${OUT_BASE}/cellsnp_per_cell_mtSNP"
CALL_OUT="${OUT_BASE}/cellsnp_bulk_mtSNP"

mkdir -p "${CELL_SNP_OUT}" "${CALL_OUT}"

command -v cellsnp-lite >/dev/null 2>&1 || { echo "ERROR: cellsnp-lite not found in PATH"; exit 1; }
command -v bcftools   >/dev/null 2>&1 || { echo "ERROR: bcftools not found in PATH"; exit 1; }
command -v tabix      >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH"; exit 1; }

echo "[INFO] BAM      : ${BAM}"
echo "[INFO] EPI_BARCODES   : ${EPI_BARCODES}"
echo "[INFO] DEPTH_BARCODES : ${DEPTH_BARCODES}"
echo "[INFO] VARFILT_TSV    : ${VARFILT_TSV}"
echo "[INFO] OUT_BASE : ${OUT_BASE}"
echo "[INFO] MT_CONTIG: ${MT_CONTIG}"
echo "[INFO] THREADS  : ${THREADS}"

############################################
# Step 1: cellsnp-lite mode 2b (site discovery)
############################################
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

# Ensure index exists
if [[ ! -f "${CALL_OUT}/cellSNP.base.vcf.gz.tbi" && ! -f "${CALL_OUT}/cellSNP.base.vcf.gz.csi" ]]; then
  echo "[INFO] Indexing step-1 VCF: ${CALL_OUT}/cellSNP.base.vcf.gz"
  tabix -f -p vcf "${CALL_OUT}/cellSNP.base.vcf.gz"
fi

echo "ending step 1 (mode 2b)"

############################################
# Step 1.5: intersect barcodes (epi ∩ depth)
############################################
BARCODES="${OUT_BASE}/barcodes_epi_x_depth.tsv"

comm -12 \
  <(sort -u "${EPI_BARCODES}") \
  <(tail -n +2 "${DEPTH_BARCODES}" | sort -u) \
  > "${BARCODES}"

echo "[INFO] Intersected barcodes written: ${BARCODES}"
echo "[INFO] n_barcodes: $(wc -l < "${BARCODES}")"

############################################
# Step 1.75: restrict VCF sites to mgatk filtered positions
############################################
VCF_IN="${CALL_OUT}/cellSNP.base.vcf.gz"
REGIONS_BED="${CALL_OUT}/mgatk_pass_sites.bed"
VCF_KEEP="${CALL_OUT}/cellSNP.base.filtered_to_mgatk.vcf.gz"

awk -v OFS="\t" -v CONTIG="${MT_CONTIG}" 'NR>1 {print CONTIG, $1-1, $1}' "${VARFILT_TSV}" > "${REGIONS_BED}"

bcftools view -R "${REGIONS_BED}" -Oz -o "${VCF_KEEP}" "${VCF_IN}"
tabix -f -p vcf "${VCF_KEEP}"

[[ -f "${VCF_KEEP}" ]] || { echo "ERROR: Missing filtered VCF: ${VCF_KEEP}"; exit 1; }

############################################
# Step 2: cellsnp-lite mode 1a (per-cell genotyping)
############################################
echo "starting step 2 (mode 1a)"

cellsnp-lite \
  -s "${BAM}" \
  -b "${BARCODES}" \
  -O "${CELL_SNP_OUT}" \
  -p "${THREADS}" \
  -R "${VCF_KEEP}" \
  --UMItag None \
  --cellTAG CB \
  --chrom "${MT_CONTIG}" \
  --minMAF 0.001 \
  --minCOUNT 10 \
  --gzip \
  > "${CELL_SNP_OUT}/cellsnp-lite.geno.log" 2>&1

echo "ending step 2 (mode 1a)"

# Ensure genotype VCF exists and is indexed (some versions produce it)
if [[ -f "${CELL_SNP_OUT}/cellSNP.base.vcf.gz" ]]; then
  if [[ ! -f "${CELL_SNP_OUT}/cellSNP.base.vcf.gz.tbi" && ! -f "${CELL_SNP_OUT}/cellSNP.base.vcf.gz.csi" ]]; then
    tabix -f -p vcf "${CELL_SNP_OUT}/cellSNP.base.vcf.gz" || true
  fi
fi

echo "[DONE] cellsnp pipeline finished"
