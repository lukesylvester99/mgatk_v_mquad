#!/usr/bin/env bash
set -euo pipefail

# Args:
#   1 BAM
#   2 BARCODES (no header; epi ∩ depth list)
#   3 VCF_MGATK (cellSNP.base.filtered_to_mgatk.vcf.gz)
#   4 MQUAD_PASSED (passed_variant_names.txt; lines like 14783T>C)
#   5 CLONE_ASSIGN (mquad vireo/clone_assignments.tsv)
#   6 OUTDIR (workdir for trunk analysis; will contain vireo/ outputs)
#   7 VIREO_TRUNK_PY (path to vireo_trunk.py)
#   8 THREADS (optional; default 8)
#   9 MT_CONTIG (optional; default chrM)

BAM="$1"
BARCODES="$2"
VCF_MGATK="$3"
MQUAD_PASSED="$4"
CLONE_ASSIGN="$5"
OUTDIR="$6"
VIREO_TRUNK_PY="$7"
THREADS="${8:-8}"
MT_CONTIG="${9:-chrM}"

command -v cellsnp-lite >/dev/null 2>&1 || { echo "ERROR: cellsnp-lite not found in PATH"; exit 1; }
command -v bcftools    >/dev/null 2>&1 || { echo "ERROR: bcftools not found in PATH"; exit 1; }
command -v tabix       >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH"; exit 1; }

mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/vireo"

# Where we’ll write the new per-cell matrices
CELL_SNP_TRUNK_OUT="${OUTDIR}/cellsnp_per_cell_mtSNP_trunk"
mkdir -p "${CELL_SNP_TRUNK_OUT}"

echo "==== vireo_trunk.sh START $(date) ===="
echo "[INFO] BAM         : ${BAM}"
echo "[INFO] BARCODES    : ${BARCODES}"
echo "[INFO] VCF_MGATK   : ${VCF_MGATK}"
echo "[INFO] MQUAD_PASSED: ${MQUAD_PASSED}"
echo "[INFO] CLONE_ASSIGN: ${CLONE_ASSIGN}"
echo "[INFO] OUTDIR      : ${OUTDIR}"
echo "[INFO] THREADS     : ${THREADS}"
echo "[INFO] MT_CONTIG   : ${MT_CONTIG}"
echo

############################################
# Step 1: Copy mgatk-filtered VCF to OUTDIR
############################################
VCF_LOCAL="${OUTDIR}/cellSNP.base.filtered_to_mgatk.vcf.gz"
VCF_LOCAL_TBI="${VCF_LOCAL}.tbi"

cp -f "${VCF_MGATK}" "${VCF_LOCAL}"
if [[ -f "${VCF_MGATK}.tbi" ]]; then
  cp -f "${VCF_MGATK}.tbi" "${VCF_LOCAL_TBI}"
elif [[ -f "${VCF_MGATK}.csi" ]]; then
  cp -f "${VCF_MGATK}.csi" "${VCF_LOCAL}.csi"
else
  echo "[INFO] Indexing copied VCF: ${VCF_LOCAL}"
  tabix -f -p vcf "${VCF_LOCAL}" || true
fi

[[ -s "${VCF_LOCAL}" ]] || { echo "ERROR: Missing/empty copied VCF: ${VCF_LOCAL}"; exit 2; }

############################################################
# Step 2: Remove mquad-passed variants -> new trunk VCF
# (remove by chrM position; passed list is like 14783T>C)
############################################################
EXCLUDE_BED="${OUTDIR}/mquad_pass_sites.bed"
VCF_TRUNK="${OUTDIR}/cellSNP.base.filtered_to_mgatk_and_mquad.vcf.gz"

# Build BED: chrM  (pos-1)  pos
# Example line: 14783T>C  => pos=14783
awk -v OFS="\t" -v CONTIG="${MT_CONTIG}" '
  {
    gsub(/\r/, "", $1)
  }
  NF>=1 && $1 ~ /^[0-9]+[ACGT]>[ACGT]$/ {
    pos = $1
    sub(/[^0-9].*$/, "", pos)
    if (pos ~ /^[0-9]+$/) {
      p = pos + 0
      if (p > 0) print CONTIG, p-1, p
    }
  }
' "${MQUAD_PASSED}" | sort -u > "${EXCLUDE_BED}"

echo "[INFO] Exclusion BED written: ${EXCLUDE_BED}"
echo "[INFO] n_exclude_sites: $(wc -l < "${EXCLUDE_BED}")"

# If exclusion bed is empty, just copy VCF forward
if [[ "$(wc -l < "${EXCLUDE_BED}")" -eq 0 ]]; then
  echo "[WARN] Exclusion BED is empty; no sites removed. Copying VCF forward."
  cp -f "${VCF_LOCAL}" "${VCF_TRUNK}"
  if [[ -f "${VCF_LOCAL}.tbi" ]]; then
    cp -f "${VCF_LOCAL}.tbi" "${VCF_TRUNK}.tbi"
  else
    tabix -f -p vcf "${VCF_TRUNK}" || true
  fi
else
  # Exclude positions using bcftools -T ^BED
  bcftools view -T ^"${EXCLUDE_BED}" -Oz -o "${VCF_TRUNK}" "${VCF_LOCAL}"
  tabix -f -p vcf "${VCF_TRUNK}"
fi

[[ -s "${VCF_TRUNK}" ]] || { echo "ERROR: Missing/empty trunk VCF: ${VCF_TRUNK}"; exit 2; }

############################################
# Step 3: cellsnp-lite mode 1a (per-cell genotyping) on trunk VCF
############################################
echo
echo "starting trunk step (mode 1a genotyping)"

cellsnp-lite \
  -s "${BAM}" \
  -b "${BARCODES}" \
  -O "${CELL_SNP_TRUNK_OUT}" \
  -p "${THREADS}" \
  -R "${VCF_TRUNK}" \
  --UMItag None \
  --cellTAG CB \
  --chrom "${MT_CONTIG}" \
  --minMAF 0.001 \
  --minCOUNT 10 \
  --gzip \
  > "${CELL_SNP_TRUNK_OUT}/cellsnp-lite.trunk.geno.log" 2>&1

echo "ending trunk step (mode 1a genotyping)"

# Require expected matrices
[[ -s "${CELL_SNP_TRUNK_OUT}/cellSNP.tag.AD.mtx" ]] || { echo "ERROR: missing AD mtx"; exit 2; }
[[ -s "${CELL_SNP_TRUNK_OUT}/cellSNP.tag.DP.mtx" ]] || { echo "ERROR: missing DP mtx"; exit 2; }
[[ -s "${CELL_SNP_TRUNK_OUT}/cellSNP.samples.tsv" ]] || { echo "ERROR: missing samples.tsv"; exit 2; }

# Some cellSNP-lite builds do NOT emit cellSNP.variants.tsv.
# But mode 1a typically emits a genotyped VCF with the variants used.
GENO_VCF="${CELL_SNP_TRUNK_OUT}/cellSNP.base.vcf.gz"
[[ -s "${GENO_VCF}" ]] || {
  echo "ERROR: missing ${GENO_VCF} (needed to derive variant list)"; 
  ls -lh "${CELL_SNP_TRUNK_OUT}" || true
  exit 2
}

# Create variants.tsv in the expected format: CHROM POS REF ALT
VARIANTS_TSV="${CELL_SNP_TRUNK_OUT}/cellSNP.variants.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${GENO_VCF}" > "${VARIANTS_TSV}"
[[ -s "${VARIANTS_TSV}" ]] || { echo "ERROR: failed to create variants.tsv"; exit 2; }

############################################
# Step 4: Stage expected filenames + run vireo_trunk.py
############################################
echo
echo "staging inputs for vireo_trunk.py"

cd "${OUTDIR}"

ln -sf "${CELL_SNP_TRUNK_OUT}/cellSNP.tag.AD.mtx" "cellSNP.tag.AD.mtx"
ln -sf "${CELL_SNP_TRUNK_OUT}/cellSNP.tag.DP.mtx" "cellSNP.tag.DP.mtx"
ln -sf "${CELL_SNP_TRUNK_OUT}/cellSNP.samples.tsv" "cellSNP.samples.tsv"
ln -sf "${CELL_SNP_TRUNK_OUT}/cellSNP.variants.tsv" "cellSNP.variants.tsv"
ln -sf "${CLONE_ASSIGN}" "clone_assignments.tsv"

# Keep this around for provenance (optional, but helps debugging)
ln -sf "${MQUAD_PASSED}" "mquad_passed_variant_names.txt"

echo
echo "running vireo_trunk.py"
python -u "${VIREO_TRUNK_PY}"

echo
echo "==== vireo_trunk.sh END $(date) ===="
