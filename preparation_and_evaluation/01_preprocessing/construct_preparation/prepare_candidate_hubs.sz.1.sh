#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Check for Command-Line Argument ---
if [ -z "$1" ]; then
  echo "Usage: $0 <FactorID>"
  echo "Example: $0 MYB"
  echo "Example: $0 bHLH"
  exit 1
fi

# Store the factor ID from the first argument
FACTOR_ID="$1"
echo "Running pipeline for Factor ID: ${FACTOR_ID}"
echo # Blank line for readability

# --- Input Files (Constructed using FACTOR_ID) ---
# GWAS VCF file: e.g., GWASc1-5_MYB_snps.vcf
INPUT_VCF="GWASc1-5_${FACTOR_ID}_snps.vcf"
# Reference Genome FASTA: (static)
REF_FASTA="Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
# List of sequence IDs/regions: e.g., MYBc1-5_snp_windows.ids.txt
ID_LIST="${FACTOR_ID}c1-5_snp_windows.ids.txt"
# GFF annotation file: (static)
GFF_FILE="Arabidopsis_thaliana.TAIR10.55.gff3"
# Python scripts: (static paths assumed)
EXTRACT_SEQ_SCRIPT="python extract_seq_ranges_from_list.py"
ADD_SUFFIX_SCRIPT="python add_header_suffix.py"
FIND_ORIENT_SCRIPT="python find_orientation_for_promoter.py"
REV_COMP_SCRIPT="python rev_comp_minus_seq.py"
PREPARE_FASTA_SCRIPT="python prepare_fasta_for_TJscreening.py"
FILTER_OVERLAP_SCRIPT="python filter_fasta_overlaps.py"
REMOVE_IDENTICAL_SCRIPT="python remove_identical_sequences.py" # Assuming script name based on Step 7
KEEP_PAIRS_SCRIPT="python keep_pairs.py" # Assuming script name based on Step 7

# --- Intermediate and Output Filenames (Derived from FACTOR_ID) ---
SORTED_VCF="GWASc1-5_${FACTOR_ID}_snps.0.vcf"
BGZIPPED_VCF="${SORTED_VCF}.gz"
CONSENSUS_FASTA="Atha-GWASc1-5_${FACTOR_ID}_snps.0.fa"
CTRL_BASE="Athal_${FACTOR_ID}c1-5_dCIS-ctrl" # Base name for control files
VAR_BASE="Athal_${FACTOR_ID}c1-5_dCIS-var"   # Base name for variant files
COMBINED_BASE="Athal_${FACTOR_ID}c1-5_dCIS_ctrl-var" # Base name for combined files

# --- Ensure Input Files Exist ---
echo "Checking for required input files..."
if [ ! -f "${INPUT_VCF}" ]; then echo "Error: Input VCF not found: ${INPUT_VCF}"; exit 1; fi
if [ ! -f "${REF_FASTA}" ]; then echo "Error: Reference FASTA not found: ${REF_FASTA}"; exit 1; fi
if [ ! -f "${ID_LIST}" ]; then echo "Error: ID list not found: ${ID_LIST}"; exit 1; fi
if [ ! -f "${GFF_FILE}" ]; then echo "Error: GFF file not found: ${GFF_FILE}"; exit 1; fi
# Add checks for python scripts if they aren't guaranteed to be in PATH
echo "Input files checked."

# --- VCF Processing ---
echo "Initializing GWAS: Sorting VCF file..."
bcftools sort "${INPUT_VCF}" -o "${SORTED_VCF}"

echo "Compressing VCF file..."
bgzip "${SORTED_VCF}"

echo "Indexing compressed VCF file..."
tabix -p vcf "${BGZIPPED_VCF}"

# --- Consensus and Initial Extraction ---
echo "Generating consensus FASTA sequence..."
bcftools consensus -f "${REF_FASTA}" "${BGZIPPED_VCF}" > "${CONSENSUS_FASTA}"

echo "Extracting candidate sequence ranges"
${EXTRACT_SEQ_SCRIPT} "${CONSENSUS_FASTA}" "${ID_LIST}" "${VAR_BASE}.0.fa"
${EXTRACT_SEQ_SCRIPT} "${REF_FASTA}" "${ID_LIST}" "${CTRL_BASE}.0.fa"
echo "done..."

# --- Sequence Preparation ---
echo "Start preparing sequences for submission:"
echo "Step 1: Adding header suffix..."
${ADD_SUFFIX_SCRIPT} "${CTRL_BASE}.0.fa" "${CTRL_BASE}.1.fa"
${ADD_SUFFIX_SCRIPT} "${VAR_BASE}.0.fa" "${VAR_BASE}.1.fa"

echo "Step 2: Finding orientation..."
${FIND_ORIENT_SCRIPT} "${GFF_FILE}" "${CTRL_BASE}.1.fa" -d 2000 -o "${CTRL_BASE}.2.fa"
${FIND_ORIENT_SCRIPT} "${GFF_FILE}" "${VAR_BASE}.1.fa" -d 2000 -o "${VAR_BASE}.2.fa"

echo "Step 3: Reverse complementing minus-strand sequences..."
${REV_COMP_SCRIPT} "${CTRL_BASE}.2.fa" -o "${CTRL_BASE}.3.fa"
${REV_COMP_SCRIPT} "${VAR_BASE}.2.fa" -o "${VAR_BASE}.3.fa"

echo "Step 4: Preparing FASTA for TJscreening..."
${PREPARE_FASTA_SCRIPT} "${CTRL_BASE}.3.fa" "${CTRL_BASE}.4.fa"
${PREPARE_FASTA_SCRIPT} "${VAR_BASE}.3.fa" "${VAR_BASE}.4.fa"

# --- Filtering and Joining ---
echo "Step 5: Remove overlapping sequences"
${FILTER_OVERLAP_SCRIPT} "${CTRL_BASE}.4.fa" "${CTRL_BASE}.5.fa" -t 40
${FILTER_OVERLAP_SCRIPT} "${VAR_BASE}.4.fa" "${VAR_BASE}.5.fa" -t 40

echo "Step 6: Join control and variant multi-fasta"
cat "${CTRL_BASE}.5.fa" "${VAR_BASE}.5.fa" > "${COMBINED_BASE}.0.fa"
echo "Combined file created: ${COMBINED_BASE}.0.fa"

# --- Final Filtering Steps ---
echo "Step 7: Remove identical sequences and keep pairs"
${REMOVE_IDENTICAL_SCRIPT} "${COMBINED_BASE}.0.fa" "${COMBINED_BASE}.1.fa"
${KEEP_PAIRS_SCRIPT} "${COMBINED_BASE}.1.fa" "${COMBINED_BASE}.final.fa"
echo "Final filtered pairs file created: ${COMBINED_BASE}.final.fa"


# --- Cleanup ---
echo "Step 8: Cleaning up intermediate files..."
# This list is based on the user's original script, parameterizing the names.
# Note: Other intermediates like .3.fa, VCFs, Consensus might still exist.
rm -f \
    "${SORTED_VCF}" \
    "${BGZIPPED_VCF}" \
    "${BGZIPPED_VCF}.tbi" \
    "${CONSENSUS_FASTA}" \
    "${CTRL_BASE}.0.fa" \
    "${CTRL_BASE}.1.fa" \
    "${CTRL_BASE}.2.fa" \
    "${CTRL_BASE}.3.fa" \
    "${CTRL_BASE}.4.fa" \
    "${VAR_BASE}.0.fa" \
    "${VAR_BASE}.1.fa" \
    "${VAR_BASE}.2.fa" \
    "${VAR_BASE}.3.fa" \
    "${VAR_BASE}.4.fa" \
    "${COMBINED_BASE}.0.fa" \
    "${COMBINED_BASE}.1.fa"
    # Note: Keeping .5.fa files as they were explicitly kept in the user's prompt context previously,
    # although they are technically intermediates before the final .final.fa file. Remove if not needed.
    # Also keeping original input VCF and ID list.

echo "Cleanup complete. Kept files include:"
echo "- ${CTRL_BASE}.5.fa" # Kept based on previous context
echo "- ${VAR_BASE}.5.fa"   # Kept based on previous context
echo "- ${COMBINED_BASE}.final.fa" # The final output
echo "(And original input files like ${INPUT_VCF}, ${ID_LIST} etc.)"

echo "Pipeline finished successfully for Factor ID: ${FACTOR_ID}."

# --- End of Script ---