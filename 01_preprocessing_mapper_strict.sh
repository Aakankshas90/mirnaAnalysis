#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# miRNA-seq preprocessing and mapping (strict parameters)
# Includes FastQC, fastp, and miRDeep2 mapper.pl

# directories
case_dir="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/case"
control_dir="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/control"
bowtie_index="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/hg38_reference/hg38_index"

file_extension=".fastq.gz"
log_file="mirna_pipeline.log"

echo "Starting miRNA-seq pipeline (strict mode)" | tee -a "$log_file"

log_failure() {
    local step="$1"
    local input_file="$2"
    echo "Error: $step failed for $input_file" | tee -a "$log_file"
}

process_file() {
    local input_file="$1"
    local subfolder="$2"

    local results_qc="$subfolder/results_qc"
    local results_mapper="$subfolder/results_mapper_strict"

    local file_base
    file_base=$(basename "$input_file" "$file_extension")

    local fastp_output="$results_qc/${file_base}_trimmed.fastq.gz"

    mkdir -p "$results_qc" "$results_mapper"

    # FastQC
    fastqc -o "$results_qc" "$input_file" >> "$log_file" 2>&1 || { log_failure "FastQC" "$input_file"; return 1; }

    # fastp
    fastp -i "$input_file" -o "$fastp_output" \
        -h "$results_qc/${file_base}_fastp.html" \
        -j "$results_qc/${file_base}_fastp.json" >> "$log_file" 2>&1 || { log_failure "fastp" "$input_file"; return 1; }

    gunzip "$fastp_output"
    fastp_uncompressed="${fastp_output%.gz}"

    # mapper.pl (strict)
    mapper.pl "$fastp_uncompressed" -e -h -i -m -j -l 18 -v \
        -s "$results_mapper/${file_base}_collapsed.fa" \
        -t "$results_mapper/${file_base}_reads_vs_ref.arf" \
        -p "$bowtie_index" >> "$log_file" 2>&1 || { log_failure "mapper.pl" "$input_file"; return 1; }

    gzip "$fastp_uncompressed"
}

process_directory() {
    local dir="$1"
    for subfolder in "$dir"/*/; do
        for file in "$subfolder"/*"$file_extension"; do
            process_file "$file" "$subfolder"
        done
    done
}

process_directory "$case_dir"
process_directory "$control_dir"

echo "Pipeline completed (strict mode)" | tee -a "$log_file"
