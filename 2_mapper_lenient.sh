#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# miRNA mapping with relaxed parameters (for comparison)

case_dir="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/case"
control_dir="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/control"
bowtie_index="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/hg38_reference/hg38_index"

file_extension=".fastq.gz"
log_file="mirna_pipeline.log"

echo "Starting miRNA mapping (lenient mode)" | tee -a "$log_file"

log_failure() {
    local step="$1"
    local input_file="$2"
    echo "Error: $step failed for $input_file" | tee -a "$log_file"
}

process_file() {
    local input_file="$1"
    local subfolder="$2"

    local results_qc="$subfolder/results_qc"
    local results_mapper="$subfolder/results_mapper_lenient"

    local file_base
    file_base=$(basename "$input_file" "$file_extension")

    local fastp_output="$results_qc/${file_base}_trimmed.fastq.gz"

    mkdir -p "$results_mapper"

    gunzip "$fastp_output"
    fastp_uncompressed="${fastp_output%.gz}"

    mapper.pl "$fastp_uncompressed" -e -h -i -m -j -l 18 -v -q -r 10 -o 8 \
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

echo "Pipeline completed (lenient mode)" | tee -a "$log_file"
