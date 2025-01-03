#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Paths to directories and files
case_dir="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/case/"
control_dir="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/control"
bowtie_index="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/hg38_reference/hg38_index"
file_extension=".fastq.gz"
log_file="mirna_pipeline.log"

# Log file creation
echo "Starting miRNA-seq pipeline" | tee -a "$log_file"

# Function to log the failure message
log_failure() {
    local step="$1"
    local input_file="$2"
    echo "Error: $step failed for $input_file" | tee -a "$log_file"
}

# Function to process a single file
process_file() {
    local input_file="$1"
    local subfolder="$2"
    local results_qc="$subfolder/results_qc"
    local results_mirdeep2="$subfolder/results_mapper_linient"
    local file_base
    file_base=$(basename "$input_file" "$file_extension")
    local fastp_output="$results_qc/${file_base}_trimmed.fastq.gz"
    
    mkdir -p "$results_qc" "$results_mirdeep2"

    # Uncompress fastp output
    gunzip "$fastp_output" || { log_failure "Unzipping fastp output" "$input_file"; return 1; }
    fastp_uncompressed="${fastp_output%.gz}"

    # Step 3: mapper.pl
    local mapper_output="$results_mirdeep2/${file_base}_collapsed.fa"
    if [[ ! -f "$mapper_output" ]]; then
        echo "Running mapper.pl on $fastp_uncompressed" | tee -a "$log_file"
        mapper.pl "$fastp_uncompressed" -e -h -i -m -j -l 18 -v -q -r 10 -o 8 \
              -s "$mapper_output" \
              -t "$results_mirdeep2/${file_base}_reads_vs_ref.arf" -p "$bowtie_index" >> "$log_file" 2>&1 || { log_failure "mapper.pl" "$input_file"; return 1; }
    else
        echo "mapper.pl already completed for $input_file" | tee -a "$log_file"
    fi
    
    gzip "$fastp_uncompressed" || { log_failure "Recompressing fastp output" "$input_file"; return 1; }
}

# Function to process all subfolders in a directory (case or control)
process_directory() {
    local dir="$1"
    
    for subfolder in "$dir"/*/; do
        echo "Processing subfolder $subfolder" | tee -a "$log_file"
        for file in "$subfolder"/*"$file_extension"; do
            echo "Processing file $file in $subfolder" | tee -a "$log_file"
            process_file "$file" "$subfolder" || echo "Pipeline failed for $file." | tee -a "$log_file"
        done
    done
}

# Start processing case and control directories
process_directory "$case_dir"
process_directory "$control_dir"

echo "Pipeline completed" | tee -a "$log_file"
