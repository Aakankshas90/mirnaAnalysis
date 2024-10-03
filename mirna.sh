#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Paths to directories and files
case_dir="/mnt/d/setolabo/case"
control_dir="/mnt/d/setolabo/control"
bowtie_index="/mnt/d/setolabo/hg38_reference/hg38_index"
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
    local results_mirdeep2="$subfolder/results_mapper"
    local file_base
    file_base=$(basename "$input_file" "$file_extension")
    local trimmomatic_output="$results_qc/${file_base}_trimmed.fastq.gz"
    
    mkdir -p "$results_qc" "$results_mirdeep2"
    
    # Step 1: FastQC
    local fastqc_output="$results_qc/${file_base}_fastqc.zip"
    if [[ ! -f "$fastqc_output" ]]; then
        echo "Running FastQC on $input_file" | tee -a "$log_file"
        fastqc -o "$results_qc" "$input_file" -t 4 >> "$log_file" 2>&1 || { log_failure "FastQC" "$input_file"; return 1; }
    else
        echo "FastQC already completed for $input_file" | tee -a "$log_file"
    fi
    
    # Step 2: trimmomatic
    if [[ ! -f "$trimmomatic_output" ]]; then
        echo "Running trimmomatic on $input_file" | tee -a "$log_file"
    trimmomatic SE "$input_file" "$trimmomatic_output" ILLUMINACLIP:Illumina_small_RNA_adapter.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16 >> "$log_file" 2>&1 || { log_failure "trimmomatic" "$input_file"; return 1; }
    else
        echo "trimmomatic already completed for $input_file" | tee -a "$log_file"
    fi

    # Uncompress trimmomatic output
    gunzip "$trimmomatic_output" || { log_failure "Unzipping trimmomatic output" "$input_file"; return 1; }
    trimmomatic_uncompressed="${trimmomatic_output%.gz}"

    # Step 3: mapper.pl
    local mapper_output="$results_mirdeep2/${file_base}_collapsed.fa"
    if [[ ! -f "$mapper_output" ]]; then
        echo "Running mapper.pl on $trimmomatic_uncompressed" | tee -a "$log_file"
        mapper.pl "$trimmomatic_uncompressed" -e -h -i -m -j -l 18 -v -q -r 10 -o 4 \
              -s "$mapper_output" \
              -t "$results_mirdeep2/${file_base}_reads_vs_ref.arf" -p "$bowtie_index" >> "$log_file" 2>&1 || { log_failure "mapper.pl" "$input_file"; return 1; }
    else
        echo "mapper.pl already completed for $input_file" | tee -a "$log_file"
    fi
    
    gzip "$trimmomatic_uncompressed" || { log_failure "Recompressing trimmomatic output" "$input_file"; return 1; }
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
