#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Paths to case and control datasets
case_dir="/Volumes/New Volume 1/mirnaAnalysis/case"
control_dir="/Volumes/New Volume 1/mirnaAnalysis/control"
bowtie_index="/Volumes/New Volume 1/mirnaAnalysis/hg38_reference"
ref_genome="/Volumes/New Volume 1/mirnaAnalysis/hg38_reference/hg38.fa"
miRbase_file="/Volumes/New Volume 1/mirnaAnalysis/miRNA.gff3"
mature_rna="/Volumes/New Volume 1/mirnaAnalysis/mature_human.fa"
hairpin_rna="/Volumes/New Volume 1/mirnaAnalysis/hairpin_human.fa"
mirna_gff="/Volumes/New Volume 1/mirnaAnalysis/miRNA.gff3"
file_extension=".fastq.gz"  # Adjust if using different extension

# Define log file
log_file="pipeline.log"

# Clear the log file if it exists
true > "$log_file"

# Check if required tools are in PATH
for tool in fastqc fastp mapper.pl miRDeep2.pl; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool is not installed or not in PATH." | tee -a "$log_file" >&2
        exit 1
    fi
done

# Function to run the pipeline on each dataset
process_dataset() {
    local input_dir=$1

    for subfolder in "$input_dir"/*/; do
        if [ -d "$subfolder" ]; then
            local results_dir="${subfolder}/results"
            mkdir -p "$results_dir"
            local mirdeep_out="${results_dir}/mirdeep_out"
            mkdir -p "$mirdeep_out"
            echo "Processing subfolder: $subfolder" | tee -a "$log_file"

            # Quality Control using FastQC and MultiQC
            echo "Running FastQC on $subfolder" | tee -a "$log_file"
            fastqc -o "$results_dir" "${subfolder}"/*"${file_extension}" >> "$log_file" 2>&1 || {
                echo "FastQC failed for $subfolder" | tee -a "$log_file"
                continue
            }

            # Run adapter trimming and alignment in parallel
            for fastq_file in "${subfolder}"/*"${file_extension}"; do
                if [ -e "$fastq_file" ]; then
                    local sample_name
                    sample_name=$(basename "$fastq_file" "$file_extension")
                    local trimmed_file="${results_dir}/${sample_name}_trimmed.fastq.gz"
                    local arf_file="${results_dir}/${sample_name}.arf"
                    local mapped_file="${results_dir}/${sample_name}.fa"

                    # Adapter trimming using fastp
                    (
                        echo "Running fastp on $fastq_file" | tee -a "$log_file"
                        fastp -i "$fastq_file" -o "$trimmed_file" >> "$log_file" 2>&1 || {
                            echo "fastp failed for $fastq_file" | tee -a "$log_file"
                            exit 1
                        }

                        # miRNA alignment and quantification using miRDeep2
                        echo "Running miRDeep2 (mapper.pl) on $trimmed_file" | tee -a "$log_file"
                        mapper.pl "$trimmed_file" -e -h -i -j -m -p "$bowtie_index" -s "$mapped_file" -t "$arf_file" >> "$log_file" 2>&1 || {
                            echo "miRDeep2 mapper failed for $trimmed_file" | tee -a "$log_file"
                            exit 1
                        }

                        # Run miRDeep2 for miRNA prediction and quantification
                        echo "Running miRDeep2 on $arf_file" | tee -a "$log_file"
                        miRDeep2.pl "$mapped_file" "$ref_genome" "$miRbase_file" "$mature_rna" "$hairpin_rna" "$arf_file" "$mirna_gff" "${mirdeep_out}" >> "$log_file" 2>&1 || {
                            echo "miRDeep2 failed for $arf_file" | tee -a "$log_file"
                            exit 1
                        }

                        echo "Successfully processed $sample_name" | tee -a "$log_file"
                    ) & # Run in background

                else
                    echo "No files found with extension $file_extension in $subfolder" | tee -a "$log_file"
                fi
            done

            # Wait for all background processes to complete
            wait
        else
            echo "No subfolders found in $input_dir" | tee -a "$log_file"
        fi
    done
}

# Data Preprocessing for Case and Control datasets
echo "Processing case datasets..." | tee -a "$log_file"
process_dataset "$case_dir"

echo "Processing control datasets..." | tee -a "$log_file"
process_dataset "$control_dir"
