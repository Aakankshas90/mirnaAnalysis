#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Define paths
case_dir="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/case/"
control_dir="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/control/"
ref_genome="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/hg38_reference/hg38.fa"
mature_rna="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/mature_human.fa"
hairpin_rna="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/hairpin_human.fa"
mirna_gff="/Users/balgovindyadav/Downloads/Aakanksha/miRNAmetaAnalysis/miRNA.gff3"
log_file="mirdeep2_pipeline.log"

# Function to run miRDeep2 analysis
run_mirdeep2() {
    local dataset_dir=$1
    local collapsed_fa=$2
    local arf_file=$3

    # Extract dataset identifier
    dataset_id=$(basename "$collapsed_fa" | cut -d'_' -f1)
    
    # Create result directory
    result_dir="${dataset_dir}/${dataset_id}_mirdeep2_results"
    mkdir -p "$result_dir"
    
    echo "Processing dataset $dataset_id in $dataset_dir" | tee -a "$log_file"

    # Run the miRDeep2 analysis
    miRDeep2.pl "$collapsed_fa" "$ref_genome" "$arf_file" "$mature_rna" none "$hairpin_rna" "$mirna_gff" -t hsa -v | tee -a "$log_file"

    # Move output files to the result directory
    mv -v ./expression_analyses* ./pdfs ./result* ./tmp "${result_dir}/" | tee -a "$log_file"
    echo "Results moved to ${result_dir}" | tee -a "$log_file"
}

# Process case and control directories
for parent_dir in "$case_dir" "$control_dir"; do
    # Iterate over all subdirectories
    for dataset_dir in "$parent_dir"*/; do
        echo "Entering directory $dataset_dir" | tee -a "$log_file"

        # Check for collapsed.fa and .arf files
        for collapsed_fa in "${dataset_dir}"*_collapsed.fa; do
            # Ensure corresponding ARF file exists
            arf_file="${collapsed_fa/_collapsed.fa/_reads_vs_ref.arf}"
            if [[ -f "$collapsed_fa" && -f "$arf_file" ]]; then
                run_mirdeep2 "$dataset_dir" "$collapsed_fa" "$arf_file"
            else
                echo "Error: Missing required files for dataset in $dataset_dir" | tee -a "$log_file"
            fi
        done
    done
done

echo "miRDeep2 analysis complete for all datasets." | tee -a "$log_file"
