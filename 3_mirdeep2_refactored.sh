#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# miRDeep2 analysis for miRNA identification and quantification

# directories (consistent with previous scripts)
case_dir="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/case"
control_dir="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/control"

ref_genome="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/hg38_reference/hg38.fa"
mature_rna="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/mature_human.fa"
hairpin_rna="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/hairpin_human.fa"
mirna_gff="/RAID1/working/A423/aakanksha/miRNAmetaAnalysis/miRNA.gff3"

log_file="mirdeep2_pipeline.log"

echo "Starting miRDeep2 analysis" | tee -a "$log_file"

run_mirdeep2() {
    local dataset_dir="$1"
    local collapsed_fa="$2"
    local arf_file="$3"

    # safer sample name extraction
    local sample_name
    sample_name=$(basename "$collapsed_fa" "_collapsed.fa")

    local result_dir="${dataset_dir}/${sample_name}_mirdeep2_results"
    local temp_dir="${dataset_dir}/${sample_name}_mirdeep2_tmp"

    mkdir -p "$result_dir" "$temp_dir"

    echo "Processing sample: $sample_name" | tee -a "$log_file"

    # Run miRDeep2 inside temp directory (prevents overwrite issues)
    (
        cd "$temp_dir"

        miRDeep2.pl "$collapsed_fa" "$ref_genome" "$arf_file" \
            "$mature_rna" none "$hairpin_rna" "$mirna_gff" \
            -t hsa -v >> "$log_file" 2>&1

        # Move outputs to result directory
        mv -v expression_analyses* pdfs result* tmp "$result_dir/" 2>> "$log_file"
    )

    echo "Finished sample: $sample_name" | tee -a "$log_file"
}

# Process both case and control datasets
for parent_dir in "$case_dir" "$control_dir"; do
    for dataset_dir in "$parent_dir"/*/; do
        echo "Entering directory: $dataset_dir" | tee -a "$log_file"

        for collapsed_fa in "${dataset_dir}"*_collapsed.fa; do
            arf_file="${collapsed_fa/_collapsed.fa/_reads_vs_ref.arf}"

            if [[ -f "$collapsed_fa" && -f "$arf_file" ]]; then
                run_mirdeep2 "$dataset_dir" "$collapsed_fa" "$arf_file"
            else
                echo "Missing files for $collapsed_fa" | tee -a "$log_file"
            fi
        done
    done
done

echo "miRDeep2 analysis completed" | tee -a "$log_file"
