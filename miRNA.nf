// test workflow for miRNA analysis
#!/usr/bin/env nextflow

// Define workflow parameters with meaningful defaults
params.input = params.input ?: './data/**/*.fastq.gz'
params.outdir = params.outdir ?: './results'
params.index = params.index ?: './reference/genome_index' // Path to the Bowtie2 index
params.gtf = params.gtf ?: './reference/annotation.gtf'   // Path to the GTF annotation file

// Workflow definition
workflow {
    // Create a channel from the input files, pairing them if necessary
    Channel
        .fromFilePairs(params.input, size: -1)
        .map { file, dir -> tuple(file, dir, file.nameWithoutExtension) }
        .set { input_file_info }

    // Run FastQC on the input files
    qc_results = input_file_info
        .process(runFastQC)

    // Align reads using Bowtie2
    bam_files = qc_results
        .process(runBowtie2)

    // Perform differential expression analysis (placeholder process)
    runDifferentialExpression(bam_files)

    // Aggregate quality control results (optional)
    aggregateQCResults(qc_results)
}

// FastQC Process
process runFastQC {
    input:
    tuple path(file), path(dir), val(baseName)

    output:
    path("${baseName}_qc_results")

    script:
    """
    mkdir -p ${dir}/${baseName}_qc_results
    fastqc ${file} -o ${dir}/${baseName}_qc_results &> ${dir}/${baseName}_qc_results/fastqc.log

    if [ $? -ne 0 ]; then
        echo "Error: FastQC failed for ${file}" >&2
        exit 1
    fi
    """
}

// Bowtie2 Alignment Process
process runBowtie2 {
    input:
    tuple path(file), path(dir), val(baseName)

    output:
    path("${baseName}_aligned.bam")

    script:
    """
    bowtie2 -x ${params.index} -U ${file} | samtools view -Sb - > ${dir}/${baseName}_aligned.bam

    if [ $? -ne 0 ]; then
        echo "Error: Bowtie2 alignment failed for ${file}" >&2
        exit 1
    fi
    """
}

// Differential Expression Process with DESeq2
process runDifferentialExpression {
    input:
    path(bam_files)

    output:
    path("${bam_files}_DE_results")

    script:
    """
    mkdir -p ${bam_files}_DE_results

    // Simulate a differential expression analysis (this is the placeholder part)
    Rscript --vanilla -e "
    library(DESeq2)
    // Add your DESeq2 code here to read count data from the BAM file and perform differential expression analysis
    // Save the results to ${bam_files}_DE_results
    " > ${bam_files}_DE_results/deseq2_results.txt

    echo "Differential expression completed for ${bam_files}" > ${bam_files}_DE_results/de_log.txt

    if [ $? -ne 0 ]; then
        echo "Error: Differential expression analysis failed" >&2
        exit 1
    fi
    """
}

// Optional: Aggregate QC results (e.g., using MultiQC)
process aggregateQCResults {
    input:
    path(qc_output_dir.collect())

    output:
    path("${params.outdir}/multiqc_report.html")

    script:
    """
    multiqc ${qc_output_dir.join(' ')} -o ${params.outdir}/multiqc &> ${params.outdir}/multiqc/multiqc.log

    if [ $? -ne 0 ]; then
        echo "Error: MultiQC aggregation failed" >&2
        exit 1
    fi
    """
}

