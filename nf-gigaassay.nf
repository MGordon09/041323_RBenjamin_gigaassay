#!/usr/bin/env nextflow


// WOrkflow parameter definition

params.reads =  
params.subsample = 'TRUE' // use conditional in workflow execution
params.samplesize = 0.1 // 0.1 random sample of raw file
params.adapters = // file with adapter sequences 
// params.reference = 
// params.adapters =
// params.multiqc = 
// params.outdir = 

log.info """\
    G I G A A S S A Y   P I P E L I N E
    ===================================
    reference: ${params.reference}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)


/*
 * subsample reads (for testing)
 */


process SEQTK_SAMPLE {
    tag "Subsampling $sample_id raw reads"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    val(sample_size)

    output:
    path "${sample_id}__logs"
    
    script:
    """
    seqtk sample -s100 ${reads[0]} $sample_size > ${sample_id}_${sample_size}_1.fastq.gz
    seqtk sample -s100 ${reads[1]} $sample_size > ${sample_id}_${sample_size}_2.fastq.gz
    """
    }

/*
 * qc raw reads
 */

process FASTQC_RAW {
    tag "Running FASTQC on $sample_id raw fastq files"
    publishDir "${params.outdir}/preprocessing/fastqc_raw", mode:'copy'
    
    input: 
    tupple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc_raw_logs" 

    script:
    """
    fastqc -o ${sample_id}_fastqc_logs -q ${reads}
    """
}


/*
 * trim adapters
 */

process BBDUK {
    tag "Running bbduk trimming $sample_id reads"
    publishDir "${params.outdir}/preprocessing/fastqc_raw", mode:'copy'

    input: 
    tupple val(sample_id), path(reads)
    path(adapters)

    output:
    path "${sample_id}_fastqc_raw_logs" 

    script:
    """
    mkdir -p preprocessed/adapter_removal/${sample_id}_
    fastqc -f fastq -o fastqc_raw/${sample_id}_fastqc_logs -q ${reads}
    """


}

/*
 * merge reads
 */


BBMERGE
    tag "Merging overlapping $sample_id reads"
    publishDir params.outdir, mode:'copy'


    input:
    tuple val(sample_id), path(reads)

    output:






/*
 * raw file multiqc report
 * just combine the two
 */

process MULTIQC_RAW {
    tag "Running MULTIQC on $sample_id raw fastqc files"
    publishDir params.outdir, mode:'copy'
    
    input: 
    path(*) // all input supplied in workflow
    // multiqc config???

    output:
    path "multiqc_raw.html" 

    script:
    """
    mkdir -p multiqc
    multiqc -o multiqc -n multiqc_raw.html
    """
}


workflow {
    Channel 
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    fastqc_raw_ch = FASTQC_RAW(read_pairs_ch)
    MULTIQC_RAW(fastqc_raw_ch)
}