#!/usr/bin/env nextflow

// WOrkflow parameter definition

params.reads =  "$projectDir/data/fastq/SRR*_{1,2}.fastq.gz"
params.subsample = 'FALSE' // use conditional in workflow execution
params.samplesize = 0.1 // 0.1 random sample of raw file
params.adapters = "$projectDir/flanks.fa" // file with sequences to trim surrounding Tat contig
params.tat_flank_5prime = "GAATTC"
params.tat_flank_3prime = "GCGATCGC"
params.umi_5prime = 'TGGATCCGGTACCGAGGAGATCTG'
params.umi_3prime = 'GCGATCGC'
params.l_distance = 2 // max maximum Levenshtein distance for starcode clustering 
// params.reference = 
// params.multiqc = 
// params.outdir = 

log.info """\
    G I G A A S S A Y   P I P E L I N E
    ===================================
    reference: ${params.reference}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    adapters     : ${params.adapters}
    umi_5'_flank : ${params.umi_5prime}
    umi_3'_flank : ${params.umi_3prime}
    clustering_distance_threshold : ${params.l_distance}
    """
    .stripIndent(true)


/*
 * subsample reads (for testing)
 */

process SEQTK_SAMPLE {
    tag "Subsampling $sample_id raw reads"
    publishDir "${params.outdir}/preprocessing/subsample", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val(samplesize)

    output:
    tuple val(sample_id), path("*sub*.fastq.gz")
    
    """
    seqtk sample -s100 ${reads[0]} $samplesize | gzip --no-name > ${sample_id}_sub.R1.fastq.gz
    seqtk sample -s100 ${reads[1]} $samplesize | gzip --no-name > ${sample_id}_sub.R2.fastq.gz
    """
    }

    //try later as more robust 
    /* printf "%s\\n" $reads | while read f;
    do 
        seqtk sample -s100 \$f $samplesize | gzip --no-name > \$f
    done */



/*
 * qc raw reads
 */

process FASTQC_RAW {
    tag "Running FASTQC on $sample_id fastq files"
    publishDir "${params.outdir}/preprocessing/fastqc_raw", mode: 'copy'
    
    input: 
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc_raw_logs" 

    script:
    """
    mkdir ${sample_id}_fastqc_raw_logs
    fastqc -o ${sample_id}_fastqc_raw_logs -q ${reads}
    """
}

/*
 * merge reads - not impacted by artifical seqs as on 3' ends
 * Losing lots of reads, may need to loosen params; consult Ronald
 */

process BBMERGE { //optimise memory usage
    tag "Merging overlapping $sample_id reads"
    publishDir "${params.outdir}/preprocessing/merged_reads", mode: 'copy'


   input:
   tuple val(sample_id), path(reads)

   output:
   tuple val(sample_id), path("*_merge.fastq.gz"), emit: merged_reads
   tuple val(sample_id), path("*_U{1,2}merge.fastq.gz")
   path("*.txt")
   path("*.log")

   script:
   """
   bbmerge.sh \
   loose \
   in1=${reads[0]} in2=${reads[1]} \
   out=${sample_id}_merge.fastq.gz \
   outu1=${sample_id}_U1merge.fastq.gz outu2=${sample_id}_U2merge.fastq.gz \
   outinsert=${sample_id}_insertinfo.txt ihist=insert-histogram.txt \
   &> ${sample_id}_bbbmerge.log    
   """
}

//todo check read proportion
// trim linked adapters using bbduk (trim to digestion sites), could also try use bbduk to extract UMI
// maybe try filter reads based on insert size??


/*    bbmerge-auto.sh -Xmx26g \
   in1=${reads[0]} in2=${reads[1]} \
   adapter1=$adapters1 adapter2=$adapters2 \
   out=${sample_id}_merge.fastq.gz \
   outu1=${sample_id}_U1_merge.fastq.gz outu2=${sample_id}_U2_merge.fastq.gz \
   outinsert=${sample_id}_insertinfo.txt outadapter=${sample_id}_adapterinfo.txt \
   verystrict=t minoverlap=30 
 */

/*
 * trim adapters only to leave biological sequence
 * maybe try minlen filtering but can only do after quality trimming
 */

/* process BBDUK_ADAPTER {
    tag "Running bbduk trimming $sample_id reads"
    publishDir "${params.outdir}/preprocessing/trimmed_reads", mode:'copy'

    input: 
    tuple val(sample_id), path(reads)
    path(adapters)

    output:
    tuple val(sample_id), path("*trim.fastq.gz"), emit: trimmed_reads
    path("*.txt")
    path("*.log")

    script:
    """
    bbduk.sh \
    in=${reads} out=${sample_id}_trim.fastq.gz \
    ref=${adapters} \
    ktrim=rl  k=5 hdist=0 tbo \
    maq=10 \
    bhist=${sample_id}_bhist.txt qhist=${sample_id}_qhist.txt gchist=${sample_id}_gchist.txt aqhist=${sample_id}_aqhist.txt lhist=${sample_id}_lhist.txt gcbins=auto \
    &> ${sample_id}_bbbduk.log
    """
} */

process CUTADAPT_TRIM {
    tag "Trimming reads from $sample_id"
    publishDir "${params.outdir}/preprocessing/trimmed_reads", mode:'copy'

    input: 
    tuple val(sample_id), path(reads)
    val tat_flank_5prime
    val tat_flank_3prime

    output:
    tuple val(sample_id), path("*trim.fastq.gz"), emit: trimmed_reads
    path("*.log")

    script:
    """
    cutadapt \
    -g ${tat_flank_5prime}...${tat_flank_3prime} \
    ${reads} -o ${sample_id}_trim.fastq.gz \
    --discard-untrimmed \
    --report minimal \
    > ${sample_id}_cutadapt.log
    """
}


// length filtering seems to be quality filtering the reads?

/* process REFORMAT_FILTER {
    tag "Length filtering $sample_id reads"
    publishDir "${params.outdir}/preprocessing/filtered_reads", mode:'copy'

    input: 
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*filtered.fastq.gz")
    path("*.txt")
    path("*.log")

    script:
    """
    reformat.sh \
    in=${reads} out=${sample_id}_filtered.fastq.gz \
    minlength=461 \
    bhist=${sample_id}_bhist.txt qhist=${sample_id}_qhist.txt qchist=${sample_id}_qchist.txt aqhist=${sample_id}_aqhist.txt bqhist=${sample_id}_bqhist.txt lhist=${sample_id}_lhist.txt \
    &> ${sample_id}_reformat.log
    """
} */

/*
 * extract UMIs using flanking sequences and place in read header
 * -g both adapters required, filtering by read length to ensure they are equal
 */

/* process CUTADAPT_UMI {
    tag "Extracting UMIs from $sample_id"
    publishDir "${params.outdir}/preprocessing/umi_fastq", mode:'copy'

    input: 
    tuple val(sample_id), path(reads)
    val umi_5
    val umi_3

    output:
    tuple val(sample_id), path("*_umi.fastq.gz"), emit: umi_reads
    path("*.log")

    script:
    """
    cutadapt \
    -g "${umi_5}...${umi_3}" \
    ${reads} -o ${sample_id}_umi.fastq.gz \
    --minimum-length 32 --maximum-length 32 \
    --discard-untrimmed \
    --report minimal \
    > ${sample_id}_cutadapt.log
    """
} */

process CUTADAPT_UMI {
    tag "Extracting UMIs from $sample_id"
    publishDir "${params.outdir}/umi/umi_barcodes", mode:'copy'

    input: 
    tuple val(sample_id), path(reads)
    val umi_5prime

    output:
    tuple val(sample_id), path("*_umi.fastq.gz"), emit: umi_reads
    path("*.log")

    script:
    """
    cutadapt \
    -g ${umi_5prime} \
    ${reads} -o ${sample_id}_umi.fastq.gz \
    --discard-untrimmed \
    --report minimal \
    > ${sample_id}_cutadapt.log
    """
}


/*
 * extract UMI using flanking sequence and incorporate into the read header
 * regex pattern for umi capture: match flanking regions (allow 2 mismatches ~0.1e) and extract 32bp UMI
 */


process UMITOOLS_EXTRACT {
    tag "Extracting UMIs from $sample_id"
    publishDir "${params.outdir}/umi/umi_reads", mode:'copy'

    input: 
    tuple val(sample_id), path(reads)
    val umi_5prime

    output:
    tuple val(sample_id), path("*_umi.header.fastq.gz"), emit: umi_reads
    path("*.log")

    script:
    """
    umi_tools extract \
    --stdin=${reads} \
    --stdout=${sample_id}_umi.header.fastq.gz \
    --extract-method=regex \
    --bc-pattern='.+(?P<discard_1>${umi_5prime}){s<=2}(?P<umi_1>.{32})' \
    --log=${sample_id}_processed.log
    """
}
/*
backup umi-tools extract with full lenght reads
    script:
    """
    umi_tools extract \
    --stdin=${reads} \
    --stdout=${sample_id}_umi.header.fastq.gz \
    --extract-method=regex \
    --bc-pattern='.+(?P<discard_1>${umi_5}){s<=2}(?P<umi_1>.{32})(?P<discard_2>${umi_3}.+){s<=2}' \
    --log=${sample_id}_processed.log
    """
} */

/*
 * merge files with barcodes, unzip & cluster
 * Levenshtein distance=4 to form clusters
 * Using spheres algorithm  intra-cluster dist<=4 (message passing builds by merging clusters)
 * Use seq-id in output to demultiplex input files based on barcode

 */


process STARCODE_CLUSTERING { 
    tag "Creating global clusters from UMI barcodes"
    publishDir "${params.outdir}/clustering/starcode", mode:'copy'

    input: 
    path(reads)
    val l_distance

    output:
    path("*clusters.txt")
    path("*clusters.log")

    script:
    """
    starcode \
    -i <(cat ${reads} | gzip -cd) \
    --output starcode_umi_clusters.txt \
    --sphere \
    --print-clusters \
    --dist $l_distance \
    2> starcode_umi_clusters.log
    """
 }

/* process DEMULTIPLEX_BARCODES {
    tag "Creating global clusters from UMI barcodes"
    publishDir "${params.outdir}/clustering/starcode_reads", mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    path clusters_file

    output:
    


} */
/* 
    starcode \
    --input merged_umi.fastq.gz \
    --output starcode_umi_clusters.txt \
    --sphere \
    --dist $l_distance \
    --non-redundant \
    --seq-id \
    2> starcode_umi_clusters.log
 */


// The sequence reads are demultiplexed into subsets of read sequences for each cell clone based on UMI-barcode groups
// so output {cell_clone}_{barcode}_fastq.gz

/*
 * Starcode read clustering
 */


/*
 * raw file multiqc report
 * just combine the two
 */

process MULTIQC_RAW {
    tag "Running MULTIQC on raw fastqc files"
    publishDir "${params.outdir}/preprocessing/multiqc_raw", mode: 'copy'
    
    input: 
    path("*") // all input supplied in workflow multiqc config???

    output:
    path "multiqc_raw.html" 

    script:
    """
    multiqc . -n multiqc_raw.html
    """
}


workflow {
    Channel // read input channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    Channel //subsample val input channel
        .value(params.samplesize)
        .set { sample_size_ch }

    Channel 
        .value(params.samplesize)
        .set { sample_size_ch }

    Channel 
        .value(params.adapters)
        .set { adapter_ch }

    Channel 
        .value(params.tat_flank_5prime)
        .set { flank_5_ch }
    
    Channel 
        .value(params.tat_flank_3prime)
        .set { flank_3_ch }
    
    Channel 
        .value(params.umi_5prime)
        .set { umi_5_ch }
    
    Channel 
        .value(params.umi_3prime)
        .set { umi_3_ch } 

    Channel 
        .value(params.l_distance)
        .set { lev_dist_ch }

    //workflow

    read_input_ch    = SEQTK_SAMPLE(read_pairs_ch, sample_size_ch)
    multiqc_input_ch = FASTQC_RAW(read_input_ch)
    merge_reads_ch   = BBMERGE(read_input_ch)
    trimmed_reads_ch = CUTADAPT_TRIM(merge_reads_ch.merged_reads, flank_5_ch, flank_3_ch)
    //trimmed_reads_ch = BBDUK_ADAPTER(merge_reads_ch.merged_reads, adapter_ch)
    //filtered_read_ch = REFORMAT_FILTER(trimmed_reads_ch.trimmed_reads)
    umi_reads_ch     = CUTADAPT_UMI(trimmed_reads_ch.trimmed_reads, umi_5_ch)
    umi_h_reads_ch   = UMITOOLS_EXTRACT(trimmed_reads_ch.trimmed_reads, umi_5_ch)
    //umi_h_reads_ch   = UMITOOLS_EXTRACT(trimmed_reads_ch.trimmed_reads, umi_5_ch, umi_3_ch)
    //clustered_umi_ch = STARCODE_CLUSTERING(umi_reads_ch.umi_reads.flatten().collect(), lev_dist_ch)
    umi_reads_ch.umi_reads.flatten().filter( ~/^\/.*\.fastq\.gz$/ ).collect().view()
    clustered_umi_ch = STARCODE_CLUSTERING(umi_reads_ch.umi_reads.flatten().filter( ~/^\/.*\.fastq\.gz$/ ).collect(), lev_dist_ch) //flattens tuple & extracts elements ending in *.fastq.gz
    
    
    //umi_h_reads_ch   = UMITOOLS_EXTRACT(merge_reads_ch.merged_reads, umi_5_ch, umi_3_ch)
    //clusters_ch      = STARCODE_CLUSTERING(umi_reads_ch.)
    MULTIQC_RAW(multiqc_input_ch.collect())


}