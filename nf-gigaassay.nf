#!/usr/bin/env nextflow


/*
========================================================================================
                         nf-gigaassay
========================================================================================
 pipeline to call variants in data produced using the GigaAssay protocol 
 doi: https://www.sciencedirect.com/science/article/pii/S0888754322001847
 #### Homepage / Documentation
 https://github.com/MGordon09/041323_RBenjamin_gigaassay
----------------------------------------------------------------------------------------
*/

//TODO package each of the scripts into a script in bind and execute with find xargs 
// look at top answer below for guidance; find each file and execute seperately
//https://stackoverflow.com/questions/40700230/find-xargs-execute-chain-of-commands-for-each-file

params.caller = 'freebayes' // use freebayes as default variant caller for MNP handling
params.adapters = "$projectDir/docs/flanks.fa" // file with sequences to trim surrounding Tat contig
params.tat_flank_5prime = "GAATTC"
params.tat_flank_3prime = "GCGATCGC"
params.umi_5prime = 'TGGATCCGGTACCGAGGAGATCTG'
params.umi_3prime = 'GCGATCGC'
params.l_distance = 2 // max Levenshtein distance for clustering 
params.scripts = "$projectDir/bin" // scripts for pipeline here
params.reference = "$projectDir/docs/AF324493.2.fa" // maybe look at building the reference 
params.reference_gbk = "AF324493.2" //reference genbank accession no
params.bed ="$projectDir/docs/intervals.bed" // genomic interval (0-based) for tat 
params.min_cluster_size = 3 // minimum number of reads per cluster.. need to find sensible number. For this amount of data is 5 good? To use all just input 1 here
params.snpeff_db = "$projectDir/docs/AF324493.2"
params.snpeff_config = "$projectDir/docs/snpEff.config"
// params.multiqc = 
//params.outdir = "$projectDir/test-run'



def helpMessage() {

    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-gigaassay  -profile conda --reads 'path/to/*.fastq.gz' --skip_subsample
    Mandatory arguments:
      --reads [str]                 Full path to directory containing the input reads for analysis
      -profile [str]                Configuration profile to use. Currently supports conda, docker & singularity 
                            
    Options:
      --skip_subsample [boolean]     For testing
      --samplesize [int/float]       Number or proportion of reads to subsample
      --l_distance [int]             Levenstein distance used for starcode clustering
      --adapters [file]              Path to fasta for adapter sequences to be trimmed

    References:
      --reference [file]             Fasta file used to perform read alignment       
      --adapter [file]               Path to fasta for adapter sequences to be trimmed
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


log.info """\
    G I G A A S S A Y   P I P E L I N E
    ===================================
    reference    : ${params.reference}
    reference_genbank_accession : ${params.reference_gbk}
    reads        : ${params.reads}
    samplesize   : ${params.samplesize}
    outdir       : ${params.outdir}
    adapters     : ${params.adapters}
    umi_5'_flank : ${params.umi_5prime}
    umi_3'_flank : ${params.umi_3prime}
    clustering_distance_threshold : ${params.l_distance}
    min_cluster_size : ${params.min_cluster_size}
    """
    .stripIndent(true)


/*
 * index reference
 */

process BWA_MEM_INDEX {
    tag "Indexing reference genome"
    label 'process_high'
    publishDir "${params.outdir}/reference", mode: 'copy'

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"

    input:
    path fasta

    output:
    path("bwa"), emit: index 

    script:
    """
    mkdir bwa
    bwa index -a is -p bwa/${fasta.baseName} $fasta
    """
}


/*
 * Build SNPEff database using reference files
 */

// process SNPEFF_BUILD {
//     tag "Building SNPeff Annotation Database for $reference_gbk"
//     label 'process_low'
//     publishDir "${params.outdir}/annotation/snpeff", mode:'copy'

//     conda "bioconda::snpeff=5.1"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2' :
//         'quay.io/biocontainers/snpeff:5.1--hdfd78af_2' }"

//     input:
//     val reference_gbk
//     path scripts

//     output:
//     path("data/$reference_gbk"), emit: snpeff_db // dir containing output of SNPEff build process
//     path("snpEff.config"), emit: snpeff_config

//     script:
//     """
//     bash ${scripts}/buildDbNcbi.sh $reference_gbk
//     """
// }


/*
 * subsample reads (for testing)
 */

process SEQTK_SAMPLE {
    tag "Subsampling $sample_id raw reads"
    label 'process_low'
    publishDir "${params.outdir}/preprocessing/subsample", mode: 'copy'

    conda "bioconda::seqtk=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

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


/*
 * qc raw reads
 */

process FASTQC_RAW {
    tag "Running FASTQC on $sample_id fastq files"
    label 'process_medium'
    publishDir "${params.outdir}/preprocessing/fastqc_raw", mode: 'copy'

    conda "bioconda::fastqc=0.11.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"
    
    input: 
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc_raw_logs" 

    script:
    """
    mkdir ${sample_id}_fastqc_raw_logs
    
    fastqc \\
        --threads $task.cpus \\
        -o ${sample_id}_fastqc_raw_logs \\
        -q ${reads}
    """
}

/*
 * Losing lots of reads, may need to loosen params; consult Ronald
 */

process BBMERGE { //optimise memory usage
    tag "Merging overlapping $sample_id reads"
    label 'process_high'
    publishDir "${params.outdir}/preprocessing/merged_reads", mode: 'copy'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_merge.fastq.gz"), emit: merged_reads
    tuple val(sample_id), path("*_U{1,2}merge.fastq.gz")
    path("*.txt")
    path("*.log")

    script:
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g') 

    bbmerge.sh \\
        -Xmx\$maxmem \\
        in1=${reads[0]} in2=${reads[1]} \\
        out=${sample_id}_merge.fastq.gz \\
        outu1=${sample_id}_U1merge.fastq.gz outu2=${sample_id}_U2merge.fastq.gz \\
        outinsert=${sample_id}_insertinfo.txt ihist=insert-histogram.txt \\
        minoverlap=20 maxratio=0.15 \\
        &> ${sample_id}_bbbmerge.log    
    """
}

/*
 * trim adapters only to leave biological sequence
 * maybe try minlen filtering but can only do after quality trimming
 */

process CUTADAPT_TRIM {
    tag "Trimming reads from $sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/preprocessing/trimmed_reads", mode:'copy'

    conda "bioconda::cutadapt=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
        'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"

    input: 
    tuple val(sample_id), path(reads)
    val tat_flank_5prime
    val tat_flank_3prime

    output:
    tuple val(sample_id), path("*trim.fastq.gz"), emit: trimmed_reads
    path("*.log")

    script:
    """
    cutadapt \\
        --cores $task.cpus \\
        -g ${tat_flank_5prime}...${tat_flank_3prime} \\
        ${reads} -o ${sample_id}_trim.fastq.gz \\
        -q 10 \\
        --minimum-length 320 \\
        --discard-untrimmed \\
        --report minimal \\
        > ${sample_id}_cutadapt.log
    """
}


/*
 * extract UMIs using flanking sequences and place in read header
 * -g both adapters required, filtering by read length to ensure they are equal
 */

process CUTADAPT_UMI {
    tag "Extracting UMIs from $sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/umi/umi_barcodes", mode:'copy'

    conda "bioconda::cutadapt=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
        'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"

    input: 
    tuple val(sample_id), path(reads)
    val umi_5prime

    output:
    tuple val(sample_id), path("*_umi.fastq.gz"), emit: umi_reads
    path("*.log")

    script:
    """
    cutadapt \\
        --cores $task.cpus \\
        -g ${umi_5prime} \\
        ${reads} \\
        -o ${sample_id}_umi.fastq.gz \\
        --report minimal \\
        > ${sample_id}_cutadapt.log
    """
}

/*
 * Create barcode sequence clusters
 * Default algorithm with Levenshtein distance=2
 * Use seq-id in output to demultiplex input files based on barcode
 */

process STARCODE_CLUSTERING { 
    tag "Creating global clusters from UMI barcodes"
    label 'process_medium'
    publishDir "${params.outdir}/clustering/starcode", mode:'copy'

    conda "bioconda::starcode==1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/starcode:1.4--hec16e2b_3' :
        'quay.io/biocontainers/starcode:1.4--hec16e2b_3' }"

    input: 
    tuple val(sample_id), path(reads)
    val l_distance

    output:
    tuple val(sample_id), path("*clusters.txt"), emit: clusters_file
    path("*clusters.log")

    script:
    """
    starcode \\
        --threads $task.cpus \\
        -i <(cat ${reads} | gzip -cd) \\
        --output ${sample_id}_starcode_umi_clusters.txt \\
        --sphere \\
        --seq-id \\
        --dist ${l_distance} \\
        2> ${sample_id}_starcode_umi_clusters.log

    #gzip ${sample_id}_starcode_umi_clusters.txt
    """
 }

/*
 * Filter out noisy clusters with less than N reads (less than 5 reads for now)
 */

process FILTER_CLUSTERS { 
    tag "Removing clusters with few reads"
    label 'process_single'
    publishDir "${params.outdir}/clustering/starcode", mode:'copy'

    conda "bioconda gzip==1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11' :
        'quay.io/biocontainers/gzip:1.11' }"

    input: 
    tuple val(sample_id), path(clusters_file)
    val cluster_size

    output:
    tuple val(sample_id), path("*clusters.clean.txt"), emit: clusters_file
    path("*clusters.clean.log")


    script:
    """
    #zless ${clusters_file} | awk 'NR==1 { print } NR != 1 && \$2 >= ${cluster_size} { print }' | gzip > ${sample_id}_starcode_umi_clusters.clean.txt.gz
    cat ${clusters_file} | awk 'NR==1 { print } NR != 1 && \$2 >= ${cluster_size} { print }' > ${sample_id}_starcode_umi_clusters.clean.txt

    start_seq=`wc -l  ${clusters_file} | awk '{print \$1}'`
    end_seq=`wc -l ${sample_id}_starcode_umi_clusters.clean.txt | awk '{print \$1}'`

    echo "Initial clusters: \$start_seq Remaining clusters: \$end_seq" >  ${sample_id}_starcode_umi_clusters.clean.log
    """
 }


/*
 * Demultiplex samples fastq and assign reads to barcode groups generated using starcode clustering
 * output should be 1 file per sample per cluster: ${sample}_${cluster}.fastq
 */

process DEMULTIPLEX_BARCODES {
    tag "Demultiplex $sample_id fastq files"
    label 'process_medium'
    publishDir "${params.outdir}/clustering/reads/$sample_id", mode:'copy'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(sample_id), path(reads), path(clusters_file)
    path scripts

    output:
    //tuple val(sample_id), path("*.demux.fastq.gz"), emit: demux_reads
    tuple val(sample_id), path("demux.files.tar.gz"), emit: demux_folder
    path("*demux.log")

    script:
    """
    # create output dir
    mkdir -p ./demux.files

    python ${scripts}/demux_index_optim_v2.py \\
        --fastq_path ${reads} \\
        --barcode_path ${clusters_file} \\
        --sample_name ${sample_id} \\
        --output_dir './demux.files/' \\
        > ${sample_id}_demux.log


    # compress files
    tar -vcf demux.files.tar.gz ./demux.files/ 
    """
}

/*
 * Align each sequence to its reference using  BWA-MEM aligner
 * Assigning read groups with -R indicating $barcode & $sample origin origin.
 * Can use this info to manipulate sequences downstream
 */

process BWA_MEM_ALIGN {
    tag"Aligning $sample_id fastqs against $index"
    label 'process_medium'
    publishDir "${params.outdir}/alignment/bwa-mem/$sample_id", mode:'copy'

    conda "bioconda::mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40==c56a3aabc8d64e52d5b9da1e8ecec2031668596d-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:c56a3aabc8d64e52d5b9da1e8ecec2031668596d-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:c56a3aabc8d64e52d5b9da1e8ecec2031668596d-0' }"

    input:
    path index
    tuple val(sample_id), path(demux_folder)
    path scripts

    output:
    tuple val(sample_id), path("bam.files.tar.gz"), emit: bam_folder

    shell:
    """
    tar -xf ${demux_folder}

    # iterate through each file in the tar directory 
    find ./demux.files -type f -name "*.demux.fastq.gz" -exec ${scripts}/run_bwa_mem_align.sh ${task.cpus} \\{\\} \\;

    # cleanup
    mkdir ./bam.files
    find . -type f \\( -name "*.bwa.err" -o -name "*.sorted.bam" \\) -exec mv \\{\\} ./bam.files \\;
    #mv *{bwa.err,sorted.bam} ./bam.files
    
    # compress files
    tar -vcf bam.files.tar.gz ./bam.files/
    """
}

/*
 * Variant Calling
 * Using  a combination of freebayes & bcftoolsq
 *   freebayes --ploidy 2 --targets $bed --min-alternate-fraction 0.5 --min-alternate-count 1 --min-mapping-quality 1  --min-base-quality 3 --use-best-n-alleles=1 -f $index $bam > ${sample_id}_${barcode}.freebayes.vcf
 * --use-best-n-alleles as only one site is mutated only take best site (SNPS vs MNPs? how does freebayes distinguish take 3 for now)
 *    find ./demux.files -type f -name "*.demux.fastq.gz" -exec bash -c "ls  \"\$0\"" " {} \; 

 */

process FREEBAYES {

    tag "Calling variants on $sample_id demultiplexed files"
    label 'process_medium'
    publishDir "${params.outdir}/calling/${params.caller}/$sample_id", mode:'copy'

    conda "bioconda::freebayes=1.3.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2' :
        'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2' }"

    input:
    path index
    tuple val(sample_id), path(bam_folder)
    path scripts

    output:
    tuple val(sample_id), path("variant.files.tar.gz"), emit: variants_folder

    script:
    """
    tar -xf ${bam_folder}

    # iterate through each file in the tar directory 
    find ./bam.files -type f -name "*.sorted.bam" -exec ${scripts}/run_freebayes.sh ${index} \\{\\} \\;

    # cleanup
    mkdir ./variant.files
    find . -type f -name "*.vcf*" -exec mv \\{\\} ./variant.files \\;

    #mv *{.vcf,.vcf.err} ./variant.files
    
    # compress files
    tar -vcf variant.files.tar.gz ./variant.files/
    """
}

//bcftools variant calling options
// mpileup:
//-I skip indelcalling, -Ou output uncompressed bcf
// -R limit pileups to target region in bed file (genomic coordinates of tat)
// -Ou stream output in uncompressed bcf format to speed up process (avoid conversion between steps)
// call:
// skip indels and output variants only
// use old consensus caller algorithm as only expecting one allele per site in this experiment. -m handles multi-allelic sites better
//


process BCFTOOLS_MPILEUP_CALL {

    tag "Calling variants on $sample_id demultiplexed files"
    label 'process_medium'
    publishDir "${params.outdir}/calling/${params.caller}/$sample_id", mode:'copy'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    path index
    tuple val(sample_id), path(bam_folder)
    path scripts 

    output:
    tuple val(sample_id), path("variant.files.tar.gz"), emit: variants_folder

    script:
    """
    tar -xf ${bam_folder}

    # iterate through each file in the tar directory 
    find ./bam.files -type f -name "*.sorted.bam" -exec ${scripts}/run_bcftools.sh ${index} \\{\\} \\;

    # cleanup
    mkdir ./variant.files
    find . -type f -name "*.bcf*" -exec mv \\{\\} ./variant.files \\;
    #mv *{.bcf,.bcf.err} ./variant.files
    
    # compress files
    tar -vcf variant.files.tar.gz ./variant.files/
    """
}

//bcftools variant calling options
// view:
// filter variants outside the target sequence (tat) & those with low supporting information: require at least 60% of reads to support call
// dropping this filtering step above; not needed as require 60% of reads to map anyway; bcftools view --targets "AF324493.2:5829-6044,AF324493.2:8368-8414" -q 0.6:nref $vcf -Ou |
// annotate vcf id column with sample_barcode info 
// normalise the variant calls (lelft-align & max parisomony) 
// output is compressed bcf for feeding into concat input
// improve target option here by modifying above

process BCFTOOLS_VIEW_ANNOTATE_NORM_INDEX {
    tag "Removing low-support & normalising variants"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/bcftools/norm/${params.caller}/$sample_id", mode:'copy'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"
    
    input:
    path fasta
    tuple val(sample_id), path(variants_folder)
    path scripts 

    output:
    tuple val(sample_id), path("norm.files.tar.gz"), emit: norm_variants_folder

    script:
    """
    tar -xf ${variants_folder}

    # iterate through each vcf/bcf file 
    find variant.files/ -type f \\( -name "*.bcf" -o -name "*.vcf" \\) -exec ${scripts}/run_bcftools_norm.sh ${fasta} \\{\\} \\;

    # cleanup
    mkdir ./norm.files
    find . -type f -name "*.norm.bcf*" -exec mv \\{\\} ./norm.files \\;
    #mv *{.norm.bcf,.norm.bcf.csi} ./norm.files
    
    # compress files
    tar -vcf norm.files.tar.gz ./norm.files/
    """

}

// bcftools concat to concatenate vcfs per sample
// vcf files

process BCFTOOLS_CONCAT {
    tag " Concatenating sample $sample_id vcfs"
    label 'process_low'
    publishDir "${params.outdir}/annotation/bcftools/concat/${params.caller}/$sample_id", mode:'copy'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(sample_id), path(norm_variants_folder)

    output:
    tuple val(sample_id), path("*.combined.vcf"), emit: vcf

    script:
    """
    tar -xf ${norm_variants_folder}

    # create list of input files
    find -L ./ -name "*.norm.bcf" > vcfs_list
    
    bcftools concat \\
        --allow-overlaps \\
        --file-list ./vcfs_list \\
        --output ${sample_id}.norm.combined.vcf
    """
}

/*
 * Annotate variants in VCF using SNPEffect 
 * HGVS.c and HGVS.p annotation for codon and aa substitutions
 */

process SNPEFF_ANNO {
    tag "Annotating variants in $sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/snpeff/${params.caller}/$sample_id", mode:'copy'

    conda "bioconda::snpeff=5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2' :
        'quay.io/biocontainers/snpeff:5.1--hdfd78af_2' }"

    input:
    path snpeff_db
    path snpeff_config //path snpeff_config make optional as already available through db path
    val reference_v
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("*.ann.vcf"), emit: ann_vcf
    tuple val(sample_id), path("*.snpEff_summary.html"), emit: ann_html

    script:
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """  
    snpEff ann  \\
        -Xmx${avail_mem}M \\
        -nodownload -v ${reference_v} \\
        -c ${snpeff_config} \\
        -csvStats -hgvs \\
        -noStats \\
     ${vcf} > ${sample_id}.ann.vcf

    mv ./snpEff_summary.html ./${sample_id}.snpEff_summary.html
    """
}

/*
 * raw file multiqc report
 * just combine the two
 */

process MULTIQC_RAW {
    tag "Running MULTIQC on raw fastqc files"
    publishDir "${params.outdir}/preprocessing/multiqc_raw", mode: 'copy'
    
    conda "bioconda::multiqc=1.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0' }"

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

    Channel //reference fasta
        .value(params.reference)
        .set { reference_ch }      

    Channel
        .value(params.reference_gbk)
        .set { reference_gbk_ch }

    Channel
        .value(params.bed)
        .set { bed_ch }

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

    Channel
        .value(params.min_cluster_size)
        .set { min_cluster_size_ch }

    Channel
        .value(params.snpeff_config)
        .set { snpeff_config_ch }
    
    Channel
        .value(params.snpeff_db)
        .set { snpeff_db_ch }



    // adding snpeff db and config
    Channel
        .value(params.scripts)
        .set { scripts_ch }
    
    
    //workflow

    if (!params.skip_subsample) {
        read_input_ch    = SEQTK_SAMPLE(read_pairs_ch, sample_size_ch)
    } else {
        read_input_ch    = read_pairs_ch
    }

    // Build reference index & snpeff db
    index_ch         = BWA_MEM_INDEX(reference_ch)
    //snpeffdb_ch      = SNPEFF_BUILD(reference_gbk_ch, scripts_ch)

    multiqc_input_ch = FASTQC_RAW(read_input_ch)
    merge_reads_ch   = BBMERGE(read_input_ch)
    trimmed_reads_ch = CUTADAPT_TRIM(merge_reads_ch.merged_reads, flank_5_ch, flank_3_ch)
    umi_reads_ch     = CUTADAPT_UMI(trimmed_reads_ch.trimmed_reads, umi_5_ch)
    
    // cluster umi for each sample
    clustered_umi_ch = STARCODE_CLUSTERING(umi_reads_ch.umi_reads, lev_dist_ch)

    // filter the clusters to only keep clusters with >= n
    clean_clustered_umi_ch = FILTER_CLUSTERS(clustered_umi_ch.clusters_file, min_cluster_size_ch)


    //combine input channels by sample keys
    demulti_ch = DEMULTIPLEX_BARCODES(trimmed_reads_ch.trimmed_reads.combine(clean_clustered_umi_ch.clusters_file, by:0), scripts_ch)

    // take the output dir 

    //grouped_by_sample_barcode = demulti_ch.demux_reads.flatMap { sample, files -> files.collect { [sample, it] } } //flatMap transform ch into flattened [sample, file_path] tuples chs.
    //    .map { sample, file_path -> //create new tuple using map
    //    def matcher = (file_path =~ /.*\/(.*)_(.*)\.demux\.fastq\.gz/) //capture both barcode & sample groups
    //  def (barcode, sample_id) = matcher ? [matcher[0][1], matcher[0][2]] : null
    //return tuple(sample_id, barcode, file_path)
    //}
    
    //bam_ch = BWA_MEM_ALIGN(index_ch,grouped_by_sample_barcode)
    bam_ch = BWA_MEM_ALIGN(index_ch,demulti_ch.demux_folder,scripts_ch)

    // calling
    if (params.caller == 'freebayes') {
        vcf_ch = FREEBAYES(reference_ch, bam_ch.bam_folder, scripts_ch)
    } else if (params.caller == 'bcftools') {
        vcf_ch = BCFTOOLS_MPILEUP_CALL(reference_ch, bam_ch.bam_folder, scripts_ch)
    }

    // filter off-target variants those with lower read support. annotate ID col with sample_barcode and normalise variants (left align, max parisomony representation) 
    norm_vcf_ch = BCFTOOLS_VIEW_ANNOTATE_NORM_INDEX(reference_ch, vcf_ch.variants_folder, scripts_ch)


    // concatenate variants in vcf files
    vcf_concat_ch = BCFTOOLS_CONCAT(norm_vcf_ch.norm_variants_folder)

    // annotate comb variants file w snpEff
    ann_ch = SNPEFF_ANNO(snpeff_db_ch, snpeff_config_ch, reference_gbk_ch, vcf_concat_ch.vcf)

    MULTIQC_RAW(multiqc_input_ch.collect())

}
