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

// notes: reason for such long barcode is to ensure different cells/tat constructs don't have the same barcode by chance; allows us to count
// need to identify 'true' barcodes from sequencing artifacts
// need to find how many things truly sequenced!
// plot cumulative fraction of reads and find the knee
// dont think indexing file is necessary

//  write a python script to merge adjacent SNPs into a MNP or look at using
//Read the input VCF file using a VCF parsing library or custom VCF parsing code in Python.
//Iterate over the variants in the VCF file.
//Identify adjacent SNPs based on their positions and alleles. You may consider a specific distance threshold or genomic context for merging.
//Merge the adjacent SNPs into MNPs by combining their alleles and adjusting the position information.
//Update the INFO, FORMAT, and genotype fields as necessary to reflect the merged MNPs.
//Write the merged MNPs into a new VCF file or overwrite the input file, depending on your requirements

//# ######## TODO #############
// merge vcfs using bcftools merge; sort first?
// bcftools csq wont work as need an ensmbl specific gff file format; lok 
// counting reads using RPM; how does this control for sequencing depth??
// extract the aa substitution, the supporting read depth (number of counts in each barcode),
// identify the barcode groups for common variants and multiplex this data 
// look at freebayes for calling MNPs
// use sabre https://github.com/najoshi/sabre https://astrobiomike.github.io/amplicon/demultiplexing for demultiplexing as likely a lot faster as written in C; adv is it also removes the barcode from the file; easier for mapping etc and i get output of n reads in each barcode
// for starcode, i need to run cutadapt twice (one for umi, one for reads).. if I use  umi-tools whitelist, I can find the UMIs and error-matched UMI, & reformat the output for sabre input    
// just need a tool to identify the barcodes; could use starcode or umi-tools whitelist; for whitelist, extract first col, 
// maybe look at cutadapt for demultiplexing (maximum default error rate is 10%, so 2-3 bp error allowed) https://cutadapt.readthedocs.io/en/stable/guide.html#named-adapters
// use umi-tools white-list to extract the likely cell barcodes & write to file; just take first col
// cutadapt alt: extract barcode and write to file (maybe umi-tools)
// take starcode clustering output with sequence & sorted cluster, split into fasta format each fasta record has multiple rows) and then
// add snpcconfig to bin (maybe add both to see what works) and specify with -c on the command line
// big q: do we perform clustering after merging or not? Clustering with plasmid file
// Filter the variants by the intersection between plasmid and other file; use bcftools for this
// this is the 'final' output; all barcodes should be present in both groups, so any that are not can be dropped from further analysis
// Add resource requirements per process;
// create singularity & conda env to run pipeline  on HPC
// Inspect length of umi fastq
// output needs to be: sample cell-line gft-exp barcode variant call support(n reads Q. what does this mean when there are multiple? codon substitutions? are they all 100%?) codon_substitutions aa_substitutions n_codon substitions n_aa_substitutions
//filtering: remove low-freq calls as should be present in all samples? either i) don't call below a certain fraction or ii) look at intersection of each barcode for all samples and only extract the correct one
// To the vcf file, add HGVS.c and HGVS.p information


/// Notes
//names: 1_Tatlib_JLTRG_GFP_high 1_ replicate
// Bloom lab and FOwler group: review 
//  plasmid shoudl have all the barcodes; mother 

//Custom Python scripts are used to identify the amino acid substitution for the VCFs, the number of reads for each UMI-barcode in each sample (this will just be read depth as a proportion of total fastq file depth), and the barcodes groups for cells with the same amino acid substitutions. The PyVCF library was used in scripts that gathered the information for each variant from the VCF files. å

// ####################### Default Parameter Settings #######################

params.help = false
params.reads =  "$projectDir/data/fastq/SRR*_{1,2}.fastq.gz"
params.skip_subsample = '' // option to skip
params.samplesize = 0.1 // 0.1 random sample of raw file
params.caller = 'freebayes' // use freebayes as default variant caller as better at handling MNPs
params.adapters = "$projectDir/docs/flanks.fa" // file with sequences to trim surrounding Tat contig
params.tat_flank_5prime = "GAATTC"
params.tat_flank_3prime = "GCGATCGC"
params.umi_5prime = 'TGGATCCGGTACCGAGGAGATCTG'
params.umi_3prime = 'GCGATCGC'
params.l_distance = 2 // max Levenshtein distance for clustering 
params.scripts = "$projectDir/bin" // scripts for pipeline here
params.reference = "$projectDir/docs/AF324493.2.fa" // maybe look at building the reference 
params.reference_gbk = "AF324493.2" //reference genbank accession no
params.bed ="$projectDir/bin/intervals.bed" // genomic interval (0-based) for tat 

// params.multiqc = 
// params.outdir = 



def helpMessage() {

    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-gigaassay  -profile conda --reads 'path/to/*.fastq.gz' --skip_subsample
    Mandatory arguments:
      --reads [str]                 Full path to directory containing the input reads for analysis
      -profile [str]                Configuration profile to use. Currently supports conda & singularity 
                            
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
    """
    .stripIndent(true)


/*
 * index reference
 */

process BWA_MEM_INDEX {
    tag "Indexing reference genome"
    publishDir "${params.outdir}/reference", mode: 'copy'

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

process SNPEFF_BUILD {
    tag "Building SNPeff Annotation Database for $reference_gbk"
    publishDir "${params.outdir}/annotation/snpeff", mode:'copy'

    input:
    val reference_gbk
    path scripts

    output:
    path("data/$reference_gbk"), emit: snpeff_db // dir containing output of SNPEff build process
    path("snpEff.config"), emit: snpeff_config

    script:
    """
    bash ${scripts}/buildDbNcbi.sh $reference_gbk
    """
}


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
   in1=${reads[0]} in2=${reads[1]} \
   out=${sample_id}_merge.fastq.gz \
   outu1=${sample_id}_U1merge.fastq.gz outu2=${sample_id}_U2merge.fastq.gz \
   outinsert=${sample_id}_insertinfo.txt ihist=insert-histogram.txt \
   minoverlap=20 maxratio=0.15 \
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
    -q 10 \
    --minimum-length 300 \
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
    minl  ength=461 \
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
    --report minimal \
    > ${sample_id}_cutadapt.log
    """
}

/*
 * extract UMI using flanking sequence and incorporate into the read header
 * regex pattern for umi capture: match flanking regions (allow 2 mismatches ~0.1e) and extract 32bp UMI
 */

// process UMITOOLS_EXTRACT {
//     tag "Extracting UMIs from $sample_id"
//     publishDir "${params.outdir}/umi/umi_reads", mode:'copy'

//     input: 
//     tuple val(sample_id), path(reads)
//     val umi_5prime

//     output:
//     tuple val(sample_id), path("*_umi.header.fastq.gz"), emit: umi_reads
//     path("*.log")

//     script:
//     """
//     umi_tools extract \
//     --stdin=${reads} \
//     --stdout=${sample_id}_umi.header.fastq.gz \
//     --extract-method=regex \
//     --bc-pattern='.+(?P<discard_1>${umi_5prime}){s<=2}(?P<umi_1>.{32})' \
//     --log=${sample_id}_processed.log
//     """
// }
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
 * Levenshtein distance=2 to form clusters
 * Using spheres algorithm  intra-cluster dist<=4 (message passing builds by merging clusters)
 * Use seq-id in output to demultiplex input files based on barcode

 */

/// TODO this is merged clustering approach - try if ours doesnt work

// process CLUSTER_DEMULTIPLEX { 
//     tag "Creating global UMI clusters & demultiplexing reads"
//     publishDir "${params.outdir}/clustering_demultiplexing", mode:'symlink'

//     input: 
//     path umi
//     path reads
//     val l_distance
//     path scripts 

//     output:
//     path("*clusters.txt.gz"), emit: clusters_file
//     path("*clusters.log")
//     path("*.demux.fastq.gz"), emit: demux_reads


//     script:
//     """
//     #merge reads & fastq
//     cat ${umi} > comb_umi
//     wc -l comb_umi
//     cat ${reads} > comb_reads
//     wc -l comb_reads

//     #create clusters


//     starcode \
//     -i <(cat comb_umi | gzip -cd) \
//     --output starcode_umi_clusters.txt \
//     --sphere \
//     --seq-id \
//     --dist $l_distance \
//     2> starcode_umi_clusters.log

//     gzip starcode_umi_clusters.txt
    

//     #Demultiplex combined fastq file

//     python $scripts/demux_index_optim.py \
//     --fastq_path comb_reads \
//     --barcode_path starcode_umi_clusters.txt.gz \
//     --sample_name 'combined' \
//     --output_dir './'
//     """

//  }

/*
 * Create barcode sequence clusters
 * Default algorithm with Levenshtein distance=2
 * Use seq-id in output to demultiplex input files based on barcode
 */

// cluster 

process STARCODE_CLUSTERING { 
    tag "Creating global clusters from UMI barcodes"
    publishDir "${params.outdir}/clustering/starcode", mode:'copy'

    input: 
    tuple val(sample_id), path(reads)
    val l_distance

    output:
    tuple val(sample_id), path("*clusters.txt.gz"), emit: clusters_file
    path("*clusters.log")

    script:
    """
    starcode \
    -i <(cat ${reads} | gzip -cd) \
    --output ${sample_id}_starcode_umi_clusters.txt \
    --seq-id \
    --dist ${l_distance} \
    2> ${sample_id}_starcode_umi_clusters.log

    gzip ${sample_id}_starcode_umi_clusters.txt
    """
 }

/*
 * Demultiplex samples fastq and assign reads to barcode groups generated using starcode clustering
 * output should be 1 file per sample per cluster: ${sample}_${cluster}.fastq
 */

process DEMULTIPLEX_BARCODES {
    tag "Demultiplex $sample_id fastq files"
    publishDir "${params.outdir}/clustering/reads/$sample_id", mode:'symlink'

    input:
    tuple val(sample_id), path(reads), path(clusters_file)
    path scripts

    output:
    tuple val(sample_id), path("*.demux.fastq.gz"), emit: demux_reads
    path("*demux.log")

    script:
    """
    python ${scripts}/demux_index_optim_v2.py \
    --fastq_path ${reads} \
    --barcode_path ${clusters_file} \
    --sample_name ${sample_id} \
    --output_dir './'
    > ${sample_id}_demux.log
    """
}

/*
 * Align each sequence to its reference using sensitive and rapid BWA-MEM aligner
 * Assigning read groups with -R indicating $barcode & $sample origin origin.
 * Can use this info to manipulate sequences downstream
 */

process BWA_MEM_ALIGN {
    tag"Aligning $sample_id fastqs against $index"
    publishDir "${params.outdir}/alignment/bwa-mem/$sample_id", mode:'symlink'

    input:
    path index
    tuple val(sample_id), val(barcode), path(reads)

    output:
    tuple val(sample_id), val(barcode), path("*.sorted.bam"), emit: sorted_bam
    path("*.bwa.err")

    script:
    """  
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'` #find indexed ref files and strip suffix

    bwa mem -R '@RG\\tID:"${barcode}"\\tSM:${sample_id}\\tPL:Illumina' \$INDEX $reads 2> ${sample_id}_${barcode}.bwa.err \
        | samtools sort -O bam -o ${sample_id}_${barcode}.sorted.bam
    """

}

/*
 * Build BAM index for alignment visualisation & collect mapping statistics
 */

//process SAMTOOLS_INDEX_FLAGSTAT {
    // tag "Indexing $sample_id bam files & collecting alignment statistics"
    // publishDir "${params.outdir}/alignment/samtools/$sample_id", mode:'symlink'


    // input:
    // tuple val(sample_id), val(barcode), path(bam) 

    // output:
    // tuple val(sample_id), val(barcode), path("index/*.bam.bai"), emit:bam_indx
    // tuple val(sample_id), val(barcode), path("flagstat/*flagstat.out") 


    // script:
    // """
    // mkdir -p ./index
    // samtools index $bam > ./index/${sample_id}_${barcode}.sorted.bam.bai

    // mkdir -p ./flagstat
    // samtools flagstat $bam > ./flagstat/${sample_id}_${barcode}.flagstat.out
    // """

//}

/*
 * Variant Calling
 * Using  a combination of freebayes & bcftools
 * trade-off between sensitivity (bcftools) and specificity (freebayes)
 * min-alt-count just for testing for now for testing; disable for final setup (also, targets does not seem to be working..)
 *   freebayes --ploidy 2 --targets $bed --min-alternate-fraction 0.5 --min-alternate-count 1 --min-mapping-quality 1  --min-base-quality 3 --use-best-n-alleles=1 -f $index $bam > ${sample_id}_${barcode}.freebayes.vcf
*       dont worry about filtering as performed downstream using bcftools
 */

process FREEBAYES {
  tag "Calling variants on  $sample_id $barcode file"
  publishDir "${params.outdir}/calling/${params.caller}/$sample_id", mode:'symlink'

  input:
  path index
  path bed
  tuple val(sample_id), val(barcode), path(bam)   

  output:
  tuple val(sample_id), val(barcode), path("*.vcf"), emit: vcf


  script:
  """
  freebayes --min-alternate-count 1 --min-mapping-quality 1  --min-base-quality 3 --use-best-n-alleles 1 -f $index $bam > ${sample_id}_${barcode}.vcf
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
  tag "Calling variants on demultiplexed files"
  publishDir "${params.outdir}/calling/${params.caller}/$sample_id", mode:'copy'

  input:
  path index
  tuple val(sample_id), val(barcode), path(bam)   

  output:
  tuple val(sample_id), val(barcode), path("*.bcf"), emit: vcf
  

  script:
  """
  bcftools mpileup -I -Ou -f $index $bam | bcftools call --ploidy 2 --skip-variants indels --consensus-caller --variants-only -Ou -o ${sample_id}_${barcode}.bcf
  """
}

//bcftools variant calling options
// view:
// filter variants outside the target sequence (tat) & those with low supporting information: require at least 60% of reads to support call
// annotate vcf id column with sample_barcode info 
// normalise the variant calls (lelft-align & max parisomony) 
// output is compressed bcf for feeding into concat input

process BCFTOOLS_VIEW_ANNOTATE_NORM_INDEX {
  tag "Removing off-target & low-support variants. Normalizing variants in files"
  publishDir "${params.outdir}/annotation/bcftools/norm/${params.caller}/$sample_id", mode:'copy'

  input:
  path fasta
  path bed
  tuple val(sample_id), val(barcode), path(vcf) 

  output:
  tuple val(sample_id), path("*.norm.bcf"), emit: vcf
  tuple val(sample_id), path("*.csi"), emit: idx
  

  script:
  """
  bcftools view --targets "AF324493.2:5829-6044,AF324493.2:8368-8414" -q 0.6:nref $vcf -Ou | bcftools annotate --set-id "${sample_id}_${barcode}" |  bcftools norm --fasta-ref $fasta -Ob -o ${sample_id}_${barcode}.norm.bcf
  #create tbi index of file
  bcftools index ${sample_id}_${barcode}.norm.bcf
  """
}

// bcftools concat to concatenate vcfs per sample
// vcf files

process BCFTOOLS_CONCAT {
    tag " Concatenating sample $sample_id vcfs"
    publishDir "${params.outdir}/annotation/bcftools/concat/${params.caller}/$sample_id", mode:'copy'

    input:
    tuple val(sample_id), path(vcfs), path(idxs)

    output:
    tuple val(sample_id), path("*.combined.vcf"), emit: vcf

    script:
    """
    # create list of input files
    find -L ./ -name "*.norm.bcf" > vcfs_list
    bcftools concat --allow-overlaps --file-list ./vcfs_list --output ${sample_id}.norm.combined.vcf
    """
}


/*
 * Annotate variants in VCF using SNPEffect 
 * HGVS.c and HGVS.p annotation for codon and aa substitutions
 */

process SNPEFF_ANNO {
    tag "Annotating variants in $sample_id $barcode"
    publishDir "${params.outdir}/annotation/snpeff/${params.caller}/$sample_id", mode:'copy'

    input:
    path snpeff_db
    path snpeff_config //path snpeff_config make optional as already available through db path
    path intervals_bed
    val reference_v
    tuple val(sample_id), val(barcode), path(vcf)

    output:
    tuple val(sample_id), val(barcode), path("*.ann.vcf"), emit: ann_vcf

    // using shell block; nf variables assigned !, bash variables $
    //    # ggrep installed on mac to replicate linux versions
    // slower anyway as would perform search for each process; just specify string at CL
    //    #VERSION=`find -L ./data/ -type d | ggrep -E '.*/[A-Z]{1,2}[0-9]{5,6}[.][0-9]{1}?' | xargs -I {} basename {}`
    //#echo $VERSION

    //no_stats as want to prevent snpeff filtering my variants
    // other option emit the data dir and config file as seperate channels and then can easily specify
    """  
    snpEff ann -nodownload -v ${reference_v} \\
     -c ${snpeff_config} \\
     -csvStats -hgvs \\
     -noStats \\
     ${vcf} > ${sample_id}_${barcode}.ann.vcf
    """
}

///*    // #GENOME_V=`find . -type d | grep -E '.*/[A-Z]{1,2}[0-9]{5,6}\.[0-9]?' | xargs -I {} basename {}`

 /* To verify variant calls for a specific codon, we compared each barcode group among all sample VCFs. (Create channels of overlapping barcodes)
 * The minor fraction of variant calls in a particular barcode group that did not agree with a designed codon substitution was filtered and discarded.
 *
 */

//process BCFTOOLS_
//     -dataDir $snpeff_db \
//      -filterInterval ${intervals_bed}/intervals.bed \

/*
 * Annotate variants in VCF using VEP
 * HGVS.c and HGVS.p annotation for codon and aa substitutions
 */

// build SNPEff anotation dataabse given reference in genbank format

// process SNPEFF_ANN {
//     tag "Building SNPeff Annotation DB for $reference gbk file"
//     publishDir "${params.outdir}/annotation/snpeff/db", mode:'copy'

//    input:
//    path reference
//    tuple val(sample_id), val(barcode), path(bam)   

//    output:
//    tuple val(sample_id), val(barcode), path("mpileup/*.bcftools.raw.bcf")
//    tuple val(sample_id), val(barcode), path("call/*.bcftools.vcf"), emit: vcf
  


// }

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
    snpeffdb_ch      = SNPEFF_BUILD(reference_gbk_ch, scripts_ch)

    multiqc_input_ch = FASTQC_RAW(read_input_ch)
    merge_reads_ch   = BBMERGE(read_input_ch)
    trimmed_reads_ch = CUTADAPT_TRIM(merge_reads_ch.merged_reads, flank_5_ch, flank_3_ch)
    //trimmed_reads_ch = BBDUK_ADAPTER(merge_reads_ch.merged_reads, adapter_ch)
    //filtered_read_ch = REFORMAT_FILTER(trimmed_reads_ch.trimmed_reads)
    umi_reads_ch     = CUTADAPT_UMI(trimmed_reads_ch.trimmed_reads, umi_5_ch)
    //umi_h_reads_ch   = UMITOOLS_EXTRACT(trimmed_reads_ch.trimmed_reads, umi_5_ch)
    //umi_h_reads_ch   = UMITOOLS_EXTRACT(trimmed_reads_ch.trimmed_reads, umi_5_ch, umi_3_ch)
    // cluster umi for each sample
    clustered_umi_ch = STARCODE_CLUSTERING(umi_reads_ch.umi_reads, lev_dist_ch)

    //combine input channels by sample keys
    demulti_ch = DEMULTIPLEX_BARCODES(trimmed_reads_ch.trimmed_reads.combine(clustered_umi_ch.clusters_file, by:0), scripts_ch)

    //demulti_ch.demux_reads.view()

    // grouped_by_sample_barcode = demulti_ch.demux_reads.map { file_path ->
    //     def matcher = (file_path =~ /.*\/(.*)_(.*)\.demux\.fastq\.gz/)
    //     def (barcode, sample) = matcher ? [matcher[0][1], matcher[0][2]] : null
    // return tuple("${sample}_${barcode}", file_path)
    // }

    // grouped_by_sample_barcode.view()


    grouped_by_sample_barcode = demulti_ch.demux_reads.flatMap { sample, files -> files.collect { [sample, it] } } //flatMap transform ch into flattened [sample, file_path] tuples chs.
        .map { sample, file_path -> //create new tuple using map
        def matcher = (file_path =~ /.*\/(.*)_(.*)\.demux\.fastq\.gz/) //capture both barcode & sample groups
        def (barcode, sample_id) = matcher ? [matcher[0][1], matcher[0][2]] : null
    //return tuple("${sample_id}_${barcode}", file_path)
    return tuple(sample_id, barcode, file_path)
    }
    
    bam_ch = BWA_MEM_ALIGN(index_ch,grouped_by_sample_barcode)
    //samtools_idx_ch = SAMTOOLS_INDEX_FLAGSTAT(bam_ch.sorted_bam)

    // calling; depends on the caller specified

    if (params.caller == 'freebayes') {
        vcf_ch = FREEBAYES(reference_ch, bed_ch, bam_ch.sorted_bam)
    } else if (params.caller == 'bcftools') {
        vcf_ch = BCFTOOLS_MPILEUP_CALL(reference_ch, bam_ch.sorted_bam)
    }

    // filter off-target variants those with lower read support. annotate ID col with sample_barcode and normalise variants (left align, max parisomony representation) 
    norm_vcf_ch = BCFTOOLS_VIEW_ANNOTATE_NORM_INDEX(reference_ch, bed_ch, vcf_ch.vcf)

    
    // collapse vcf & indx files into a single channel per sample
    // tuple w 3 values: samplename, vcf, idx
    comb_ch = norm_vcf_ch.vcf.groupTuple()
    .combine(norm_vcf_ch.idx.groupTuple(), by:0).view()

    // concatenate variants in vcf files
    vcf_concat_ch = BCFTOOLS_CONCAT(comb_ch)

    // annotate variants w snpEff
    //ann_ch = SNPEFF_ANNO(snpeffdb_ch.snpeff_db, snpeffdb_ch.snpeff_config, scripts_ch, reference_gbk_ch, norm_vcf_ch.vcf)

    // .map { sample, file_path ->
    //     def umi_id    = file.name.toString().tokenize('_').get(1).tokenize('.demux.fastq.gz').get(0)
    //     def sample_umi_id = sample + '_' + umi_id
    // return tuple(sample_umi_id, file) }.view()

    //     read_pairs_ch.view()
//     read_pairs_ch.collaps.map { sample_id,file -> 
//           def umi_id = file.toString().tokenize('/').get(9).tokenize('_').get(1).tokenize('.').get(0) //use barcode as sample identifier will need to change when running on larger system
//         //  def single_file = file.filter( ~/^\/.*\.fastq\.gz$/ )
//           return tuple(sample_id, umi_id, file) }.view()

    // this command was used to combine umi files to generate clusters.. wont work with retrieving index
    //clustered_umi_ch = STARCODE_CLUSTERING(umi_reads_ch.umi_reads.flatten().collect(), lev_dist_ch)
    //umi_reads_ch.umi_reads.flatten().filter( ~/^\/.*\.fastq\.gz$/ ).collect().view()
    //clustered_umi_ch = STARCODE_CLUSTERING(umi_reads_ch.umi_reads.flatten().filter( ~/^\/.*\.fastq\.gz$/ ).collect(), lev_dist_ch) //flattens tuple & extracts elements ending in *.fastq.gz
    
    //need to maintain ordering for demux process as depends on fq index position. Combine UMI & Read channel and then 'split'
    //comb_ch = umi_reads_ch.umi_reads.combine(trimmed_reads_ch.trimmed_reads, by:0)
    //comb_umi_ch = comb_ch.flatten().filter( ~/^\/.*umi\.fastq\.gz$/ ).collect()
    //comb_fq_ch = comb_ch.flatten().filter( ~/^\/.*trim\.fastq\.gz$/ ).collect()

    //combined clustering & demux
    //demulti_ch = CLUSTER_DEMULTIPLEX(comb_umi_ch, comb_fq_ch, lev_dist_ch, scripts_ch) 

    // demulti_ch.demux_reads.map { sample_id,file -> 
    //     def umi_id = file.toString().tokenize('/').get(9).tokenize('_').get(0) //use barcode as sample identifier will need to change when running on larger system
    //     return tuple(sample_id, umi_id, file) }.groupTuple()// .groupTuple() would be a way to merge files by cell types

    // demulti_ch.demux_reads.map { sample_id,file -> 
    //      def umi_id = file.toString().tokenize('/').get(9).tokenize('_').get(0) //use barcode as sample identifier will need to change when running on larger system
    //      return tuple(sample_id, umi_id, file) }.groupTuple().view()

    // bam_ch = BWA_MEM_ALIGN(index_ch.index, 
    //     demulti_ch.demux_reads.map { sample_id,file -> 
    //     def umi_id = file.toString().tokenize('/').get(9).tokenize('_').get(0) //use barcode as sample identifier will need to change when running on larger system
    //     return tuple(sample_id, umi_id, file) }.groupTuple(),
    //     scripts_ch)
    //clustered_umi_ch = STARCODE_CLUSTERING(umi_reads_ch.umi_reads, lev_dist_ch)
    
    //demux_reads_ch   = DEMULTIPLEX_BARCODES(trimmed_reads_ch.trimmed_reads,clustered_umi_ch.clusters_file,scripts_ch) 
    
    //umi_h_reads_ch   = UMITOOLS_EXTRACT(merge_reads_ch.merged_reads, umi_5_ch, umi_3_ch)
    //clusters_ch      = STARCODE_CLUSTERING(umi_reads_ch.)
    MULTIQC_RAW(multiqc_input_ch.collect())

//     read_pairs_ch.view()
//     read_pairs_ch.collapsmap { sample_id,file -> 
//           def umi_id = file.toString().tokenize('/').get(9).tokenize('_').get(1).tokenize('.').get(0) //use barcode as sample identifier will need to change when running on larger system
//         //  def single_file = file.filter( ~/^\/.*\.fastq\.gz$/ )
//           return tuple(sample_id, umi_id, file) }.view()

}