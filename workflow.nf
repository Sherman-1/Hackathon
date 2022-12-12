#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastqDump {
    publishDir "data/fastq/"

    input:
        val samples
    
    output:
        tuple val(samples), file("*1.fastq"), file("*2.fastq")
    
    """    
    fasterq-dump --split-files ${samples}
    """
}

process extractAllGenome {
    publishDir "data/genome/"

    input:
        val chromosome
        val genome_URL

    output:
        path "${chromosome}.fa.gz"  //ch_chr
    
    """
    wget -O ${chromosome}.fa.gz ${genome_URL}.${chromosome}.fa.gz
    """
}

process assembleGenome {
    publishDir "data/genome/"

    input:
    path fasta //file('*.fa.gz') from ch_chr.collect()
		
    output: 
    path "genome.fa" // ch_genome
    
    script:
    """
    gunzip -c ${fasta} > genome.fa
    """
}

process index {
    publishDir "results/genome_index/"
    cpus params.index_cpus

    input:
    path genome //from ch_genome.collect()

    output:
    path "ref" //into ch_ref

    script:
    """
    STAR --runMode genomeGenerate --runThreadN ${task.cpus}\
    --genomeDir ref/ \
    --genomeFastaFiles ${genome}
    """
}

process extractGFF{
    publishDir "data/genome_annotation/"
    
    input:
    val gtfURL //from gtf_URL
    
    output:
    path "annot.gtf" //into ch_annot
    
    script:
    """
    wget ${gtfURL}
    gunzip -c *.gtf.gz > annot.gtf
    """
}

process mapping {
    cpus params.mapping_cpus

    publishDir "data/bam/"

    input:
    tuple val(srr), file(r1), file(r2) //from ch_fastq2
    path index //from ch_ref
 
    output:
    path "${srr}.bam" //into(ch_bam, ch_count)
 
    script:
    """
    STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ${index} --readFilesIn ${r1} ${r2} --runThreadN ${task.cpus} --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --genomeLoad NoSharedMemory --outFileNamePrefix ${srr} > ${srr}.bam    
    """
}

process samtoolsIndex {
    publishDir "data/bam_index/"

    input:
    file bam //from ch_bam
 
    output:
    path "${bam}.bai" // into ch_samtools
 
    script:
    """
    samtools index ${bam} 
    """
}

process countingReads{
    cpus params.counting_cpus
    publishDir "results/featureCounts/"

    input:
    path indexbam //from ch_count.collect()
    path gtf // from ch_annot

    output:
    path "matrice_featureCounts.txt"

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a $gtf -o matrice_featureCounts.txt ${indexbam}

    """

}


process DESeq{
    publishDir "results/Figures/"

    input:
    path DESeq_script
    path count_matrix
    
    output:
    path Figures
    
    script:
    """
    mkdir Figures
    Rscript --vanilla ${DESeq_script} "${count_matrix}" 
    """

}




workflow {

    DESeq_script_path = Channel.fromPath("bin/DE_analysis.R")

    samples = Channel.fromList(params.samples)
    chromosome = Channel.fromList(params.chr)
    fastqDump(samples)
    extractAllGenome(chromosome, params.genome_URL)
    assembleGenome(extractAllGenome.out.collect())
    index(assembleGenome.out.collect())
    extractGFF(params.gtf_URL)
    mapping(fastqDump.out, index.out)
    // samtoolsIndex(mapping.out) 
    // samtoolsIndex not needed : BAI outputs
    countingReads(mapping.out.collect(), extractGFF.out)
    results_path = DESeq(DESeq_script_path, countingReads.out)
    results_path.view()

}
