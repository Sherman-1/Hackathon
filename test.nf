#!/usr/bin/env nextflow

//Données externes

sample = Channel.of("SRR628588", "SRR628589")

chr = Channel.of("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16" , "17", "18", "19", "20", "21", "22", "MT", "X", "Y")

gtf_URL="ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz"

/*
process downloadFastq{
	publishDir "data/SRA/"

	input:
	val sraid from sample

	output:
	tuple val(sraid), file("${sraid}.sra") into ch_fastq1

	
	"""
    wget -O ${sraid}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/srapub/${sraid}
	"""
}

process fastqDump {
    publishDir "data/fastq/"
    
    input:
    tuple val(sraid), file("${sraid}.sra") from ch_fastq1
    
    output:
    tuple val(sraid), file("*1.fastq.gz"), file("*2.fastq.gz") into ch_fastq2
    
    
    """    
    fastq-dump --gzip --split-files ${sraid}.sra
    """
}
*/


process fastqDump {
    publishDir "data/fastq/"
    
    input:
    val sraid from sample
    
    output:
    tuple val(sraid), file("*1.fastq.gz"), file("*2.fastq.gz") into ch_fastq2
    
    
    """    
    fasterq-dump --split-files ${sraid}
    """
}

process extractAllGenome {
	publishDir "data/genome/"

    input:
    val chromosome from chr
    
    output:
    file "${chromosome}.fa.gz" into ch_chr
    
    
    """
    wget -O ${chromosome}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chromosome}.fa.gz
    """
}


process assembleGenome {
    publishDir "data/genome/"

    input:
    file '*.fa.gz' from ch_chr.collect()
		
    output: 
    file "genome.fa" into ch_genome
    
    script:
    """
    gunzip -c *.fa.gz > genome.fa
    """
}

process index {
    publishDir "results/genome_index/"
    
    input:
    file (genome) from ch_genome

    output:
    path "ref" into ch_ref

    script:
    """
    STAR --runMode genomeGenerate --runThreadN 4\
    --genomeDir ref/ \
    --genomeFastaFiles ${genome}
    """
}

process extractGFF{
    publishDir "data/genome_annotation/"
    
    input:
    val gtfURL from gtf_URL
    
    output:
    file "annot.gtf" into ch_annot
    
    script:
    """
    wget ${gtfURL}
    gunzip -c *.gtf.gz > annot.gtf
    """
}

process mapping {
    publishDir "data/bam/"

    input:
    tuple val(srr), file(r1), file(r2) from ch_fastq2
    path index from ch_ref
 
    output:
    file "${srr}.bam" into ch_bam
 
    script:
    """
    STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ${index} \
        --readFilesIn <(gunzip -c ${r1}) <(gunzip -c ${r2}) \
        --runThreadN 4 \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --outFileNamePrefix ${srr} 
        > ${srr}.bam    
    """
}

process samtoolsIndex {
    publishDir "data/bam_index/"

    input:
    file bam from ch_bam
 
    output:
    file("${bam}.bai")
 
    script:
    """
    samtools index ${bam} 
    """
}
