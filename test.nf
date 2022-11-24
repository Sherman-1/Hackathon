#!/usr/bin/env nextflow

//Données externes

sample = Channel.of("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588", "SRR628589")

chr = Channel.of("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16" , "17", "18", "19", "20", "21", "22", "MT", "X", "Y")

gtf_URL="ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz"

process fastqDump {
    publishDir "data/fastq/"
    
    input:
    val sraid from sample
    
    output:
    tuple val(sraid), file("*1.fastq"), file("*2.fastq") into ch_fastq2
    
    
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
    file (genome) from ch_genome.collect()

    output:
    path "ref" into ch_ref

    script:
    """
    STAR --runMode genomeGenerate --runThreadN 20\
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
    file "${srr}.bam" into(ch_bam, ch_count)
 
    script:
    """
    STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ${index} --readFilesIn ${r1} ${r2} --runThreadN 20 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --genomeLoad NoSharedMemory --outFileNamePrefix ${srr} > ${srr}.bam    
    """
}

process samtoolsIndex {
    publishDir "data/bam_index/"

    input:
    file bam from ch_bam
 
    output:
    file("${bam}.bai") into ch_samtools
 
    script:
    """
    samtools index ${bam} 
    """
}


process countingReads{
    publishDir "results/featureCounts/"

    input:
    file indexbam from ch_count.collect()
    file gtf from ch_annot

    output:
    file "matrice_featureCounts.txt"

    script:
    """
    featureCounts -T 20 -t gene -g gene_id -s 0 -a $gtf -o matrice_featureCounts.txt ${indexbam}

    """

}
