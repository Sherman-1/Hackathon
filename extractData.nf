#!/usr/bin/env nextflow
nextflow.enable.dsl=2

sample = Channel.of("SRR628589","SRR628588")

process extractFastq {
    /* blabla extrait des samples  blabla... */
    input:
    val Srr from sample

    output:
    tuple val Srr, file("*1.fastq.gz"), file("*2.fastq.gz") into ch_fastq

    script:
    """
    prefetch $Srr && fastq-dump --gzip --split-file $Srr
    """
}

chr = Channel.of('1','2' ,'3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')

process extractAllGenome {
    /* chaque chr */

    input:
    val chromosome from chr

    output:
    file "${chromosome}.fa.gz" into ch_chr
    script:
    """
    wget -O ${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    """
}

process assembleGenome {
    /* on fait en sorte que le génome soit dans le même fichier et pas en 24 trucs" */

    input:
    file chromosomes from ch_chr.collect()

    output:
    file "genome.fa" into ch_genome

    script:
    """
    gunzip -c ${chromosomes} > genome.fa
    """
}


process index{
    /* index tu connais */
    input:
    file (genome) from ch_genome.collect()

    output:
        path "ref" into ch_ref

    script:
    """
    STAR --runMode genomeGenerate --runThreadN 4 \
     --genomeDir ./index_genome \
     --genomeFastaFiles $genome \
    """
}

process extractGFF{
    /* extract gff c explicit */

    output:
    file "annot.gtf" into ch_annot

    script:
    """
    wget -O annot.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip -c annot.gtf.gz > annot.gtf
    """
}


process mapping{
    /* On mappe les zamis*/
    input:
        val Srr, file(r1), file (r2) from ch_fastq
        path index from ch_ref
    
    output:
    file "${Srr}.bam" into ch_bam

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
    > ${Srr}.bam
    """

}

process samtoolsIndex{
    /* oui */

    input:
    file bam from ch_bam

    output:
    file "${bam}.bai"

    script:
    """
    samtools index ${bam}
    """

}