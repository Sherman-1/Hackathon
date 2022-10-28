nextflow.enable.dsl=2

process extractFastq {
    input:
        val Srr

    script:
    """
    prefetch $Srr && fastq-dump $Srr
    """
}

process extractgenome {

    output:
        path 'ref_genome'

    script:
    """
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz
    gunzip -c chr3.fa.gz > ref.fa  
    """
}

process extractGFF{

    output:
        path annot

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    """
}

process index{
    input:
        path ref, path annot

    output:
        path index_genome

    script:
    """
    STAR --runMode genomeGenerate --runThreadN 4 \
     --genomeDir ./index_genome \
     --genomeFastaFiles $ref \
     --sjdbGTFfile $annot
    """
}

process mapping{
    input:
        path fastqfile, path index_genome

    script:
    """
    STAR --runThreadN 4 --outFilterMultimapNmax 10  --outFilterMismatchMax 4\
       --genomeDir $index_genome \
       --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix '_mapping' \
       --readFilesIn $fastqfile
  
    """

}

workflow{
    sample = channel.of("SRR628589","SRR628588")
    extractFastq(sample)
    genome = extractgenome()
    annot = extractGFF()
    index(genome,annot)
}