process extractFastq {

    """
    prefetch SRR636533 && fastq-dump SRR636533
    """
}

process extractgenome {
    """
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz
    gunzip -c chr3.fa.gz > ref.fa  
    """
}

process extractGFF{
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    """
}
