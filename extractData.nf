process extractData {

    """
    prefetch -v SRR636533 --force all
    fastq-dump --gzip --split-files ./SRR636533.sra

    prefetch -v SRR636567 
    fastq-dump --gzip --split-files ./SRR636533.sra
    """
}