FROM ubuntu:22.04

# Dépendances de base
RUN apt-get update && apt-get install -y wget\
	software-properties-common\
    unzip\
    curl\
    gcc\
    build-essential\
    git

#Logiciels d'analyse BASH
RUN git clone https://github.com/alexdobin/STAR.git
RUN cd STAR/source\
    make STAR
RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
RUN tar -vxzf sratoolkit.tar.gz
RUN export PATH=$PATH:$PWD/sratoolkit.3.0.0-mac64/bin


# DESeq2 et FeaturesCount à installer