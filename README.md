# **Projet Repro-Hackathon 2022-2023** : ATIA Safiya, BOSSUT Noémie et HERMAN Simon

Ce projet vise à reproduire une partie des résultats de deux articles:
- [Furney *et al.*, Cancer Discovery (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321577/)
- [Harbour *et al.*, Nature Genetics (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789378/)

Les données RNA-seq de ces deux papiers sont disponibles en open-access : [**Données NCBI**](https://www.ncbi.nlm.nih.gov/sra?term=SRA062359). Dans un premier temps, seul l'étude de transcriptome est étudié. </br> 
L'objectif de ces deux article est d'étudier les expression de gènes, et notamment le gène SF3B1, d'individus atteint de mélanome ulvéal. 
# ECRIRE ICI CONCLU DES ARTICLES

A l'aide d'un workflow Nextflow et de containers Singularity, notre groupe a tenté de comprendre pourquoi les résultats des deux articles divergent, et quelles sont nos propres observations sur le sujet. 
# NOTRE CONCLU A NOUS

## **Pré-requis:** Nextflow & Singularity
Afin de faire tourner notre pipeline, 1000000Gb DE RAM ET 300 COEURS, ainsi que deux logiciels sont nécessaires:  
    - Nextflow (version 21.10.6.5660)
    - Singularity (version 3.8.7)