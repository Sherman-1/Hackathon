# **Projet Repro-Hackathon 2022-2023** : ATIA Safiya, BOSSUT Noémie et HERMAN Simon

Ce projet vise à reproduire une partie des résultats de deux articles:
- [Furney *et al.*, Cancer Discovery (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321577/)
- [Harbour *et al.*, Nature Genetics (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789378/)

Les données RNA-seq de ces deux papiers sont disponibles en open-access : [**Données NCBI**](https://www.ncbi.nlm.nih.gov/sra?term=SRA062359). Dans un premier temps, seul l'étude de transcriptome est étudiée. </br> 
L'objectif de ces deux articles est d'étudier les expression de gènes, et notamment le gène SF3B1, d'individus atteint de mélanome uvéal. Harbour *et al.* affirment que le gène SF3B1 est muté, pour des données exomiques, mais est généralement présente dans des tumeurs bénignes avec un bon taux de survie des patients.  Furney *et al.* montrent que, en réutilisant le même jeu de données qu'il existe une association entre la mutation SF3B1 et plusieurs types d'épissage alternatif : rétention d'introns, épissage sur sites cryptiques et epissage en site 3' alternatif, ce qui rendrait cette mutation potentiellement à l'origine du développement cancéreux. 

A l'aide d'un workflow Nextflow et de *containers* Singularity, notre groupe a tenté de comprendre pourquoi les résultats des deux articles divergent, et quelles sont nos propres observations sur le sujet. 

## **Pré-requis:** Nextflow & Singularity
Afin de faire tourner notre *pipeline*, **64 Gb de RAM et 16 coeurs**, ainsi que deux logiciels sont nécessaires:  
- Nextflow (version 21.10.6.5660) [*installation*](https://www.nextflow.io/docs/latest/getstarted.html)
- Singularity (version 3.8.7) [*installation*](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
- uidmap (pour la construction des containers Singularity)
  
```bash
    sudo apt-get install uidmap
```

## **Le pipeline:**
 ![image](https://user-images.githubusercontent.com/75751225/206907177-60bd1d6f-84ae-4c55-a9e2-80cc1f44a2a7.png)


1) **Téléchargement des données** : chromosomes humains (dont chromosome mitochondrial), annotation du génome et données RNA-seq des 8 individus (***sratoolkit***). 
2) **Indexation** du génome (***STAR***).
3) **Alignement** et **Tri** des données RNA-seq sur le génome. Obtention de fichiers *.bam* triés en sortie (***STAR***).
4) **Indexation** des fichiers *.bam*. en *.bai* (***samtools***). (optionnel)
5) **Comptage** des séquences exprimées (***featureCounts***).
6) **Analyse statistique** des résultats (***DESeq2***).

L'ensemble des données et des résultats peuvent être retrouvés dans l'arborescence ci-dessous:

*Avant execution du workflow:*
``` bash
.
├── bin
│   └── DE_analysis.R
├── containers
│   ├── featurecounts
│   │   └── Singularity.featcount
│   ├── samtools
│   │   └── Singularity.samtools
│   ├── Singularity.R
│   ├── sratoolkit
│   │   └── Singularity.sratoolkit
│   └── star
│       └── Singularity.star_nb
├── nextflow.config
├── README.md
├── run.sh
└── workflow.nf
```


## **Exécution du workflow**
Si les pré-requis sont bien satisfaits, placez-vous dans le répertoire voulu et récupérez le projet
```bash
    git clone https://github.com/Sherman-1/Hackathon
    cd Hackathon
```
Le fichier `run.sh` permet d'initialiser votre environnement, ainsi que de créer les images singularity
```bash
    bash run.sh
```
Vérifiez que les répertoires ont bien été créés. Si c'est le cas, vous pouvez lancer le pipeline :
```bash
    nextflow run workflow.nf
```
### Options d'execution
En plus de la commande par défaut, vous pouvez utiliser les paramètres suivants :
* `-resume` si votre *pipeline* a été interrompu et que vous souhaitez le reprendre
* `-with-trace [nom du fichier]` pour obtenir le DAG correspondant
* `-with-report [nom du fichier]` pour obtenir un rapport complet et de nombreuses metadatas sur le *pipeline*
