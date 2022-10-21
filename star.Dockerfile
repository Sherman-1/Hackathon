# DockerFile pour STAR
FROM ubuntu:22.04

RUN apt-get update --fix-missing
RUN apt-get install -y rna-star

#Commande juste pour vérifier que l'image tourne avec STAR
CMD STAR
