#DockerFile pour sra-toolkit
FROM ubuntu:22.04

RUN apt-get update

RUN apt-get install -y sra-toolkit

CMD which fastq-dump
