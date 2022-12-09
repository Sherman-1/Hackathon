echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;103;30m|   Initialisation du pipeline ATIA - BOSSUT - HERMAN                           |\e[0m"
echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"

echo ""


echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;103;30m|   1 : Verification des pré-requis                                             |\e[0m"
echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"

if ! command -v nextflow &> /dev/null
then
    echo -e "Nextflow n'a pas été trouvé, verifiez que les pré-requis sont bien installés et relancez :)"
    exit
fi

if ! command -v singularity &> /dev/null
then
    echo -e "Singularity n'a pas été trouvé, verifiez que les pré-requis sont bien installés et relancez :)"
    exit
fi

echo -e "OK"

echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;103;30m|   2 : Creation des images                                             |\e[0m"
echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"

cd containers

cd featurecounts
echo -e "FeatureCounts:"
singularity build featcount.sif Singularity.featcount --fakeroot
cd ..

cd samtools
echo -e "Samtools:"
singularity build samtools.sif Singularity.samtools --fakeroot
cd ..

cd sratoolkit
echo -e "SRAtoolkit:"
singularity build sra.sif Singularity.sratoolkit --fakeroot
cd ..

cd Star
echo -e "STAR:"
singularity build star.sif Singularity.star_nb --fakeroot
cd ..


echo "Vive les DOGGOS"