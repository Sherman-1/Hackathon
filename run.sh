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
echo -e "\e[1;103;30m|   2 : Creation des images                                                     |\e[0m"
echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"

cd containers

cd featurecounts
echo -e "FeatureCounts:"
singularity build --fakeroot featcount.sif Singularity.featcount 
cd ..

cd samtools
echo -e "Samtools:"
singularity build --fakeroot samtools.sif Singularity.samtools 
cd ..

cd sratoolkit
echo -e "SRAtoolkit:"
singularity build --fakeroot sra.sif Singularity.sratoolkit 
cd ..

cd star
echo -e "STAR:"
singularity build --fakeroot star.sif Singularity.star_nb 
cd ..


echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;103;30m|   Terminé!                                                                    |\e[0m"
echo -e "\e[1;103;30m---------------------------------------------------------------------------------\e[0m"


echo -e "VOs images ont bien été crées, vous pouvez lancer le pipeline!"
