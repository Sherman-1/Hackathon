Bootstrap: docker
From: ubuntu:20.04



%environment 

	export LC_ALL=C





%post 

	# Install compiling libraries and R dependencies
	apt-get update --fix-missing -qq 
	apt-get install -y --no-install-recommends software-properties-common dirmngr
	apt-get install -y wget tzdata build-essential libssl-dev libz-dev libpng-dev libblas-dev \
       			liblapack-dev libcurl4-openssl-dev libxml2-dev pkg-config gfortran

	# Get latest R version as shown on CRAN website
	wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc |  tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
	add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" 
	apt-get install -y --no-install-recommends r-base 
	
	
	# Install DE analysis tools
	R -e "install.packages('BiocManager')"
	R -e "BiocManager::install('DESeq2')" 

%test 


	R -e "1+1"
	R -e "library('ggplot2')"


