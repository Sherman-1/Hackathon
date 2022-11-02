Bootstrap: docker
From: ubuntu:20.04


%post
	apt update -qq -y
	
	# Install R dependencies
	apt install --no-install-recommends -y software-properties-common dirmngr 
	apt install -y gnupg apt-transport-https ca-certificates 
	
	# Add CRAN to apt known repos
	apt-key adv --keyserver keyserver.ubuntu.com \
	--recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
	add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
	
	apt -y install r-base

	# Install DESeq2 dependencies
	apt install -y libxml2 libxml2-dev curl libssl-dev libcurl4-openssl-dev
