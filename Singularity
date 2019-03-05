Bootstrap:docker
From:nfcore/base

%labels
	MAINTAINER Björn Hallström <bjorn.hallstrom@skane.se>
	DESCRIPTION Singularity container for CMD twist pipeline
	VERSION 0.0.1

%environment
	PATH=/opt/conda/envs/CMD-twist/bin:$PATH
	PICARD_HOME=/opt/conda/envs/CMD-twist/share/picard-2.18.26-0/


%files
        environment.yml /
        /data/bnf/scripts/postaln_qc.pl /usr/local/bin
	

%post
	/opt/conda/bin/conda env create -f /environment.yml
	/opt/conda/bin/conda clean -a
