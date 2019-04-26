Bootstrap:docker
From:nfcore/base

%labels
	MAINTAINER Björn Hallström <bjorn.hallstrom@skane.se>
	DESCRIPTION Singularity container for CMD twist pipeline
	VERSION 0.0.1

%environment
	PATH=/opt/conda/envs/CMD-twist/bin:/opt/sentieon-genomics-201808.05/bin/:$PATH
	PICARD_HOME=/opt/conda/envs/CMD-twist/share/picard-2.18.26-0/


%files
        environment.yml /
        /data/bnf/scripts/postaln_qc.pl /usr/local/bin
	/data/bnf/sw/sentieon/sentieon-genomics-201808.05 /opt
        data/GenomeAnalysisTK-3.8.tar.bz2 /opt

%post
	/opt/conda/bin/conda env create -f /environment.yml
	/opt/conda/bin/conda clean -a
	#gatk3-register /opt/sentieon-genomics-201808.01