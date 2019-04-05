/data/bnf/sw/nextflow/nextflow run main.nf \
			       --fasta /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
			       --read1 /data/bnf/dev/bjorn/nextflow_test/data/supersmall_R1_001.fastq.gz \
			       --read2 /data/bnf/dev/bjorn/nextflow_test/data/supersmall_R2_001.fastq.gz \
			       --name TEST \
			       --bed /data/bnf/proj/twist/myeloid_panel/final_design/LundUni_NGSTECustom_0001506/all_target_segments_covered_by_probes_LundUni-Targets+Genes_NGSTECustom_0001506c_hg19.bed \
			       -with-singularity container_2019-04-03.sif \
			       -with-dag test.dag.png
 
