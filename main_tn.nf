#!/usr/bin/env nextflow

genome_file = file(params.fasta)
regions_bed = file(params.bed)
name        = params.name
threads     = 10

fastq = Channel.from( [['T', file(params.fastq_T_R1), file(params.fastq_T_R2)],
		       ['N', file(params.fastq_N_R1), file(params.fastq_N_R2)]] )

if(params.fasta ){
    bwaId = Channel
            .fromPath("${params.fasta}.bwt")
            .ifEmpty { exit 1, "BWA index not found: ${params.fasta}.bwt" }
}	     

Channel
    .fromPath("${params.bed}")
    .ifEmpty { exit 1, "Regions bed file not found: ${params.bed}" }
    .splitText( by: 1000, file: 'bedpart.bed' )
    .into { beds_mutect; beds_freebayes; beds_tnscope }



process bwa_align {

  input: 
	set val(type), file(r1), file(r2) from fastq

  output:
	set val(type), file("${type}_bwa.sort.bam") into bwa_bam

  """
     bwa mem -R '@RG\\tID:${name}_${type}\\tSM:${name}_${type}\\tPL:illumina' -M -t $threads $genome_file $r1 $r2 | samtools view -Sb - | samtools sort -o ${type}_bwa.sort.bam -
     samtools index ${type}_bwa.sort.bam
  """
}

process markdup {
  input:
  set val(type), file(sorted_bam) from bwa_bam

  output:
	set val(type), file("${type}.markdup.bam"), file("${type}.markdup.bam.bai") into bams

  """
    sambamba markdup --tmpdir /data/tmp -t 6 $sorted_bam ${type}.markdup.bam
  """
}

// Split tumor and normal bams into different channels
bamT = Channel.create()
bamN = Channel.create()
bams.choice(bamT, bamN) {it[0] == "T" ? 0 : 1}

// Send them to the different variant callers
(bamT_freebayes, bamT_mutect, bamT_tnscope) = bamT.into(3)
(bamN_freebayes, bamN_mutect, bamN_tnscope) = bamN.into(3)



process freebayes {
    input:
	set val(typeT), file(bamT), file(baiT) from bamT_freebayes
        set val(typeN), file(bamN), file(baiN) from bamN_freebayes
        each file(bed) from beds_freebayes

    output:
	set val("freebayes"), file("freebayes_${bed}.vcf") into vcfparts_freebayes

    when:
	params.freebayes
    
    """
    freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bamT $bamN > freebayes_${bed}.vcf
    """
}



process mutect {
    input:
	set val(typeT), file(bamT), file(baiT) from bamT_mutect
        set val(typeN), file(bamN), file(baiN) from bamN_mutect
        each file(bed) from beds_mutect

    output:
	set val("mutect"), file("mutect_${bed}.vcf") into vcfparts_mutect

    
    when:
	params.mutect
    
    """
    gatk --java-options "-Xmx2g" Mutect2 -R $genome_file -I $bamT -I $bamN -tumor ${name}_T -normal ${name}_N -L $bed -O mutect_${bed}.vcf
    """
}



process sentieon_preprocess_bam {
    cpus 8

    input:
	set val(typeT), file(bamT), file(baiT) from bamT_tnscope
        set val(typeN), file(bamN), file(baiN) from bamN_tnscope

    output:
	set file(bamT), file(baiT), file("T_recal.table") into processed_bamT_tnscope
	set file(bamN), file(baiN), file("N_recal.table") into processed_bamN_tnscope

    when:
	params.tnscope
    
    """
    /opt/sentieon-genomics-201808.01/bin/sentieon driver -t ${task.cpus} -r $genome_file -i ${bamT} --algo QualCal T_recal.table
    /opt/sentieon-genomics-201808.01/bin/sentieon driver -t ${task.cpus} -r $genome_file -i ${bamN} --algo QualCal N_recal.table
    """
}

process sentieon_tnscope {
  input:
    set file(bamT), file(baiT), file(recal_tableT) from processed_bamT_tnscope
    set file(bamN), file(baiN), file(recal_tableN) from processed_bamN_tnscope
    each file(bed) from beds_tnscope

  output:
  set val("tnscope"), file("tnscope_${bed}.vcf") into vcfparts_tnscope

  """
    /opt/sentieon-genomics-201808.01/bin/sentieon driver -t ${task.cpus} -r $genome_file -i $bamT -q $recal_tableT -i $bamN -q $recal_tableN --algo TNscope --tumor_sample ${name}_T --normal_sample ${name}_N --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.00005 tnscope_${bed}.tmp.vcf
    /opt/sentieon-genomics-201808.01/bin/sentieon driver -t ${task.cpus} -r $genome_file --algo TNModelApply --model /data/bnf/ref/sentieon/Sentieon_GiAB_HighAF_LowFP_201711.05.model -v tnscope_${bed}.tmp.vcf tnscope_${bed}.tmp2.vcf
    bcftools filter -s "ML_FAIL" -i "INFO/ML_PROB > 0.81" tnscope_${bed}.tmp2.vcf -m x -o tnscope_${bed}.vcf
  """  
}


// Prepare vcf parts for concatenation
vcfparts_freebayes = vcfparts_freebayes.groupTuple()
vcfparts_tnscope   = vcfparts_tnscope.groupTuple()
vcfparts_mutect    = vcfparts_mutect.groupTuple()
vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_mutect, vcfparts_tnscope)

process concatenate_vcfs {
    input:
	set vc, file(vcfs) from vcfs_to_concat

    output:
	set val("sample"), file("${vc}.vcf.gz") into concatenated_vcfs

    """
    vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
    vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
    vt normalize ${vc}.decomposed.vcf.gz -r $genome_file | vt uniq - -o ${vc}.vcf.gz
    """
}



process aggregate_vcfs {
    input:
        set sample, file(vcfs) from concatenated_vcfs.groupTuple()
    
    output:
	file 'all.vcf.gz' into result

    """
    cat $vcfs > all.vcf.gz
    """
  
}
