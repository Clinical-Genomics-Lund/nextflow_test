#!/usr/bin/env nextflow

genome_file = file(params.fasta)
fq_read1    = file(params.read1)
fq_read2    = file(params.read2)
regions_bed = file(params.bed)
name        = params.name
threads     = 10

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
  file r1 from fq_read1
  file r2 from fq_read2

  output:
  set val(name), file("${name}_bwa.sort.bam") into bwa_bam

  """
     bwa mem -R '@RG\\tID:${name}\\tSM:${name}\\tPL:illumina' -M -t $threads $genome_file $r1 $r2 | samtools view -Sb - | samtools sort -o ${name}_bwa.sort.bam -
     samtools index ${name}_bwa.sort.bam
  """
}

process markdup {
  input:
  set val(name), file(sorted_bam) from bwa_bam

  output:
  set file("markdup.bam"), file("markdup.bam.bai") into bam_freebayes, bam_mutect, bam_tnscope

  """
    sambamba markdup --tmpdir /data/tmp -t 6 $sorted_bam markdup.bam
  """
}

process freebayes {
  input:
  set file(sorted_bam), file(sorted_bam_bai) from bam_freebayes
  file bed from beds_freebayes

  output:
  set val("freebayes"), file("freebayes_${bed}.vcf") into vcfparts_freebayes

  """
    freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $sorted_bam > freebayes_${bed}.vcf
  """
}

process mutect {
  input:
  set file(sorted_bam), file(sorted_bam_bai) from bam_mutect
  file bed from beds_mutect

  output:
  set val("mutect"), file("mutect_${bed}.vcf") into vcfparts_mutect

  """
    gatk --java-options "-Xmx2g" Mutect2 -R $genome_file -I $sorted_bam -L $bed -O mutect_${bed}.vcf
  """
}

process sentieon_preprocess_bam {
  cpus 8

  input:
	set file(sorted_bam), file(sorted_bam_bai) from bam_tnscope

  output:
	set file("realn.bam"), file("realn.bam.bai"), file("recal.table") into processed_bam_tnscope

  """
    /opt/sentieon-genomics-201808.01/bin/sentieon driver -t ${task.cpus} -r $genome_file -i $sorted_bam --algo Realigner realn.bam
    /opt/sentieon-genomics-201808.01/bin/sentieon driver -t ${task.cpus} -r $genome_file -i realn.bam --algo QualCal recal.table
  """
}


process sentieon_tnscope {
  input:
    set file(bam), file(bai), file(recal_table) from processed_bam_tnscope
    file bed from beds_tnscope

  output:
  set val("tnscope"), file("tnscope_${bed}.vcf") into vcfparts_tnscope

  """
    /opt/sentieon-genomics-201808.01/bin/sentieon driver -t ${task.cpus} -r $genome_file -i $bam -q $recal_table --algo TNscope --tumor_sample ${name} tnscope_${bed}.vcf
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
	set val("sample"), file("*.vcf.gz") into concatenated_vcfs

    """
    vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
    vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
    vt normalize ${vc}.decomposed.vcf.gz -r $genome_file | vt uniq - -o ${vc}.vcf.gz
    """
}




process aggregate_vcfs {
    input:
        set sample, file(vcfs) from concatenated_vcfs.toList()
    
    output:
	file 'all.vcf.gz' into result

    """
    cat $vcfs > all.vcf.gz
    """
  
}
