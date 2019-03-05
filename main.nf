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
  .into { beds_mutect; beds_freebayes }



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
  set file("markdup.bam"), file("markdup.bam.bai") into bam_freebayes, bam_mutect

  """
    sambamba markdup --tmpdir /data/tmp -t 6 $sorted_bam markdup.bam
  """
}

process freebayes {
  input:
  set file(sorted_bam), file(sorted_bam_bai) from bam_freebayes
  file bed from beds_freebayes

  output:
  set val("freebayes"), file("freebayes_${bed}.vcf") into vcf_parts_freebayes

  """
    freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $sorted_bam > freebayes_${bed}.vcf
  """
}

process mutect {
  input:
  set file(sorted_bam), file(sorted_bam_bai) from bam_mutect
  file bed from beds_mutect

  output:
  set val("mutect"), file("mutect_${bed}.vcf") into vcf_parts_mutect

  """
    gatk --java-options "-Xmx2g" Mutect2 -R $genome_file -I $sorted_bam -L $bed -O mutect_${bed}.vcf

  """
}


process merge_freebayes_vcfs {
  input:
  set val(vc), file(vcfs) from vcf_parts_freebayes.groupTuple()

  output:
  file "${vc}.vcf.gz" into merged_vcf_freebayes

  """
    vcf-concat $vcfs | vcf-sort | gzip -c > ${vc}.vcf.gz
  """
}

process merge_mutect_vcfs {
  tag "mutect"
  input:
  set val(vc), file(vcfs) from vcf_parts_mutect.groupTuple()

  output:
  file "${vc}.vcf.gz" into merged_vcf_mutect

  """
    vcf-concat $vcfs | vcf-sort | gzip -c > ${vc}.vcf.gz
  """
}


process aggregate_vcfs {
  input:
  file(vcf1) from merged_vcf_mutect
  file(vcf2) from merged_vcf_freebayes

  output:
    file 'all.vcf' into result

  """
    cat $vcf1 $vcf2 > 'all.vcf'
  """
  
}

