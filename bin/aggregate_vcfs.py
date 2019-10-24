#!/usr/bin/env python2
import vcf
import sys
import argparse
import pprint

def main():
    opt = parse_arguments()

    aggregated_variants = {}

    for vcf_fn in opt.vcfs:
        vcf_reader = vcf.Reader(open(vcf_fn))

        variant_caller = which_variantcaller(vcf_reader.metadata)
        for var in vcf_reader:

            # Check if multi-allelic site
            if len(var.ALT) > 1:
                raise NameError('Split and normalize you variants!')

            calc_GT_fields(var, variant_caller)
                
            # Save variant in aggregated_variants if it hasn't been found before
            var_id = str(var.CHROM)+"_"+str(var.POS)+"_"+str(var.REF)+"_"+str(var.ALT[0])
            if not aggregated_variants.get(var_id):
                aggregated_variants[var_id] = var

            # Add variant caller information to an INFO field
            info_field = aggregated_variants[var_id].INFO
            if not info_field.get("variant_callers"):
                info_field["variant_callers"] = variant_caller
            else:
                info_field["variant_callers"] += "|"+variant_caller

                
    # Output aggregated vcf
    vcf_template = vcf.Reader(open(opt.vcfs[0]))    
    vcf_writer = vcf.Writer(open('outtest.vcf', 'w'), vcf_template)
    for key, record in aggregated_variants.items():
        vcf_writer.write_record(record)
        vcf_writer.flush()
                
    vcf_writer.close()

def parse_arguments():
    parser = argparse.ArgumentParser(description='Aggregate vcfs from different variant callers.')
    parser.add_argument('vcfs', metavar='vcfs', type=str, nargs='+', help='vcf files')

    opt = parser.parse_args()
    return opt


def which_variantcaller(meta):

    source  = meta.get("source", [])
    cmdline = meta.get("commandline", [])
    mutect2_cmdline = meta.get("GATKCommandLine.MuTect2", [])
    if any("freeBayes" in s for s in source) or any("freebayes" in s for s in cmdline):
        return "freebayes"

    if any("Mutect2" in s for s in source) or any("Mutect2" in s for s in cmdline) or mutect2_cmdline:
        return "mutect2"

    if meta.get("SentieonCommandLine.TNscope"):
        return "tnscope"
    
    return "unknown"


def calc_GT_fields(var, vc):
    if vc == "mutect2":
#        for c in var.samples:
            #c.add_field("VAF", c.AF)
#            print(c.CallData)
        for i in range(len(var.samples)):
            #var.samples[i] = var.samples[i]._replace(CallData = "a")
            #var.samples[i].data.gt_alleles = "1/1"
            var.samples[i].data._replace(GT="1/1")
            print(var.samples[i].data)
            print(var.samples[i])
            


if __name__ == '__main__':
    main()




