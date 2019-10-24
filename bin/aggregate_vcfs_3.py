#!/usr/bin/env python2
from cyvcf2 import VCF
import gzip
import sys
import argparse
import pprint

def main():
    opt = parse_arguments()

    aggregated_variants = {}

    for vcf_fn in opt.vcfs:

        vcf_reader = VCF(vcf_fn)
        vcf_reader.add_info_to_header({'ID': 'blaha', 'Description': 'aList of variant callers which detected the variant',
                                       'Type':'Character', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'variant_callers', 'Description': 'List of variant callers which detected the variant',
                                       'Type':'Character', 'Number': '1'})

        variant_caller = which_variantcaller(vcf_reader)

        for var in vcf_reader:

            # Check if multi-allelic site
            if len(var.ALT) > 1:
                raise NameError('Split and normalize you variants!')

            #calc_GT_fields(var, variant_caller)

            # Save variant in aggregated_variants if it hasn't been found before
            var_id = str(var.CHROM)+"_"+str(var.POS)+"_"+str(var.REF)+"_"+str(var.ALT[0])

            if not aggregated_variants.get(var_id):
                aggregated_variants[var_id] = var

            # Add variant caller information to an INFO field
            vcs = ""
            if not aggregated_variants[var_id].INFO.get("variant_callers"):
                vcs = variant_caller
                aggregated_variants[var_id].INFO["variant_callers"] = vcs
            #aggregated_variants[var_id].INFO["variant_callers"] = variant_caller
            else:
                print("INTHERE")
                vcs = "TEST"
                aggregated_variants[var_id].INFO["blaha"] = vcs
                #aggregated_variants[var_id].INFO["blaha"] = "multiple"



            print var
    # Output aggregated vcf
    #vcf_template = vcf.Reader(open(opt.vcfs[0]))
    #vcf_writer = vcf.Writer(open('outtest.vcf', 'w'), vcf_template)
    #for key, record in aggregated_variants.items():
    #    vcf_writer.write_record(record)
    #    vcf_writer.flush()
    #vcf_writer.close()

def parse_arguments():
    parser = argparse.ArgumentParser(description='Aggregate vcfs from different variant callers.')
    parser.add_argument('vcfs', metavar='vcfs', type=str, nargs='+', help='vcf files')

    opt = parser.parse_args()
    return opt


def which_variantcaller(vcf):
    for line in vcf.raw_header.split("\n"):
        if line.startswith("##"):
            if "source" in line and "mutect2" in line.lower():
                return "mutect2"
            if "source" in line and "freebayes" in line.lower():
                return "freebayes"
            if "SentieonCommandLine.TNscope" in line:
                return "tnscope"

        elif line.startswith("#"):
            return "unknown"
    
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




