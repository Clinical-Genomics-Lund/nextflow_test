#!/usr/bin/env python3
import vcfpy
import sys
import argparse
import pprint

def main():
    opt = parse_arguments()

    aggregated_variants = {}

    for vcf_fn in opt.vcfs:
        vcf_reader = vcfpy.Reader.from_path(vcf_fn)

        variant_caller = which_variantcaller(vcf_reader.header)

        for var in vcf_reader:

            # Check if multi-allelic site
            if len(var.ALT) > 1:
                raise NameError('Split and normalize your variants!')

            #calc_GT_fields(var, variant_caller)
                
            # Save variant in aggregated_variants if it hasn't been found before
            var_id = str(var.CHROM)+"_"+str(var.POS)+"_"+str(var.REF)+"_"+str(var.ALT[0])
            if not aggregated_variants.get(var_id):
                aggregated_variants[var_id] = var

            #print(var)
            #print(var.keys)
            # Add variant caller information to an INFO field
            #info_field = aggregated_variants[var_id].INFO
            #if not info_field.get("variant_callers"):
            #    info_field["variant_callers"] = [variant_caller]
            #else:
            #    info_field["variant_callers"].append(variant_caller)

                
    # Output aggregated vcf
    vcf_template = vcfpy.Reader.from_path(opt.vcfs[0])
    vcf_writer = vcfpy.Writer.from_path('outtest.vcf', vcf_template.header)

    for key, record in aggregated_variants.items():
        vcf_writer.write_record(record)
                
#    vcf_writer.close()

def parse_arguments():
    parser = argparse.ArgumentParser(description='Aggregate vcfs from different variant callers.')
    parser.add_argument('vcfs', metavar='vcfs', type=str, nargs='+', help='vcf files')

    opt = parser.parse_args()
    return opt


def which_variantcaller(meta):

    for a in meta.lines:

        if a.key == "source":
            if "freeBayes" in a.value:
                return "freebayes"
            if "Mutect2" in a.value:
                return "mutect2"

        if a.key == "SentieonCommandLine.TNscope":
            return "tnscope"

    return "unkown"

#def calc_GT_fields(var, vc):
#    if vc == "mutect2":
#        for c in var.samples:
#            c.add_field("VAF", c.AF)
#            print c
            


if __name__ == '__main__':
    main()




