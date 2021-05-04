#!/usr/bin/env python


import pandas as pd
import argparse
import sys

def parse_args(args=None):
    Description = 'Check samplesheet and add files to bed file'
    Epilog = """Example usage: python modify_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet.")
    parser.add_argument('FILE_OUT', help="Output samplesheet.")
    return parser.parse_args(args)


def add_bed_file(FileIn,FileOut):
    #Open input file
    fi = open(FileIn, 'r')

    # Load mosdepth thresholds.bed.gz into a pandas dataframe
    ss = pd.read_csv(fi, delimiter=',', index_col=False, low_memory=False)

    # Open output file
    fo = open(FileOut, 'w')

    bedfolder = '/datos/ngs/dato-activo/References/Exomes/bed_files/'
    intervalfolder = '/datos/ngs/dato-activo/References/Exomes/interval_list/'
    # Dictionary for bed files
    bed = {'TWIST': bedfolder + 'Twist_Exome_RefSeq_targets_hg38.bed', 'IDT': bedfolder + 'xgen-exome-research-panel-v2-targets-hg38.bed', 'Agilent': bedfolder + 'S07604514_Regions_agilent_hg38.bed'}
    interval = {'TWIST': intervalfolder + 'TWIST_hg38.interval_list', 'IDT': intervalfolder + 'IDT_hg38.interval_list', 'Agilent': intervalfolder + 'SureSelect_v6_hg38.interval_list'}
    # Write header
    #fo.write("%s\n" %('\t'.join(l_th[1:])))

    # Compute percentages
    ss['bed'] = ss['protocol'].map(bed)
    ss['interval'] = ss['protocol'].map(interval)

    ss.to_csv(fo, index = False)
    fi.close()


def main(args=None):
    args = parse_args(args)
    add_bed_file(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
