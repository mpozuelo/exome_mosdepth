#!/usr/bin/env python


import pandas as pd
import argparse
import sys

def parse_args(args=None):
    Description = 'Check samplesheet and add files to bed file'
    Epilog = """Example usage: python modify_samplesheet.py <FILE_IN> <FILE_OUT> -r <total_reads.tsv>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet.")
    parser.add_argument('FILE_OUT', help="Output samplesheet.")
    parser.add_argument('-c', '--counts', type=str, dest="COUNTS", default='', help="Total counts per sample in bam")
    return parser.parse_args(args)


def add_percentage(FileIn,FileOut,COUNTS):
    #Open input file
    fi = open(FileIn, 'r')
    cnt = open(counts, 'r')

    # Load mosdepth thresholds.bed.gz into a pandas dataframe
    ss = pd.read_csv(fi, delimiter=',', index_col=False, low_memory=False)
    counts = pd.read_csv(cnt, delimiter='\t', index_col=False, low_memory=False)

    # Open output file
    fo = open(FileOut, 'w')

    min = counts['counts'].min()
    counts['percentage'] = round(counts['counts']/min, 2)

    # Dictionary for bed files
    # Write header
    #fo.write("%s\n" %('\t'.join(l_th[1:])))

    # Compute percentages
    ss['percentage'] = ss['sampleID'].map(counts.set_index('sampleID')['percentage'])

    ss.to_csv(fo, index = False)
    fi.close()


def main(args=None):
    args = parse_args(args)
    add_percentage(args.FILE_IN,args.FILE_OUT, args.COUNTS)


if __name__ == '__main__':
    sys.exit(main())
