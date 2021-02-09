#!/usr/bin/env python


import pandas as pd
import argparse
import sys

def parse_args(args=None):
    Description = 'Generate table for coverage and %on target'
    Epilog = """Example usage: python ontarget_coverage.py <FILE_IN> <FILE_OUT> -s <sample>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input bed file from mosdepth.")
    parser.add_argument('FILE_OUT', help="Output file.")
    parser.add_argument('-s', '--sample', type=str, dest="SAMPLE", default='', help="sample name")
    return parser.parse_args(args)


def generate_table(FileIn,FileOut,sample):
    #Open input file
    fi = open(FileIn, 'r')

    # Load mosdepth thresholds.bed.gz into a pandas dataframe
    cov = pd.read_csv(fi, delimiter='\t', header=0, index_col=False, low_memory=False)

    # Open output file
    fo = open(FileOut, 'w')

    coverage_threshold = "0,1,5,10,25,50,100,200,300,400,500"

    # Get thresholds
    l_th = list(map(lambda x: str(x) + 'X', coverage_threshold.split(',')))

    # Write header
    #fo.write("%s\n" %('\t'.join(l_th[1:])))

    # Compute percentages
    l_per = []
    for threshold in l_th[1:]:
        l_per.append((100.0*cov[threshold].sum()) / cov['0X'].sum())

            # Write them to the outputfile
        l_per_str = list(map(lambda x: str(round(x,2)) + '%', l_per))
    fo.write("%s\t%s\n" %(sample, '\t'.join(l_per_str)))
    fo.close()
    fi.close()


def main(args=None):
    args = parse_args(args)
    generate_table(args.FILE_IN,args.FILE_OUT,args.SAMPLE)


if __name__ == '__main__':
    sys.exit(main())
