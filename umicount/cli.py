"""CLI entrypoints for umicount"""

import argparse
import sys
import os

from .umiextract import process_fastx
from .umicount import process_bam


def process_smartseq():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-1", "--read1", action="store", required=True,
                        help="input smartseq3 read1 fastq")
    parser.add_argument("-2", "--read2", action="store", required=False,
                        help="input smartseq3 read2 fastq")
    parser.add_argument("--r1_out", action="store", help="read1 output file")
    parser.add_argument("--r2_out", action="store", help="read2 output file")
    parser.add_argument("--umilen", action="store", type=int, default=12,
                        help="umi length in bp")
    parser.add_argument("-u", "--only_umi", action="store_true", default=False,
                        help="only output UMI-containing PE reads")
    results = parser.parse_args()

    process_fastx((results.read1, results.read2),
                (results.r1_out, results.r2_out),
                results.umilen, results.only_umi)

def umicount():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--files", action="store", default=[], nargs="*",
                        help="input bamfiles from modified fastqs, sorted by read")
    parser.add_argument("-c", "--columns", action="store", default=None,
                        nargs="*", help="output column basenames")
    parser.add_argument("--dumpgtf", action="store",
                        default="umicounts_gtf.dump",
                        help="dump GTF parsing output to given file")
    parser.add_argument("--skipgtf", action="store",
                        help="skip GTF parsing and read input from file")
    parser.add_argument("-g", "--gtf", action="store",
                        help="input GTF file (ensembl format)")
    parser.add_argument("-d", "--nodupes", action="store_true", default=False,
                        help="dont report UMI duplicates per gene per cell")
    parser.add_argument("-o", "--output", action="store", required=True,
                        help="output counts matrix")
    parser.add_argument("--outgeneformat", action="store", default="geneid",
                        help="output counts matrix with 'geneid' or 'genename'")
    parser.add_argument("--cores", action="store", type=int, default=22,
                        help="number of cores for multiprocessing")
    results = parser.parse_args()

    if results.columns:
        assert len(results.columns) == len(results.files), "column names and file names have different lengths"
        colnames = results.columns
    else:
        colnames = list(map(os.path.basename, results.files))

    if results.gtf:
        if results.skipgtf:
            print('--gtf and --skipgtf given, will read from skipgtf')
        elif not os.path.exists(results.gtf):
            print('invalid GTF file path, exiting')
            sys.exit()

    if results.skipgtf:
        if not os.path.exists(results.skipgtf):
            print('invalid --skipgtf file path, exiting')
            sys.exit()

    if [results.gtf, results.skipgtf] == [None, None]:
        print('require one of --gtf, --skipgtf')
        sys.exit()

    if results.cores < 1:
        print('invalid number of cores, exiting')
        sys.exit()

    if not results.outgeneformat in ['geneid', 'genename']:
        print("invalid --outgeneformat, use one of 'geneid' or 'genename'")
        sys.exit()

    if len(results.files) == 0 and not results.gtf:
        print('no input files found, skipping input only valid with --gtf')
        sys.exit()
    for i in results.files:
        if not os.path.exists(i):
            print('file %s not found, exiting' %i)
            sys.exit()

    process_bam(results.files, results.gtf, results.output, results.cores,
                results.dumpgtf, results.skipgtf, results.nodupes, colnames,
                results.outgeneformat)
