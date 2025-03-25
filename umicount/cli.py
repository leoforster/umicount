"""CLI entrypoints for umicount"""

import argparse
import sys
import os

from .umiextract import process_fastq
from .umicount import process_bam

def umiextract():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-1", "--read1", action="store", required=True,
                        help="Input smartseq3 read1 fastq")
    parser.add_argument("-2", "--read2", action="store", required=False,
                        help="Input smartseq3 read2 fastq")
    parser.add_argument("--r1_out", action="store", required=True,
                        help="Output file with R1 UMI-extracted reads")
    parser.add_argument("--r2_out", action="store", 
                        help="Output file with R2 UMI-extracted reads")
    parser.add_argument("--umilen", action="store", type=int, default=12,
                        help="UMI length in bp")
    parser.add_argument("-u", "--only_umi", action="store_true", default=False,
                        help="Only output UMI-containing reads")
    parser.add_argument("-a", "--anchor", action="store", default='ATTGCGCAATG',
                        help="pre-UMI anchor sequence representing TSO")
    parser.add_argument("-t", "--trailing", action="store", default='GGG',
                        help="post-UMI trailing sequence, commonly GGG")
    parser.add_argument("-s", "--search_region", action="store", type=int, default=-1,
                        help="Sequence segment within which to search for UMI, as a bp threshold")
    parser.add_argument("-l", "--min_seqlen", action="store", type=int, default=-1,
                        help="Threshold for minimum remaining sequence length after trimming UMI")
    parser.add_argument("-f", "--fuzzy_umi", action="store_true", default=False, 
                        help=( "Enable mismatches/indels in pre/pos-UMI sequence detection, " 
                               "requires python regex library to be installed") )
    parser.add_argument("--anchor_mismatches", action="store", type=int, default=2, 
                        help="Max mismatches allowed in anchor in fuzzy UMI extraction")
    parser.add_argument("--anchor_indels", action="store", type=int, default=1, 
                        help="Max indels allowed in anchor in fuzzy UMI extraction")
    parser.add_argument("--trailing_Gs", action="store", type=int, default=2, 
                        help="Minimum number of Gs required in trailing sequence")

    results = parser.parse_args()

    process_fastq((results.read1, results.read2),
                  (results.r1_out, results.r2_out),
                  results.umilen, results.only_umi,
                  anchor_seq=results.anchor, trailing_seq=results.trailing,
                  search_region=results.search_region, min_remaining_seqlen=results.min_seqlen,
                  fuzzy_umi_extraction=results.fuzzy_umi, 
                  anchor_max_mismatch=results.anchor_mismatches,
                  anchor_max_indel=results.anchor_indels,
                  min_trailing_G=results.trailing_Gs)

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
        assert len(results.columns) == len(results.files), \
               "column names and file names have different lengths"
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
