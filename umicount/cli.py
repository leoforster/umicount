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

    r = parser.parse_args()

    process_fastq((r.read1, r.read2),
                  (r.r1_out, r.r2_out),
                  r.umilen, r.only_umi,
                  anchor_seq=r.anchor, trailing_seq=r.trailing,
                  search_region=r.search_region, min_remaining_seqlen=r.min_seqlen,
                  fuzzy_umi_params={'anchor_max_mismatch':r.anchor_mismatches,
                                    'anchor_max_indel':r.anchor_indels,
                                    'min_trailing_G':r.trailing_Gs} if r.fuzzy_umi else None
                 )

def umicount():

    def existing_file(path):
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(f"File '{path}' does not exist.")
        return path

    def existing_dir(path):
        if not os.path.exists(os.path.dirname(path)):
            raise argparse.ArgumentTypeError(f"Folder '{path}' does not exist.")
        return path

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--files", type=existing_file, default=[], nargs="*",
                        help="input bamfiles from modified fastqs, sorted by read")
    parser.add_argument("-g", "--gtf", type=existing_file, help="input GTF file (ensembl format)")
    parser.add_argument("--GTF-dump", type=existing_file, help="File path to dump parsed GTF data")
    parser.add_argument("--GTF-skip-parse", type=existing_dir, help="Path to dumped GTF data")
    parser.add_argument("-d", "--nodupes", action="store_true", default=False,
                        help="dont report UMI duplicates per gene per cell")
    parser.add_argument("-o", "--output", type=existing_dir, required=True, 
                        help="Path to output counts matrix")
    
    r = parser.parse_args()

    if r.gtf is None and r.GTF_skip_parse is None:
        print('require one of --gtf, --GTF-skip-parse')
        sys.exit()

    if len(r.files) == 0 and not r.gtf:
        print('no input files found, skipping input only valid with --gtf')
        sys.exit()

    process_bam(r.files, r.gtf, r.output,
                r.GTF_dump, r.GTF_skip_parse, r.nodupes)
