"""CLI entrypoints for umicount"""

import argparse
import sys
import os

from .umiextract import process_fastq_parallel
from .umicount import process_bam_parallel, load_gtf_data

def existing_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File '{path}' does not exist.")
    return path

def existing_dir(path):
    if not os.path.exists(os.path.dirname(path)):
        raise argparse.ArgumentTypeError(f"Folder '{path}' does not exist.")
    return path

def umiextract():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-1", "--read1", action="store", nargs="+", 
                        required=True, type=existing_file, 
                        help="Input smartseq3 read1 fastq")
    parser.add_argument("-2", "--read2", action="store", nargs="*", 
                        required=False, default=None, type=existing_file,
                        help="Input smartseq3 read2 fastq")
    parser.add_argument("-d", "--output_dir", action="store", type=existing_dir, default=os.getcwd(),
                        help="Directory to output processed fastq files")
    parser.add_argument("-c", "--cores" action="store", type=int,
                        help="Number of cores for processing R1/R2: each core processes 1 sample")
    parser.add_argument("--umilen", action="store", type=int, default=8,
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

    if r.read2:
        if len(r.read1) != len(r.read2):
            sys.exit("if giving R2, number of R2 files must match R1")

    # set up fuzzy matching
    fuzzy_umi_params={'anchor_max_mismatch':r.anchor_mismatches,
                      'anchor_max_indel':r.anchor_indels,
                      'min_trailing_G':r.trailing_Gs} if r.fuzzy_umi else None

    if fuzzy_umi_params:
        try:
            import regex
        except ModuleNotFoundError:
            sys.exit('requires regex package to enable fuzzy_umi_extraction')

    # collate R1-R2 pairs for threading
    filepairs = []
    filedir = os.path.abspath(os.path.expanduser(r.output_dir))

    def prepend_output(filename, filedir):
        filename = os.path.basename(filename)
        fbase, fext = os.path.splitext(filename)
        if fext == ".gz": 
            fbase, _ = os.path.splitext(fbase)
        fout = f"{fbase}_umiextract.fastq.gz"
        return os.path.join(filedir, fout)

    for n, read1 in enumerate(r.read1):
        read2 = r.read2[n] if r.read2 else None

        # prepend umiextract before extension on output filenames
        r1_out = prepend_output(read1, filedir)
        r2_out = prepend_output(read2, filedir) if read2 else None
        
        filepairs.append( ((r1_file, r2_file), (r1_out, r2_out)) )

    # run umiextract over parallel process_fastq wrapper
    if len(filepairs) > 0:
        process_fastq_parallel(filepairs, r.umilen, num_workers=r.cores,
                               only_umi=r.only_umi, 
                               anchor_seq=r.anchor, 
                               trailing_seq=r.trailing,
                               search_region=r.search_region,
                               min_remaining_seqlen=r.min_seqlen,
                               fuzzy_umi_params=fuzzy_umi_params):

def umicount():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--bams", type=existing_file, nargs="+",
                        help="input bamfiles from modified fastqs, sorted by read")
    parser.add_argument("-d", "--output_dir", type=existing_dir, default=os.getcwd(), 
                        help="Path to output counts matrix")
    parser.add_argument("-c", "--cores" action="store", type=int,
                        help="Number of cores for processing BAMs: each core processes 1 sample")
    parser.add_argument("-g", "--gtf", type=existing_file, help="input GTF file (ensembl format)")
    parser.add_argument("--GTF_dump", type=existing_dir, default=None, help="File path to dump parsed GTF data")
    parser.add_argument("--GTF_skip_parse", type=existing_file, default=None, help="Path to dumped GTF data")
    parser.add_argument("--no_dedup", action="store_true", default=False,
                        help="dont deduplicate UMI counts")
    parser.add_argument("--combine_unspliced", action="store_true", default=False,
                        help="dont distinguish spliced and unspliced UMI counts, instead as U")
    parser.add_argument("-u", "--UMI_correct", action="store_true", default=False, 
                        help=( "Enable gene-wise UMI correction by hamming distance " 
                               "requires python RapidFuzz library to be installed") )
    parser.add_argument("--hamming_threshold", action="store", type=int, default=1,
                        help="Hamming distance threshold for merging similar UMIs, usually 1")
    parser.add_argument("--count_ratio_threshold", action="store", type=int, default=2,
                        help=("Threshold where UMIs are only merged if one has more "
                              "counts by a factor of (threshold*counts)-1"))
    
    r = parser.parse_args()

    if r.gtf is None and r.GTF_skip_parse is None:
        sys.exit('require one of --gtf, --GTF_skip_parse')

    if (len(r.bams) == 0 or r.output_dir is None) and not (r.gtf and r.GTF_dump):
        sys.exit('no input files found, skipping input only valid with --gtf and --GTF_dump')

    # case when parsing GTF and pickling contents
    if r.gtf and r.GTF_dump:
        load_gtf_data(r.gtf, skipgtf=None, dumpgtf=r.GTF_dump)

    # populate read types to count (as columns)
    basecols = []
    if r.combine_unspliced:
        basecols += ['U', 'R']
    else:
        basecols += ['UE', 'UI', 'RE', 'RI']
    if r.no_dedup == False:
        basecols += ['D']

    # set up UMI correction
    umi_correct_params = {'countratio_threshold':r.count_ratio_threshold,
                          'hamming_threshold':r.hamming_threshold} \
                          if r.UMI_correct else None

    if umi_correct_params:
        try:
            from rapidfuzz.distance import Hamming
        except ModuleNotFoundError:
            sys.exit('requires rapidfuzz package to enable UMI correction')

    # collate input and output files in filepairs
    filepairs = []
    filedir = os.path.abspath(os.path.expanduser(r.output_dir))
    for n, bamfile in enumerate(r.bams):
        base, ext = os.path.splitext(bamfile)
        filepairs.append( (bamfile, os.path.join(filedir, f"{base}.umicounts")) )

    # run umicount on BAMs split by thread
    if len(filepairs) > 0:
        gtf_data = load_gtf_data(gtffile, skipgtf=r.GTF_skip_parse) # only need to load once
        process_bam_parallel(filepairs, gtf_data, num_workers=r.cores,
                             cols_to_use=basecols, 
                             umi_correct_params=umi_correct_params)
