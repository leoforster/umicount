"""CLI entrypoints for umicount"""

import argparse
import logging
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

def configure_logging(logfile_path):

    root = logging.getLogger()
    root.setLevel(logging.INFO)

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    # two handlers on the root logger: first for info
    if logfile_path in (sys.stdout, None, "-", ""):
        fh = logging.StreamHandler(sys.stdout)
    else:
        fh = logging.FileHandler(logfile_path, mode='w')
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    root.addHandler(fh)

    # second handler for errors
    eh = logging.StreamHandler(sys.stderr)
    eh.setLevel(logging.WARNING)
    eh.setFormatter(fmt)
    root.addHandler(eh)

    return root

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
    parser.add_argument("-c", "--cores", action="store", type=int,
                        help="Number of cores for processing R1/R2: each core processes 1 sample")
    parser.add_argument("-l", "--logfile", action="store", default=None,
                        help="Path to log file, defaults to standard output (sys.stdout)")
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
    parser.add_argument("-f", "--fuzzy_umi", action="store_true", default=False, 
                        help=( "Enable mismatches/indels in pre/pos-UMI sequence detection, " 
                               "requires python regex library to be installed") )
    parser.add_argument("--min_seqlen", action="store", type=int, default=-1,
                        help="Threshold for minimum remaining sequence length after trimming UMI")
    parser.add_argument("--anchor_mismatches", action="store", type=int, default=2, 
                        help="Max mismatches allowed in anchor in fuzzy UMI extraction")
    parser.add_argument("--anchor_indels", action="store", type=int, default=1, 
                        help="Max indels allowed in anchor in fuzzy UMI extraction")
    parser.add_argument("--trailing_hamming_threshold", action="store", type=int, default=2, 
                        help="Maximum Hamming distance in trailing sequence")

    r = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # set up logging
    logger = configure_logging(r.logfile)

    if r.read2:
        if len(r.read1) != len(r.read2):
            logger.error('umiextract requires the number of R2 files match R1 if giving R2 (-2, --read2)')
            sys.exit()

    # check implicity fuzzy_umi usage
    if r.fuzzy_umi == False:
        fuzzy_args = ['--anchor_mismatches', '--anchor_indels', '--trailing_hamming_threshold']
        fuzzy_used = [i for i in fuzzy_args if i in sys.argv]

        if len(fuzzy_used) > 0:
            logger.info(f"Note: supplied {', '.join(fuzzy_used)} without --fuzzy_umi, setting --fuzzy_umi=True")
            r.fuzzy_umi = True

    # set up fuzzy matching
    fuzzy_umi_params={'anchor_max_mismatch':r.anchor_mismatches,
                      'anchor_max_indel':r.anchor_indels,
                      'trailing_dist_thresh':r.trailing_hamming_threshold} if r.fuzzy_umi else None

    if fuzzy_umi_params:
        try:
            import regex
            from rapidfuzz.distance import Hamming
        except ModuleNotFoundError:
            logger.error('umiextract requires both regex and rapidfuzz packages for fuzzy UMI detection')
            sys.exit()

    # collate R1-R2 pairs for threading
    filedir = os.path.abspath(os.path.expanduser(r.output_dir))
    readpairs = [(r.read1[n], r.read2[n] if r.read2 else None) for n in range(len(r.read1))]
    readpairs = sorted(readpairs, key=lambda i: os.path.getsize(i[0]), reverse=True) # large files first

    def prepend_output(filename, filedir):
        filename = os.path.basename(filename)
        fbase, fext = os.path.splitext(filename)
        if fext == ".gz": 
            fbase, _ = os.path.splitext(fbase)
        fout = f"{fbase}_umiextract.fastq.gz"
        return os.path.join(filedir, fout)

    filepairs = [] # store input + output files
    for read1, read2 in readpairs:

        # prepend umiextract before extension on output filenames
        r1_out = prepend_output(read1, filedir)
        r2_out = prepend_output(read2, filedir) if read2 else None
        
        filepairs.append( ((read1, read2), (r1_out, r2_out)) )

    # run umiextract over parallel process_fastq wrapper
    if len(filepairs) > 0:
        process_fastq_parallel(filepairs, r.umilen, num_workers=r.cores,
                               only_umi=r.only_umi, 
                               anchor_seq=r.anchor, 
                               trailing_seq=r.trailing,
                               search_region=r.search_region,
                               min_remaining_seqlen=r.min_seqlen,
                               fuzzy_umi_params=fuzzy_umi_params)

def umicount():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--bams", type=existing_file, nargs="+",
                        help="input bamfiles from modified fastqs, sorted by read")
    parser.add_argument("-d", "--output_dir", type=existing_dir, default=os.getcwd(), 
                        help="Directory to output counts matrices, named as umicounts.*.tsv")
    parser.add_argument("-c", "--cores", action="store", type=int,
                        help="Number of cores for processing BAMs: each core processes 1 sample")
    parser.add_argument("-l", "--logfile", action="store", default=None,
                        help="Path to log file, defaults to standard output (sys.stdout)")
    parser.add_argument("--tmp_dir", action="store", type=existing_dir, default=None,
                        help="directory to save temporary files")
    parser.add_argument("-g", "--gtf", type=existing_file, help="input GTF file (ensembl format)")
    parser.add_argument("--GTF_dump", type=existing_dir, default=None, help="File path to dump parsed GTF data")
    parser.add_argument("--GTF_skip_parse", type=existing_file, default=None, help="Path to dumped GTF data")
    parser.add_argument("--no_dedup", action="store_true", default=False, help="dont deduplicate UMI counts")
    parser.add_argument("--mm_count_primary", action="store_true", default=False,
                        help=( "count primary alignment (BAM flag 0x100) for multimapping reads.") )
    parser.add_argument("--min_read_mapQ", action="store", type=int, default=0,
                        help="Minimum mapQ to keep read. Checks for aligned status but not whether mapQ is 255.")
    parser.add_argument("--multiple_primary_action", action="store", default="warn",
                        choices=['warn', 'raise', 'skip'],
                        help=( "how to handle cases when a read has multiple primary alignments: "
                               "warn (default) - warn user and pick random alignment to count, "
                               "raise - raise error and exit, skip - alignment is counted as multimapping") )
    parser.add_argument("--combine_unspliced", action="store_true", default=False,
                        help="dont distinguish spliced and unspliced UMI counts, instead as U")
    parser.add_argument("-u", "--UMI_correct", action="store_true", default=False, 
                        help=( "Enable gene-wise UMI correction by hamming distance " 
                               "requires python RapidFuzz library to be installed") )
    parser.add_argument("--hamming_threshold", action="store", type=int, default=1,
                        help="Hamming distance threshold for merging similar UMIs, usually 1")
    parser.add_argument("--count_ratio_threshold", action="store", type=int, default=2,
                        help=("Threshold where UMIs are only merged if one has "
                              "at least (threshold*counts)-1 as many counts as the other, usually 2"))
    
    r = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # set up logging
    logger = configure_logging(r.logfile)

    if r.gtf is None and r.GTF_skip_parse is None:
        logger.error('umicount requires one of --gtf, --GTF_skip_parse')
        sys.exit()

    if (r.bams is None) and not (r.gtf and r.GTF_dump):
        logger.error('umicount found empty input (--bams), skipping input is only valid with --gtf and --GTF_dump')
        sys.exit()

    # case when parsing GTF and pickling contents
    if r.gtf and r.GTF_dump:
        load_gtf_data(r.gtf, skipgtf=None, dumpgtf=r.GTF_dump)
        sys.exit()

    if r.tmp_dir is None:
        r.tmp_dir = r.output_dir

    # populate read types to count (as columns)
    basecols = []
    if r.combine_unspliced:
        basecols += ['U', 'R']
    else:
        basecols += ['UE', 'UI', 'RE', 'RI']
    if r.no_dedup == False:
        basecols += ['D']

    # check implicity fuzzy_umi usage
    if r.UMI_correct == False:
        umicorr_args = ['--count_ratio_threshold', '--hamming_threshold']
        umicorr_used = [i for i in umicorr_args if i in sys.argv]

        if len(umicorr_used) > 0:
            print(f"Note: supplied {', '.join(umicorr_used)} without --UMI_correct, setting --UMI_correct=True")
            r.UMI_correct = True

    # set up UMI correction
    umi_correct_params = {'countratio_threshold':r.count_ratio_threshold,
                          'hamming_threshold':r.hamming_threshold} \
                          if r.UMI_correct else None

    if umi_correct_params:
        try:
            from rapidfuzz.distance import Hamming
        except ModuleNotFoundError:
            logger.error('umicount requires the rapidfuzz package for UMI correction')
            sys.exit()

    # collate input and output files in filepairs
    filedir = os.path.abspath(os.path.expanduser(r.output_dir))
    bamfiles = sorted(r.bams, key=lambda i: os.path.getsize(i), reverse=True)

    # run umicount on BAMs split by thread
    if len(bamfiles) > 0:
        gtf_data = load_gtf_data(r.gtf, skipgtf=r.GTF_skip_parse) # only need to load once
        process_bam_parallel(bamfiles, filedir, gtf_data, num_workers=r.cores,
                             cols_to_use=basecols, 
                             tmp_dir=r.tmp_dir,
                             min_read_mapQ=r.min_read_mapQ,
                             count_primary=r.mm_count_primary,
                             multiple_primary_action=r.multiple_primary_action,
                             umi_correct_params=umi_correct_params)
