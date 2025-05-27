import sys
import os
import re
import gzip
import HTSeq as htseq
from multiprocessing import Pool
from contextlib import nullcontext

try:
    from rapidfuzz.distance import Hamming
except: # this is checked in cli.py
    pass

def slice_SequenceWithQualities(swq, start=None):
    # override SequenceWithQualities.__getitem__ appending [part] to swq.name when slicing
    assert isinstance(swq, htseq.SequenceWithQualities)

    if start:
        return htseq.SequenceWithQualities(swq.seq[start:], swq.name, swq.qualstr[start:])
    else:
        return swq

def get_fastq_str(swq):
    # replace SequenceWithQualities.get_fastq_str to not depend on cStringIO
    # otherwise code is adapted from SequenceWithQualities.write_to_fastq_file

    newline = '\n'.encode()

    fastq_str = b''
    fastq_str += ("@%s\n" % swq.name).encode()
    fastq_str += (swq.seq + newline)
    fastq_str += ('+'.encode() + newline)
    fastq_str += (swq.qualstr + newline)
    return fastq_str

def process_entry(entry, pattern, umi_len, only_umi,
                  search_region, anchor_seq, trailing_seq, 
                  fuzzy_umi_params):

    seq_str = entry.seq.decode("utf-8")
    pattern_len = len(anchor_seq) + umi_len + len(trailing_seq)

    if search_region > 0:
        if search_region + pattern_len > len(seq_str): # read to short to search
            return (None, None)
        else:
            seq_str = seq_str[:search_region + pattern_len] # trim to search_region

    match = pattern.search(seq_str) # exact matching using re
    if not match:
        if fuzzy_umi_params:
            match = fuzzy_umi_params['fuzzy_pattern'].search(seq_str) # fuzzy matching

            if match:
                _, m_end = match.span()
                trailing_hit = seq_str[m_end:m_end+len(trailing_seq)]

                # count mismatches, insertions, deletions
                cmis, cins, cdel = match.fuzzy_counts

                # check against mismatch, indel, and trailing dist thresholds
                if cmis <= fuzzy_umi_params['anchor_max_mismatch'] \
                    and (cins + cdel) <= fuzzy_umi_params['anchor_max_indel'] \
                    and Hamming.distance(trailing_seq, 
                                         trailing_hit) <= fuzzy_umi_params['trailing_dist_thresh']:
                    umi = match.group(2)
                else:
                    umi = None
            else:
                umi = None
        else:
            umi = None
    else:
        umi = match.group(0)[len(anchor_seq):len(anchor_seq) + umi_len]

    if umi:
        entry_processed = slice_SequenceWithQualities(entry, match.span()[1])
        entry_processed.name += '_' + umi
        return entry_processed, umi
    else:
        return (None, None) if only_umi else (entry, None)

def none_gen():
    while True:
        yield None

def precompile_regex(umi_len, anchor_seq, trailing_seq, fuzzy_umi_params):
    if fuzzy_umi_params:
        import regex
        anchor_max_mismatch = fuzzy_umi_params['anchor_max_mismatch']
        anchor_max_indel = fuzzy_umi_params['anchor_max_indel']

        anchor_fuzzy = rf"({anchor_seq}){{e<={anchor_max_mismatch + anchor_max_indel}}}"
        umi_capture = rf"([ACGTN]{{{umi_len}}})"

        # fuzzy doesnt use trailing seq, checked after match
        fuzzy_pattern = regex.compile(anchor_fuzzy + umi_capture, flags=regex.BESTMATCH)
        fuzzy_umi_params['fuzzy_pattern'] = fuzzy_pattern
    
    pattern = re.compile(f"({anchor_seq})[NGCAT]{{{umi_len}}}({trailing_seq})") # exact match
    return pattern

def process_fastq(paths, outnames, umi_len, 
                  only_umi=False, 
                  anchor_seq='ATTGCGCAATG', 
                  trailing_seq='GGG',
                  search_region=-1,
                  min_remaining_seqlen=-1,
                  fuzzy_umi_params=None):

    r1_path, r2_path = paths
    r1_out_path, r2_out_path = outnames

    readcount = 0
    umicount = 0
    reads_written = 0
    readnames = set()

    pattern = precompile_regex(umi_len, anchor_seq, trailing_seq, fuzzy_umi_params)

    # prepare output files
    r1_out = gzip.open(r1_out_path, 'wb', compresslevel=6) # matches bash gzip
    r2_out = gzip.open(r2_out_path, 'wb', compresslevel=6) if r2_path and r2_out_path else None

    # setup fastq readers: if not using R2 we set it to a null generator
    r1_reader = htseq.FastqReader(r1_path)
    r2_reader = htseq.FastqReader(r2_path) if r2_path else nullcontext(none_gen())

    # traverse both files in tandem
    with r1_reader as r1_in, r2_reader as r2_in:

        for entry1, entry2 in zip(r1_in, r2_in):
            readcount += 1
            entry1.name, _, _ = entry1.name.partition(' ')

            # check for duplicated R1 readname, R2 is checked for consistency below
            if entry1.name in readnames:
                sys.exit('duplicate R1 name found at read number %s' %readcount)
            else:
                readnames.add(entry1.name)

            # process R1 to search for UMI and trim
            entry1_processed, umi = process_entry(entry1, pattern, umi_len, only_umi,
                                                  search_region, anchor_seq, trailing_seq,
                                                  fuzzy_umi_params)
            
            # failure condition when sequence short or UMI not detected
            if entry1_processed is None or len(entry1_processed.seq) <= max(min_remaining_seqlen, 0):
                continue

            if umi:
                umicount += 1

            # handle R2 readname: check it matches R1, also add umi
            if entry2 is not None:
                entry2.name, _, _ = entry2.name.partition(' ')
                if entry1.name != entry2.name:
                    sys.exit('readname mismatch in R1 and R2 at read number %s' %readcount)
                if umi:
                    entry2.name += '_' + umi

            # write output
            if only_umi and umi is None:
                continue
            else:
                r1_out.write(get_fastq_str(entry1_processed))

                if r2_out and entry2 is not None:
                    r2_out.write(get_fastq_str(entry2))

                reads_written += 1

    if readcount > 0:
        print( (f"{os.path.basename(paths[0])}: {readcount} reads, " 
                f"{umicount} with UMI ({(umicount/readcount)*100:.2f}%), "
                f"{reads_written} written "
                f"({((readcount-reads_written)/readcount)*100:.2f}% skipped)") )
    else:
        print(f"empty input file")
    
    r1_out.close()
    if r2_out:
        r2_out.close()

def _fastq_worker(task_args):
    # worker function for multiprocessing
    infile, outfile, kwargs = task_args

    try:
        process_fastq(infile, outfile, **kwargs)
    except Exception as e:
        new_msg = f"in {infile}:\n{e}"
        raise type(e)(new_msg) from e
    
def process_fastq_parallel(filepairs, umi_len, num_workers=4,
                           only_umi=False, 
                           anchor_seq='ATTGCGCAATG', 
                           trailing_seq='GGG',
                           search_region=-1,
                           min_remaining_seqlen=-1,
                           fuzzy_umi_params=None):

    # wrapper around process_fastq to handle multiple FASTQ files in parallel
    # filepairs are a list of tasks as tuples of ( (r1_path, r2_path), (r1_out, r2_out) )

    # bundle constant args
    kwargs = dict(
        umi_len=umi_len,
        only_umi=only_umi,
        anchor_seq=anchor_seq,
        trailing_seq=trailing_seq,
        search_region=search_region,
        min_remaining_seqlen=min_remaining_seqlen,
        fuzzy_umi_params=fuzzy_umi_params
    )

    # map the worker over the filepairs
    tasks = [(infile, outfile, kwargs) for (infile, outfile) in filepairs]
    with Pool(num_workers) as pool:
        results = pool.map(_fastq_worker, tasks)
    