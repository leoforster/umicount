import sys
import os
import re
import gzip
import HTSeq as htseq
from contextlib import nullcontext

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
                  search_region, min_remaining_seqlen,
                  anchor_seq, trailing_seq, 
                  fuzzy_umi_params):

    seq_str = entry.seq.decode("utf-8")
    pattern_len = len(anchor_seq) + umi_len + len(trailing_seq)

    if search_region > 0:
        if search_region + pattern_len > len(seq_str):
            print( ('read length not compatible with search_region '
                    'of %s and %sbp pattern' %(search_region, pattern_len)) )
        else:
            seq_str = seq_str[:search_region + pattern_len] # trim to search_region

    match = pattern.search(seq_str)
    if not match:
        return (None, None) if only_umi else (entry, None)

    m_start, m_end = match.span()
    if len(seq_str) - (m_start + pattern_len) <= max(min_remaining_seqlen, 0):
        return (None, None) # skip reads where trimming leaves bp < threshold

    if fuzzy_umi_params:
        trailing_seq_hit = seq_str[m_end:m_end+len(trailing_seq)]

        # count mismatches, insertions, deletions
        cmis, cins, cdel = match.fuzzy_counts

        # check against mismatch, indel, and G-count thresholds
        if cmis <= anchor_max_mismatch \
            and (cins + cdel) <= anchor_max_indel \
            and trailing_seq_hit.count('G') >= min_trailing_G:
            umi = match.group(2)
    else:
        umi = match.group(0)[len(anchor_seq):len(anchor_seq) + umi_len]

    entry_processed = slice_SequenceWithQualities(entry, m_end)
    entry_processed.name += '_' + umi
    return entry_processed, umi

def none_gen():
    while True:
        yield None

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

    if fuzzy_umi_params:
        try:
            import regex
        except ModuleNotFoundError:
            print('requires regex package to enable fuzzy_umi_extraction')
            sys.exit()

        anchor_max_mismatch = fuzzy_umi_params['anchor_max_mismatch']
        anchor_max_indel = fuzzy_umi_params['anchor_max_indel']
        min_trailing_G = fuzzy_umi_params['min_trailing_G']
        anchor_edits = anchor_max_mismatch + anchor_max_indel

    # precompile regex
    if fuzzy_umi_params:
        anchor_fuzzy = rf"({anchor_seq}){{e<={anchor_edits}}}"
        umi_capture = rf"([ACGTN]{{{umi_len}}})"
        pattern_str = anchor_fuzzy + umi_capture
        pattern = regex.compile(pattern_str, flags=regex.BESTMATCH)
    else:
        pattern = re.compile(f"({anchor_seq})[NGCAT]{{{umi_len}}}({trailing_seq})")

    # prepare output files
    r1_out = gzip.open(r1_out_path, 'wb')
    r2_out = gzip.open(r2_out_path, 'wb') if r2_path and r2_out_path else None

    # setup fastq readers: if not using R2 we set it to a null generator
    r1_reader = htseq.FastqReader(r1_path)
    r2_reader = htseq.FastqReader(r2_path) if r2_path else nullcontext(none_gen())

    # traverse both files in tandem
    with r1_reader as r1_in, r2_reader as r2_in:

        for entry1, entry2 in zip(r1_in, r2_in):
            readcount += 1

            # check for duplicated R1 readname, R2 is checked for consistency below
            if entry1.name in readnames:
                print('duplicate R1 name found at read number %s' %readcount)
                sys.exit()
            else:
                readnames.add(entry1.name)

            # process R1 to search for UMI and trim
            entry1_processed, umi = process_entry(entry1, pattern, umi_len, only_umi,
                                                  search_region, min_remaining_seqlen,
                                                  anchor_seq, trailing_seq,
                                                  fuzzy_umi_params)
            
            if entry1_processed is None: # case when entire read is UMI + flanking
                continue

            if umi:
                umicount += 1

            # handle R2 readname: check it matches R1, also add umi
            if entry2 is not None:
                if entry1.name != entry2.name:
                    print('readname mismatch in R1 and R2 at read number %s' %readcount)
                    sys.exit()
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
        print( (f"{readcount} reads: {umicount} with UMI ({(umicount/readcount)*100:.2f}%), "
                f"{reads_written} written "
                f"({((readcount-reads_written)/readcount)*100:.2f}% skipped)") )
    else:
        print(f"empty input file")
    
    r1_out.close()
    if r2_out:
        r2_out.close()
