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

def process_entry(entry, umi_pattern, umi_len, only_umi, pre_umi_seq):
    # extract and trim UMI if present

    seq_str = entry.seq.decode("utf-8")
    match = umi_pattern.search(seq_str)
    if not match:
        return (None, None) if only_umi else (entry, None)

    m_start, m_end = match.span()
    if m_end == len(seq_str):
        return None, None  # entire read is just UMI + flanking, skip

    umi = match.group(0)[len(pre_umi_seq):len(pre_umi_seq) + umi_len]

    entry_processed = slice_SequenceWithQualities(entry, m_end)
    entry_processed.name += '_' + umi
    return entry_processed, umi

def none_gen():
    while True:
        yield None

def process_fastq(paths, outnames, umi_len, only_umi, 
                  pre_umi_seq='ATTGCGCAATG', 
                  post_umi_seq='GGG'):

    r1_path, r2_path = paths
    r1_out_path, r2_out_path = outnames

    readcount = 0
    umicount = 0
    reads_written = 0
    readnames = set()

    # precompile regex
    umi_pattern = re.compile(f"({pre_umi_seq})[NGCAT]{{{umi_len}}}({post_umi_seq})")

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
            entry1_processed, umi = process_entry(entry1, umi_pattern, 
                                                  umi_len, only_umi, pre_umi_seq)
            
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
        print(f"{readcount} reads: {umicount} with UMI ({(umicount/readcount)*100:.2f}%), \
              {reads_written} written ({(reads_written/readcount)*100:.2f}% skipped)")
    else:
        print(f"empty input file")
    
    r1_out.close()
    if r2_out:
        r2_out.close()
