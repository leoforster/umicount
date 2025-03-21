import sys
import os
import re
import gzip
import HTSeq as htseq

def slice_SequenceWithQualities(swq, start=None):
    # override SequenceWithQualities.__getitem__ behavior that appends [part] to swq.name when slicing
    assert isinstance(swq, htseq.SequenceWithQualities)

    if start:
        return htseq.SequenceWithQualities(swq.seq[start:], swq.name, swq.qualstr[start:])
    else:
        return swq

def get_fastq_str(swq):
    # replace SequenceWithQualities.get_fastq_str to not depend on cStringIO
    # otherwise code is adapted from SequenceWithQualities.write_to_fastq_file, which actually creates the fastq string...

    newline = '\n'.encode()

    fastq_str = b''
    fastq_str += ("@%s\n" % swq.name).encode()
    fastq_str += (swq.seq + newline)
    fastq_str += (newline)
    fastq_str += (swq.qualstr + newline)
    return fastq_str

def process_fastq(paths, outnames, umilen, only_umi, pre_umi_seq='ATTGCGCAATG', post_umi_seq='GGG'):
    r1_path, r2_path = paths
    ss3_regex = "(%s)[NGCAT]{%s}(%s)" % (pre_umi_seq, umilen, post_umi_seq)
    umicount = 0
    read_umis = {}
    empty_sequences = {}

    r1_out = gzip.open(outnames[0], 'wb')
    r2_out = gzip.open(outnames[1], 'wb') if outnames[1] else None

    r1count = 0
    with htseq.FastqReader(r1_path) as r1_in:
        for n, entry in enumerate(r1_in):
            entry_seq_str = entry.seq.decode("utf-8")

            try:
                x = read_umis[entry.name]
                print('non-unique read name: %s' % entry.name)
                print('exiting')
                sys.exit()
            except KeyError:
                read_umis[entry.name] = None

            m = re.search(ss3_regex, entry_seq_str)
            if m:
                umicount += 1
                match = m.group(0)
                mpos = m.span()[1]
                umi = match[len(pre_umi_seq):(len(pre_umi_seq) + umilen)]
                read_umis[entry.name] = umi

                if mpos == len(entry_seq_str):
                    empty_sequences[entry.name] = 1
                    continue

                entry.name = entry.name + "_" + umi
                entry = slice_SequenceWithQualities(entry, mpos) # avoid HTSeq appending '[part]'

            if only_umi and not m:
                continue
            else:
                r1_out.write(get_fastq_str(entry)) # returns a byte string
                r1count += 1

    if r2_path:
        r2count = 0
        with htseq.FastqReader(r2_path) as r2_in:
            for n, entry in enumerate(r2_in):
                if entry.name in empty_sequences:
                    continue

                if entry.name in read_umis and read_umis[entry.name]:
                    umi = read_umis[entry.name]
                    entry.name = entry.name + "_" + umi
                else:
                    if only_umi:
                        continue

                r2_out.write(get_fastq_str(entry))
                r2count += 1

    print(umicount, "umis,", r1count, "reads written,", len(empty_sequences.keys()), "reads removed")
    r1_out.close()
    if r2_out:
        r2_out.close()
