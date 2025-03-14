import sys
import os
import re
import gzip

import pysam

def process_fastx(paths, outnames, umilen, only_umi):
    r1_path, r2_path = paths
    ss3_regex = "(ATTGCGCAATG)[NGCAT]{%s}(GGG)" %umilen
    umicount = 0
    read_umis = {}
    empty_sequences = {}

    r1_out = gzip.open(outnames[0], 'wb')
    if outnames[1]:
        r2_out = gzip.open(outnames[1], 'wb')

    r1count = 0
    with pysam.FastxFile(r1_path) as r1_in:
        for n, entry in enumerate(r1_in):

            # have seen readname before?
            try:
                x = read_umis[entry.name]
                print('non-unique read name: %s' %entry.name)
                print('exiting')
                sys.exit()
            except KeyError:
                read_umis[entry.name] = None

            m = re.search(ss3_regex, entry.sequence)
            if m:
                umicount += 1
                match = m.group(0)
                mpos = m.span()[1]
                umi = match[11:(11+umilen)]
                read_umis[entry.name] = umi

                # note for deletion (for R2) and delete empty sequence
                if mpos == len(entry.sequence):
                    empty_sequences[entry.name] = 1
                    continue

                # store UMI in name, trim UMI containing region
                entry.name = entry.name + "_" + umi
                entry.sequence = entry.sequence[mpos:]
                entry.quality = entry.quality[mpos:]

            # encode and compress entry and write
            if only_umi and not m:
                continue
            else:
                r1_out.write((str(entry) + "\n").encode())
                r1count += 1

    if r2_path:
        r2count = 0
        # read R2, write those that contained UMI in R1
        with pysam.FastxFile(r2_path) as r2_in:
            for n, entry in enumerate(r2_in):
                hasumi = True

                # check if paired R1 is deleted
                try:
                    x = empty_sequences[entry.name]
                    continue
                except:
                    pass

                # check for UMI
                if read_umis[entry.name]:
                    umi = read_umis[entry.name]
                    entry.name = entry.name + "_" + umi
                else:
                    hasumi = False

                if only_umi and not hasumi:
                    continue
                else:
                    r2_out.write((str(entry) + "\n").encode())
                    r2count += 1

    print(umicount, "umis, ", r1count, "reads written, ", len(empty_sequences.keys()), "reads removed")
    r1_out.close()
    if r2_path:
        r2_out.close()
