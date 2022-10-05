import sys
import os
import multiprocessing as mp
import copy
import argparse
import pickle
from dataclasses import dataclass, field
from typing import Any

import HTSeq

BASECOLS = ['UE', 'RE', 'UI', 'RI', 'D']

def parse_gtf(gtffile, cols, exon_attr=['gene_id', 'gene_name', 'exon_number'], gene_attr=['gene_name']):
    default = {i:0 for i in cols}
    gcounts = {'_unmapped':default.copy(), '_multimapping':default.copy(),
               '_no_feature':default.copy(), '_ambiguous':default.copy(), '_FAIL':default.copy()}
    gattributes = {}; eattributes = {}
    additional_attributes = {'exon':exon_attr, 'gene':gene_attr}

    gfeatures = HTSeq.GenomicArrayOfSets("auto", False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", False)

    # parse GTF as in HTSeq count.py
    gff = HTSeq.GFF_Reader(gtffile)
    try:
        for f in gff:
            if f.type == 'exon':
                id_attribute = 'exon_id'
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError("Feature %s does not contain a '%s' attribute" %(f.name, id_attribute))

                efeatures[f.iv] += feature_id
                eattributes[f.attr[id_attribute]] = [f.attr[attr] if attr in f.attr else
                                                     '' for attr in additional_attributes[f.type]]

            if f.type == 'gene':
                id_attribute = 'gene_id'
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError("Feature %s does not contain a '%s' attribute" %(f.name, id_attribute))

                gfeatures[f.iv] += feature_id
                gcounts[f.attr[id_attribute]] = default.copy() # store raw/umi for gene counts, intron, and exon
                gattributes[f.attr[id_attribute]] = [f.attr[attr] if attr in f.attr else
                                                     '' for attr in additional_attributes[f.type]]

    except:
        sys.stderr.write(
            "Error occured when processing GFF file (%s):\n" %
            gff.get_line_number_string())
        raise

    return (gcounts, gfeatures, efeatures, gattributes, eattributes)

def dump_gtf(dumppath, items):
    with open(dumppath, 'wb') as out:
        pickle.dump(items, out)

def read_gtf_dump(dumppath):
    with open(dumppath, 'rb') as input:
        return pickle.load(input)

@dataclass
class ReadTrack:
    # track information from a PE read
    read1_almnt: Any
    read2_almnt: Any
    umi: Any = field(init=False)
    category: str = ""
    gene_overlap: list = field(default_factory=list)
    exon_overlap: list = field(default_factory=list)
    gene_to_count: str = ""
    exon_to_count: str = ""

    def __post_init__(self):
        if len(self.read1_almnt.read.name.split('_')) > 1:
            self.umi = self.read1_almnt.read.name.split('_')[-1]
        else:
            self.umi = None

    def can_do_overlap(self):
        if self.read1_almnt and self.read2_almnt and (self.category == ''):
            return True
        return False

    def find_overlap(self, gfeatures, efeatures):
        assert self.can_do_overlap()

        gene_ids = exon_ids = set()
        try:
            for almnt in [self.read1_almnt, self.read2_almnt]:
                for iv, val in gfeatures[almnt.iv].steps():
                    gene_ids = gene_ids.union(val)
                self.gene_overlap = list(gene_ids) # ids of overlapping genes

                for iv, val in efeatures[almnt.iv].steps():
                    exon_ids = exon_ids.union(val)
                self.exon_overlap = list(exon_ids) # ids of overlapping exons
        except KeyError:
            # weird bug where features doesnt contain scaffold chr
            self.category = '_FAIL'

        return self

    def evaluate_overlap(self, eattributes):
        assert self.can_do_overlap()

        if len(self.gene_overlap) == 0: # intergenic
            self.category = '_no_feature'

        elif len(self.gene_overlap) == 1:
            self.gene_to_count = self.gene_overlap[0]

            # 1 gene, how many exons ? (we skip intron checking since exon_to_count default is empty)
            if len(self.exon_overlap) > 0:
                self.exon_to_count = self.exon_overlap[0]

        else: # when multiple genes overlap
            exongenes = set([eattributes[i][0] for i in self.exon_overlap]) # gene ids from attribute
            if len(exongenes) == 1: # if exon annotated for only 1 of overlapping genes: count it
                self.gene_to_count = self.exon_to_count = list(exongenes)[0]
            else: # if exons overlap or intronic in both: its ambiguous
                self.category = '_ambiguous'

        return self

def extract_first_alignment(bundle):
    firstpair = bundle[0] # consider first reported alignment by default

    if None in firstpair:
        return ReadTrack(read1_almnt=[i for i in firstpair if i][0],
                         read2_almnt=None,
                         category='_unmapped')

    if not (firstpair[0].aligned and firstpair[1].aligned):
        return ReadTrack(read1_almnt=firstpair[0],
                         read2_almnt=None,
                         category='_unmapped')

    return ReadTrack(read1_almnt=firstpair[0], read2_almnt=firstpair[1])

def parse_bam_and_count(bamfile, skipgtf):
    read_seq_file = HTSeq.BAM_Reader(bamfile)
    gcounts, gfeatures, efeatures, gattributes, eattributes = read_gtf_dump(skipgtf)

    umicheck = {i:{} for i in gcounts.keys()} # for tracking duplicates
    for bundle in HTSeq.pair_SAM_alignments(read_seq_file, bundle=True):
        if not bundle:
            continue

        # check for unmapped, multimapped reads
        readpair = extract_first_alignment(bundle)
        if len(bundle) > 1:
            readpair.category = '_multimapping'

        # extract overlapping genes, exons
        if readpair.can_do_overlap(): readpair.find_overlap(gfeatures, efeatures)
        if readpair.can_do_overlap(): readpair.evaluate_overlap(eattributes)

        # add to counts for cell
        countskey = 'U' if readpair.umi else 'R'
        if readpair.gene_to_count == "":
            # non-gene-counts default in exon column
            gcounts[readpair.category][countskey+'E'] += 1
        else:
            # have gene-counts, so first deduplicate
            if readpair.umi:
                try:
                    umicheck[readpair.gene_to_count][readpair.umi] += 1
                    continue
                except KeyError:
                    umicheck[readpair.gene_to_count][readpair.umi] = 1

            if readpair.exon_to_count == "":
                # gene exists, no exon --> is intron
                gcounts[readpair.gene_to_count][countskey+'I'] += 1

            else:
                # gene exists, has exon --> is exon
                gcounts[readpair.gene_to_count][countskey+'E'] += 1

    # report umi duplicate counts for cell
    for i in umicheck:
        if i.startswith('_'): continue
        gcounts[i]['D'] = sum(umicheck[i].values()) - (gcounts[i]['UI'] + gcounts[i]['UE'])

    return gcounts

def process_bam(infiles, gtffile, outfile, numcores, dumpgtf, skipdup, colnames, outgeneformat):
    # parse gtf and get data structures
    global BASECOLS
    if not skipgtf:
        # parse and dump gtf contents
        print('parsing GTF file')
        gcounts, gfeatures, efeatures, gattributes, eattributes = parse_gtf(gtffile, BASECOLS)

        print('dumping GTF contents')
        dump_gtf(dumpgtf, (gcounts, gfeatures, efeatures, gattributes, eattributes))
        if len(infiles) == 0:
            sys.exit()

        skipgtf = dumpgtf
    else:
        print('reading from --skipgtf file')
        gcounts, gfeatures, efeatures, gattributes, eattributes = read_gtf_dump(skipgtf)

    # parse bam and count reads / UMIs
    print('parsing BAM')
    if len(infiles) == 1:
        results = [parse_bam_and_count(infiles[0], skipgtf)]
    else:
        pool = mp.Pool(numcores)
        results = pool.starmap(parse_bam_and_count, [(f, skipgtf) for f in infiles])
        pool.close()

    # create counts matrix
    print('writing counts')
    if skipdup:
        BASECOLS = [i for i in BASECOLS if i != 'D']

    with open(outfile, 'w') as w:
        # header
        h = 'gene'
        for i in colnames:
            for j in BASECOLS:
                h += '\t' + j + '_' + os.path.basename(i)
        w.write(h + '\n')

        # rest
        for i in gcounts.keys():
            if i.startswith('_'):
                line = str(i)
            else:
                line = str(i) if outgeneformat == 'geneid' else str(gattributes[i][0])
            for n, j in enumerate(infiles):
                for k in BASECOLS:
                    line += '\t' + str(results[n][i][k])
            w.write(line + '\n')
