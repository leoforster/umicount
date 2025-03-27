import sys
import os
import multiprocessing as mp
import argparse
import pickle
from dataclasses import dataclass, field
from typing import Any, Tuple, List

import HTSeq

BASECOLS = ['UE', 'RE', 'UI', 'RI', 'D']

def parse_gtf(gtffile, cols, 
              exon_attr=['gene_id', 'gene_name', 'exon_number'], 
              gene_attr=['gene_name']):

    default = {i:0 for i in cols}
    gcounts = {'_unmapped':default.copy(), '_multimapping':default.copy(),
               '_no_feature':default.copy(), '_ambiguous':default.copy(), '_FAIL':default.copy()}
    gattributes = {}
    eattributes = {}
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
                    raise ValueError("%s does not contain a '%s' attribute" %(f.name, id_attribute))

                efeatures[f.iv] += feature_id
                eattributes[f.attr[id_attribute]] = [f.attr.get(attr, '') \
                                                     for attr in additional_attributes[f.type]]

            if f.type == 'gene':
                id_attribute = 'gene_id'
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError("%s does not contain a '%s' attribute" %(f.name, id_attribute))

                # store raw/umi for gene counts, intron, and exon
                gfeatures[f.iv] += feature_id
                gcounts[feature_id] = default.copy() 
                gattributes[feature_id] = [f.attr.get(attr, '') \
                                           for attr in additional_attributes[f.type]]
    except:
        sys.stderr.write(
            "Error occured when processing GFF file (%s):\n" %gff.get_line_number_string())
        raise

    return (gcounts, gfeatures, efeatures, gattributes, eattributes)

def dump_gtf(dump_path, items):
    with open(dump_path, 'wb') as out:
        pickle.dump(items, out)

def read_gtf_dump(dump_path):
    with open(dump_path, 'rb') as inp:
        return pickle.load(inp)

def load_gtf_data(gtffile, 
                  skipgtf=None, dumpgtf=None,
                  cols=BASECOLS):
    """
    Load GTF data either by reading a pre-parsed dump (if skipgtf is provided)
    or by parsing the GTF file. Optionally, dump the parsed data.
    """
    if skipgtf:
        print('Reading pre-parsed GTF data from:', skipgtf)
        return read_gtf_dump(skipgtf)
    else:
        print('Parsing GTF file:', gtffile)
        gtf_data = parse_gtf(gtffile, cols)
        if dumpgtf:
            print('Dumping parsed GTF data to:', dumpgtf)
            dump_gtf(dumpgtf, gtf_data)
        return gtf_data

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
        parts = self.read1_almnt.read.name.rsplit('_', 1)
        self.umi = parts[1] if len(parts) == 2 and parts[1] else None

    def can_do_overlap(self):
        return self.read1_almnt is not None and \
               self.read2_almnt is not None and (self.category == '')

    def find_overlap(self, gfeatures, efeatures):
        if not self.can_do_overlap():
            return self

        gene_ids = set()
        exon_ids = set()
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
        if not self.can_do_overlap():
            return self

        if len(self.gene_overlap) == 0: # intergenic
            self.category = '_no_feature'
        elif len(self.gene_overlap) == 1:
            self.gene_to_count = self.gene_overlap[0]
            if len(self.exon_overlap) > 0: # exons of that gene
                self.exon_to_count = self.exon_overlap[0]
        else: # when multiple genes overlap
            exongenes = {eattributes[eid][0] for eid in self.exon_overlap \
                         if eid in eattributes} # get gene ids from attribute
            if len(exongenes) == 1: # if exon annotated for only 1 of overlapping genes: count it
                self.gene_to_count = self.exon_to_count = exongenes.pop()
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

def parse_bam_and_count(bamfile, gtf_data):
    gcounts, gfeatures, efeatures, gattributes, eattributes = gtf_data
    bam_reader = HTSeq.BAM_Reader(bamfile)
    umicheck = {key: {} for key in gcounts.keys()} # for tracking duplicates

    for bundle in HTSeq.pair_SAM_alignments(bam_reader, bundle=True):
        if not bundle:
            continue

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
            if readpair.umi: # have gene-counts, so first deduplicate
                if readpair.umi in umicheck[readpair.gene_to_count]:
                    umicheck[readpair.gene_to_count][readpair.umi] += 1
                    continue
                else:
                    umicheck[readpair.gene_to_count][readpair.umi] = 1

            if readpair.exon_to_count == "": # gene exists, no exon --> is intron
                gcounts[readpair.gene_to_count][countskey+'I'] += 1
            else: # gene exists, has exon --> is exon
                gcounts[readpair.gene_to_count][countskey+'E'] += 1

    # report umi duplicate counts for cell
    for gene, umi_dict in umicheck.items():
        if gene.startswith('_'):
            continue
        gcounts[gene]['D'] = sum(umi_dict.values()) - (gcounts[gene]['UI'] + gcounts[gene]['UE'])
    return gcounts

def write_counts(outfile, bamfile, results, gene_counts, gattributes, cols_to_use):
    with open(outfile, 'w') as out:
        header_fields = [os.path.basename(bamfile)]
        for b in cols_to_use:
            header_fields.append(b)
        out.write('\t'.join(header_fields) + '\n')
        for gene, counts in gene_counts.items():
            line_fields = [gene]
            for b in cols_to_use:
                line_fields.append(str(results[gene][b]))
            out.write('\t'.join(line_fields) + '\n')

def process_bam(bamfile, gtffile, outfile, 
                skipgtf=None, skipdup=False):
    # Load or parse the GTF data
    gtf_data = load_gtf_data(gtffile, skipgtf=skipgtf, dumpgtf=None, cols=BASECOLS)
    gcounts, gfeatures, efeatures, gattributes, eattributes = gtf_data

    results = parse_bam_and_count(bamfile, gtf_data)

    cols_to_use = BASECOLS.copy()
    if skipdup:
        cols_to_use = [col for col in cols_to_use if col != 'D']

    write_counts(outfile, bamfile, results, gcounts, gattributes, cols_to_use)
