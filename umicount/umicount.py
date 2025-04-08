import sys
import os
import multiprocessing as mp
import argparse
import pickle
from dataclasses import dataclass, field
from typing import Any, Tuple, List

import HTSeq
from rapidfuzz.distance import Hamming

def validate_cols_to_use(cols):
    # because cols_to_use values direct downstream logic inplace of function args
    if cols is None: 
        print('malformed columns: columns set is empty')
        return False

    valid_fields = ['UI', 'UE', 'RI', 'RE', 'D']
    if any(col not in valid_fields for col in cols):
        print('malformed columns contains invalid values (expected %s)' %valid_fields)
        return False

    if 'RE' not in cols and 'RI' not in cols:
        print('malformed columns supplied, need at least RE, RI and (U or UE, UI')
        return False

    has_UIE = 'UI' in cols or 'UE' in cols
    has_U = 'U' in cols
    if has_UIE and has_U:
        print('malformed columns contains (UE, UI) and U')
        return False
    if has_UIE and not ('UI' in cols and 'UE' in cols):
        print('malformed columns contains UE or UI but not both')
        return False

    return True

def parse_gtf(gtffile, cols_to_use=None, 
              exon_attr=['gene_id', 'gene_name', 'exon_number'], 
              gene_attr=['gene_name']):

    assert validate_cols_to_use(cols_to_use)

    default = {i:0 for i in cols_to_use}
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

def load_gtf_data(gtffile, skipgtf=None, dumpgtf=None, cols_to_use=None):

    if skipgtf: # load GTF from pre-parsed dump
        print('Reading pre-parsed GTF data from:', skipgtf)
        return read_gtf_dump(skipgtf)

    else: # parse anew
        print('Parsing GTF file:', gtffile)
        assert validate_cols_to_use(cols_to_use) # check cols
        gtf_data = parse_gtf(gtffile, cols_to_use)
        if dumpgtf:
            print('Dumping parsed GTF data to:', dumpgtf)
            dump_gtf(dumpgtf, gtf_data)
        return gtf_data

@dataclass
class ReadTrack:
    # track information from a paired-end read
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
            # overlap loop as in HTSeq count.py
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

def umi_correction(umicounts, countratio=2, hamming_threshold=1):
        
    # UMIs sorted by decreasing counts
    umisort = sorted(umicounts.items(), key=lambda x: x[1], reverse=True)
    corrected = {}

    while umisort:
        seed, seed_count = umisort.pop(0) # highest counts UMI
        corrected[seed] = seed_count
        if seed_count <= 0: 
            raise ValueError(f"UMI ({seed}) has 0 or fewer counts")

        # iterate backwards from low to high UMI counts
        i = len(umisort) - 1
        while i >= 0:
            candidate, candidate_count = umisort[i]
            if candidate_count <= 0:
                raise ValueError(f"UMI ({seed}) has 0 or fewer counts")

            if ((countratio * candidate_count) - 1) > seed_count:
                break # can break since remaining UMIs counts are higher
                
            if Hamming.distance(seed, candidate, pad=False) <= hamming_threshold:
                corrected[seed] += candidate_count
                umisort.pop(i)

            i -= 1

    return corrected
            
def parse_bam_and_count(bamfile, gtf_data, cols_to_use=None, umi_correct_params=None):

    assert validate_cols_to_use(cols_to_use)
    combine_unspliced = ('U' in cols_to_use)
    do_dedup = ('D' in cols_to_use)

    gcounts, gfeatures, efeatures, gattributes, eattributes = gtf_data
    bam_reader = HTSeq.BAM_Reader(bamfile)
    geneumis = {i:({'U':{}} if combine_unspliced else {'UE':{}, 'UI':{}}) \
                for i in gcounts.keys() if not i.startswith('_')} # prep corr and dedup

    # iterate over reads and assign feature overlaps
    for bundle in HTSeq.pair_SAM_alignments(bam_reader, bundle=True):
        if not bundle:
            continue

        readpair = extract_first_alignment(bundle)
        if len(bundle) > 1:
            readpair.category = '_multimapping'

        # extract overlapping genes, exons
        if readpair.can_do_overlap(): readpair.find_overlap(gfeatures, efeatures)
        if readpair.can_do_overlap(): readpair.evaluate_overlap(eattributes)

        # now count read at appropriate category based on overlap
        if readpair.gene_to_count == "": # no gene overlap
            gcounts[readpair.category]['UE' if readpair.umi else 'RE'] += 1 

        else: # have gene-counts
            if readpair.umi:
                rkey = 'U' if combine_unspliced else \
                       ('UI' if readpair.exon_to_count == "" else 'UE')
                rgene = readpair.gene_to_count
                rumi = readpair.umi
                geneumis[rgene][rkey][rumi] = geneumis[rgene][rkey].get(rumi, 0) + 1

            else:
                rkey = 'RI' if readpair.exon_to_count == "" else 'RE'
                gcounts[readpair.gene_to_count][rkey] += 1

    # parsed all reads into gene overlaps, now process UMI counts
    if umi_correct_params is None: # no UMI correction
        for g in geneumis.keys():
            if combine_unspliced:
                umisum = sum(geneumis[g]['U'].values())
                umilen = len(geneumis[g]['U'])

                if do_dedup:
                    gcounts[g]['U'] = umilen
                    gcounts[g]['D'] = umisum - umilen
                else:
                    gcounts[g]['U'] = umisum
            else: 
                for uie in ['UI', 'UE']:
                    umisum = sum(geneumis[g][uie].values())
                    umilen = len(geneumis[g][uie])

                    gcounts[g][uie] = umilen if do_dedup else umisum

                if do_dedup:
                    gcounts[g]['D'] = sum(geneumis[g]['UI'].values()) - gcounts[g]['UI'] + \
                                      sum(geneumis[g]['UE'].values()) - gcounts[g]['UE']

    else:
        # counts with UMI correction
        countratio_threshold = umi_correct_params['countratio_threshold']
        hamming_threshold = umi_correct_params['hamming_threshold']

        for g in geneumis.keys():
            if combine_unspliced:
                UIE_corrected = umi_correction(geneumis[g]['U'], 
                                               countratio=countratio_threshold, 
                                               hamming_threshold=hamming_threshold)

                if do_dedup:
                    gcounts[g]['U'] = len(UIE_corrected.keys()) # unique UMIs
                    gcounts[g]['D'] = sum(UIE_corrected.values()) - gcounts[g]['U']
                else:
                    gcounts[g]['U'] = sum(UIE_corrected.values()) # total UMI counts

            else:
                UI_corrected = umi_correction(geneumis[g]['UI'], 
                                              countratio=countratio_threshold, 
                                              hamming_threshold=hamming_threshold)
                UE_corrected = umi_correction(geneumis[g]['UE'], 
                                              countratio=countratio_threshold, 
                                              hamming_threshold=hamming_threshold)

                if do_dedup:
                    gcounts[g]['UI'] = len(UI_corrected.keys())
                    gcounts[g]['UE'] = len(UE_corrected.keys())
                    gcounts[g]['D'] = sum(UI_corrected.values()) - gcounts[g]['UI'] + \
                                      sum(UE_corrected.values()) - gcounts[g]['UE']
                else:
                    gcounts[g]['UI'] = sum(UI_corrected.values())
                    gcounts[g]['UE'] = sum(UE_corrected.values())
    return gcounts

def write_counts(outfile, bamfile, results, gene_counts, gattributes, cols_to_use):
    assert validate_cols_to_use(cols_to_use)
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

def process_bam(bamfile, gtffile, outfile, skipgtf=None, 
                cols_to_use=None, umi_correct_params=None):

    assert validate_cols_to_use(cols_to_use)

    # Load or parse the GTF data
    gtf_data = load_gtf_data(gtffile, skipgtf=skipgtf, dumpgtf=None, cols=cols_to_use)
    gcounts, gfeatures, efeatures, gattributes, eattributes = gtf_data

    results = parse_bam_and_count(bamfile, gtf_data, 
                                  cols_to_use=cols_to_use, 
                                  umi_correct_params=umi_correct_params)

    write_counts(outfile, bamfile, results, gcounts, gattributes, cols_to_use)
