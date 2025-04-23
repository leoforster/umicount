import sys
import os
import pickle
from multiprocessing import Pool
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Tuple, List

import HTSeq
try:
    from rapidfuzz.distance import Hamming
except: # this is checked in cli.py
    pass

def validate_cols_to_use(cols):
    # because cols_to_use values direct downstream logic inplace of function args
    if cols is None: 
        print('malformed columns: columns set is empty')
        return False

    valid_fields = ['UI', 'UE', 'RI', 'RE', 'D', 'R', 'U']
    if any(col not in valid_fields for col in cols):
        print('malformed columns contains invalid values (expected %s)' %valid_fields)
        return False

    has_RIE = 'RI' in cols or 'RE' in cols
    has_R = 'R' in cols
    if has_RIE and has_R:
        print('malformed columns contains (RE, RI) and R')
        return False
    if has_RIE and not ('RI' in cols and 'RE' in cols):
        print('malformed columns contains RE or RI but not both')
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

def parse_gtf(gtffile, exon_attr=['gene_id', 'gene_name', 'exon_number'], gene_attr=['gene_name']):

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
                gattributes[feature_id] = [f.attr.get(attr, '') \
                                           for attr in additional_attributes[f.type]]
    except:
        sys.stderr.write(
            "Error occured when processing GFF file (%s):\n" %gff.get_line_number_string())
        raise

    return (gfeatures, efeatures, gattributes, eattributes)

def dump_gtf(dump_path, items):
    with open(dump_path, 'wb') as out:
        pickle.dump(items, out)

def read_gtf_dump(dump_path):
    with open(dump_path, 'rb') as inp:
        return pickle.load(inp)

def load_gtf_data(gtffile, skipgtf=None, dumpgtf=None):

    if skipgtf: # load GTF from pre-parsed dump
        return read_gtf_dump(skipgtf)

    else: # parse anew
        print('Parsing GTF file:', gtffile)
        gtf_data = parse_gtf(gtffile)
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
        if self.read1_almnt:
            parts = self.read1_almnt.read.name.rsplit('_', 1)
        elif self.read2_almnt:
            parts = self.read2_almnt.read.name.rsplit('_', 1)
        else:
            parts = None

        self.umi = parts[1] if len(parts) == 2 and parts[1] else None

    def can_do_overlap(self):
        return self.read1_almnt is not None and \
               self.read2_almnt is not None and (self.category == '')

    def find_overlap(self, gfeatures, efeatures):
        assert self.can_do_overlap()

        gene_ids = set()
        exon_ids = set()
        try:
            # overlap loop as in HTSeq count.py
            for almnt in [self.read1_almnt, self.read2_almnt]:
                if almnt is None: continue

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
            if len(self.exon_overlap) > 0: # exons of that gene
                self.exon_to_count = self.exon_overlap[0]

        else: # when multiple genes overlap
            exongenes = {eattributes[eid][0] for eid in self.exon_overlap \
                         if eid in eattributes} # get gene ids from attribute

            # if exon annotated for only 1 of overlapping genes: count it
            if len(exongenes) == 1:
                self.gene_to_count = self.exon_to_count = exongenes.pop()

            # if exons overlap or intronic in both: its ambiguous
            else:
                self.category = '_ambiguous'

        return self

def extract_first_alignment(bundle, count_primary=False):

    if not bundle: raise ValueError("empty bundle: no alignments to extract")

    r1_to_count, r2_to_count = bundle[0] # default to first
    read_category = '' # default value as in ReadTrack definition

    # multimapping readpair
    if len(bundle) > 1: 

        if count_primary:
            primaries = []
            for r1, r2 in bundle:
                p1 = r1 if (r1 and r1.not_primary_alignment is False) else None
                p2 = r2 if (r2 and r2.not_primary_alignment is False) else None
                if p1 or p2: # at least one of r1/r2 has primary flag
                    primaries.append((p1, p2))

            if len(primaries) > 1:
                raise ValueError("found multiple primary alignments:\n%s" %bundle)

            if not primaries: # no primary alignments (or not set by aligner)
                read_category = '_multimapping'
            else:
                r1_to_count, r2_to_count = primaries[0]

        else: # multimapping but no count_primary: use first pair
            read_category = '_multimapping'

    # readpair has single alignment
    else:
        if (r1_to_count is None or not r1_to_count.aligned) or \
           (r2_to_count is None or not r2_to_count.aligned): # either read unaligned
            read_category = '_unmapped'

    return ReadTrack(read1_almnt=r1_to_count,
                     read2_almnt=r2_to_count,
                     category=read_category)

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
            
def parse_bam_and_count(bamfile, gtf_data, 
                        cols_to_use=None, 
                        count_primary=False,
                        umi_correct_params=None):

    assert validate_cols_to_use(cols_to_use)
    combine_unspliced = 'U' in cols_to_use
    do_dedup = 'D' in cols_to_use

    # assign GTF data and fill default columns
    gfeatures, efeatures, gattributes, eattributes = gtf_data
    gcounts = {}
    category_cols = ['_unmapped', '_multimapping', '_no_feature', '_ambiguous', '_FAIL']
    for g in category_cols + list(gattributes.keys()):
        gcounts[g] = defaultdict(int)

    # read and UMI count tracking dicts
    geneumis = {i:defaultdict(lambda: defaultdict(int)) for i in gattributes.keys()} 
    totalumis = {i:0 for i in cols_to_use + ['total', 'uncounted', 'corrected']} 

    # iterate over reads and assign feature overlaps
    bam_reader = HTSeq.BAM_Reader(bamfile)
    for bundle in HTSeq.pair_SAM_alignments(bam_reader, bundle=True):
        totalumis['total'] += 1

        if not bundle:
            totalumis['uncounted'] += 1
            continue

        readpair = extract_first_alignment(bundle, count_primary=count_primary)

        # extract overlapping genes, exons
        if readpair.can_do_overlap(): readpair.find_overlap(gfeatures, efeatures)
        if readpair.can_do_overlap(): readpair.evaluate_overlap(eattributes)

        # now count read at appropriate category based on overlap
        rkey = 'U' if readpair.umi else 'R'
        if readpair.gene_to_count == '': # no gene overlap
            gcounts[readpair.category][(rkey + 'E') if not combine_unspliced else rkey] += 1
            totalumis['uncounted'] += 1

        else: # has gene overlap
            if not combine_unspliced:
                rkey += 'I' if readpair.exon_to_count == '' else 'E'

            if readpair.umi: # store gene UMIs separately for deduplication/correction
                geneumis[readpair.gene_to_count][rkey][readpair.umi] += 1
            else:
                gcounts[readpair.gene_to_count][rkey] += 1

    # parsed all reads into gene overlaps, now process UMI counts
    if umi_correct_params is None: # no UMI correction
        if combine_unspliced:
            for g in geneumis.keys():
                umisum = sum(geneumis[g]['U'].values())
                umilen = len(geneumis[g]['U'])

                if do_dedup:
                    gcounts[g]['U'] = umilen
                    gcounts[g]['D'] = umisum - umilen
                else:
                    gcounts[g]['U'] = umisum
        else: 
            for g in geneumis.keys():
                for uie in ['UI', 'UE']:
                    umisum = sum(geneumis[g][uie].values())
                    umilen = len(geneumis[g][uie])

                    gcounts[g][uie] = umilen if do_dedup else umisum

                if do_dedup:
                    gcounts[g]['D'] = sum(geneumis[g]['UI'].values()) - gcounts[g]['UI'] + \
                                      sum(geneumis[g]['UE'].values()) - gcounts[g]['UE']

    else: # with UMI correction
        ct = umi_correct_params['countratio_threshold']
        ht = umi_correct_params['hamming_threshold']
        
        if combine_unspliced:
            for g in geneumis.keys():
                UIE_corrected = umi_correction(geneumis[g]['U'], countratio=ct, hamming_threshold=ht)

                # track sum of counts of corrected UMIs
                totalumis['corrected'] += sum([geneumis[g]['U'][i] for i in geneumis[g]['U'] \
                                               if i not in UIE_corrected.keys()])

                if do_dedup:
                    gcounts[g]['U'] = len(UIE_corrected.keys()) # unique UMIs
                    gcounts[g]['D'] = sum(UIE_corrected.values()) - gcounts[g]['U']
                else:
                    gcounts[g]['U'] = sum(UIE_corrected.values()) # total UMI counts

        else:
            for g in geneumis.keys():
                UI_corrected = umi_correction(geneumis[g]['UI'], countratio=ct, hamming_threshold=ht)
                UE_corrected = umi_correction(geneumis[g]['UE'], countratio=ct, hamming_threshold=ht)

                # track sum of counts of corrected UMIs
                totalumis['corrected'] += sum([geneumis[g]['UI'][i] for i in geneumis[g]['UI'] \
                                               if i not in UI_corrected.keys()])
                totalumis['corrected'] += sum([geneumis[g]['UE'][i] for i in geneumis[g]['UE'] \
                                               if i not in UE_corrected.keys()])

                if do_dedup:
                    gcounts[g]['UI'] = len(UI_corrected.keys())
                    gcounts[g]['UE'] = len(UE_corrected.keys())
                    gcounts[g]['D'] = sum(UI_corrected.values()) - gcounts[g]['UI'] + \
                                      sum(UE_corrected.values()) - gcounts[g]['UE']
                else:
                    gcounts[g]['UI'] = sum(UI_corrected.values())
                    gcounts[g]['UE'] = sum(UE_corrected.values())

    # update total counts
    for g in gattributes.keys():
        for i in cols_to_use:
            totalumis[i] += gcounts[g][i]

    return gcounts, totalumis

def write_counts(outfile, bamfile, gene_counts, gattributes, cols_to_use):
    assert validate_cols_to_use(cols_to_use)

    with open(outfile, 'w') as out:
        header_fields = [os.path.basename(bamfile)]
        for b in cols_to_use:
            header_fields.append(b)
        out.write('\t'.join(header_fields) + '\n')
        for gene, counts in gene_counts.items():
            line_fields = [gene]
            for b in cols_to_use:
                line_fields.append(str(gene_counts[gene][b]))
            out.write('\t'.join(line_fields) + '\n')

def process_bam(bamfile, outfile, gtf_data,
                cols_to_use=None, 
                count_primary=False,
                umi_correct_params=None):

    assert validate_cols_to_use(cols_to_use)

    # assign GTF data
    gfeatures, efeatures, gattributes, eattributes = gtf_data

    # parsing BAM and count reads
    umicount, ttl = parse_bam_and_count(bamfile, gtf_data, 
                                        cols_to_use=cols_to_use, 
                                        count_primary=count_primary,
                                        umi_correct_params=umi_correct_params)

    # output counts table to file
    write_counts(outfile, bamfile, umicount, gattributes, cols_to_use)

    # print summarized quantities of counted read categories
    if sum(ttl.values()) > 0:
        sumstr = f"{os.path.basename(bamfile)}: {ttl['total']} reads"
        sumstr += f", {ttl['uncounted']} uncounted reads ({(ttl['uncounted']/ttl['total'])*100:.2f}%)"
        for i in cols_to_use:
            sumstr += f", {ttl[i]} {i}-reads ({(ttl[i]/ttl['total'])*100:.2f}%)"
        if umi_correct_params: sumstr += f", {ttl['corrected']} counts in corrected UMIs "
        print(sumstr)
    else:
        print(f"empty BAM file")

def _bam_worker(task_args):
    # worker function for multiprocessing
    infile, outfile, gtf_data, kwargs = task_args
    process_bam(infile, outfile, gtf_data, **kwargs)

def process_bam_parallel(filepairs, gtf_data, num_workers=4,
                         cols_to_use=None, 
                         count_primary=False,
                         umi_correct_params=None):

    # wrapper around process_bam to handle multiple BAM files in parallel
    # filepairs are a list of tasks as tuples of ( bamfile, outfile )
    # gtf_data can be parsed from GTF file using load_gtf_data(gtffile)

    # bundle constant args
    kwargs = dict(
        cols_to_use=cols_to_use,
        count_primary=count_primary,
        umi_correct_params=umi_correct_params
    )

    # map the worker over the filepairs
    tasks = [(infile, outfile, gtf_data, kwargs) for (infile, outfile) in filepairs]
    with Pool(num_workers) as pool:
        results = pool.map(_bam_worker, tasks)
