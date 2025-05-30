import pytest
import HTSeq

from umicount.umicount import (
    parse_gtf,
    ReadTrack,
    set_alignment_category,
    filter_aligned_reads,
    parse_bam_and_count
)

############################
# Dummy HTSeq alignment class
############################

class DummyAlignment:
    """
    A dummy alignment object to simulate an HTSeq alignment.
    The read name encodes the UMI (if an underscore exists).
    """
    def __init__(self, name, aligned=True, not_primary_alignment=False, iv=None, qual=0):
        DummyRead = type("DummyRead", (), {})  # create an empty dummy class
        self.read = DummyRead()
        self.read.name = name
        self.aligned = aligned
        self.not_primary_alignment = not_primary_alignment
        self.iv = iv
        self.aQual = qual

def dummy_pair_SAM_alignments_factory(bundles):
    """
    Returns a generator that yields the supplied bundles one by one.
    This makes it easier to write tests that need only a couple of bundles.
    """
    def generator(reader, bundle=True):
        for b in bundles:
            yield b
    return generator

def get_dummy_gtf_dump(genes_intervals, cols_to_use=['UE', 'RE', 'UI', 'RI', 'D']):
    """
    Build a dummy GTF dump to be used by parse_bam_and_count.
    genes_intervals is a list of tuples:
       (gene_id, gene_name, exon_id, interval, include_exon)

    """
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gattributes = {}
    eattributes = {}
    for gene_id, gene_name, exon_id, interval, include_exon in genes_intervals:
        gattributes[gene_id] = [gene_name]
        gfeatures[interval] += gene_id
        if include_exon:
            efeatures[interval] += exon_id
            eattributes[exon_id] = [gene_id, gene_name, "1"]
    return (gfeatures, efeatures, gattributes, eattributes)

def create_readtrack(name, iv, aligned=True, not_primary_alignment=False, qual=0):
    """
    Helper to create a ReadTrack from two dummy alignments.
    """
    aln1 = DummyAlignment(name, aligned, not_primary_alignment, iv, qual)
    aln2 = DummyAlignment(name, aligned, not_primary_alignment, iv, qual)
    return ReadTrack(aln1, aln2)

############################
# Tests for umicount
############################

def test_parse_gtf(tmp_path):
    """
    Test the HTSeq GTF parsing function using a minimal valid GTF file.
    """
    gtf_content = (
        'chr1\tsource\tgene\t100\t200\t.\t+\t.\tgene_id "gene1"; gene_name "GeneName";\n'
        'chr1\tsource\texon\t100\t200\t.\t+\t.\texon_id "exon1"; gene_id "gene1"; gene_name "GeneName"; exon_number "1";\n'
    )
    gtf_file = tmp_path / "test.gtf"
    gtf_file.write_text(gtf_content)
    result = parse_gtf(str(gtf_file))
    gfeatures, efeatures, gattributes, eattributes = result

    # Check that gene1 is in the counts and attributes
    assert "gene1" in gattributes

    # Check that interval is present in gfeatures and efeatures
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    found_gene = False
    for iv2, val in gfeatures[iv].steps():
        if "gene1" in val:
            found_gene = True
    assert found_gene

def test_parse_malformed_gtf(tmp_path):
    """
    Test the HTSeq GTF parsing function on a malformed GTF.
    """
    gtf_content = (
        'chr1\tsourcegene\t100\t200\t.\t+\t.\tgene_id "gene1"; gene_name "GeneName";\n'
        'chr1\tsource\texon\t100\t200\t.\n+\t+\t.\texon_id "exon1";\n'
    )
    gtf_file = tmp_path / "test.gtf"
    gtf_file.write_text(gtf_content)
    
    with pytest.raises(ValueError):
        result = parse_gtf(str(gtf_file))

def test_filter_aligned_reads():
    """
    Test return of filter_aligned_reads.
    """
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    # Case 1: single readpair with mapQ pass.
    aln1 = DummyAlignment("readC_UMI", True, False, iv, 0)
    aln2 = DummyAlignment("readC_UMI", True, False, iv, 0)
    bundle = [(aln1, aln2)]
    bd = filter_aligned_reads(bundle, min_read_mapQ=0)
    assert bd == bundle

    # Case 2: multiple readpair with mapQ pass.
    aln3 = DummyAlignment("readC_UMI", True, False, iv, 0)
    aln4 = DummyAlignment("readC_UMI", True, False, iv, 0)
    bundle = [(aln1, aln2), (aln3, aln4)]
    bd = filter_aligned_reads(bundle, min_read_mapQ=0)
    assert bd == bundle

    # Case 3: multiple readpair where one fails mapQ.
    aln3 = DummyAlignment("readC_UMI", True, False, iv, 3)
    aln4 = DummyAlignment("readC_UMI", True, False, iv, 3)
    bundle = [(aln1, aln2), (aln3, aln4)]
    bd = filter_aligned_reads(bundle, min_read_mapQ=1)
    assert bd == [(aln3, aln4)]

    # Case 4: readpair fails mapQ.
    bundle = [(aln1, aln1)] # these fail as mapQ==0
    bd = filter_aligned_reads(bundle, min_read_mapQ=1)
    assert bd == []

    # Case 5: like Case 4 but checking ReadTrack.cateogory 
    rt = set_alignment_category(bundle, bamfile="test", min_read_mapQ=1)
    assert rt.read1_almnt is not None
    assert rt.read2_almnt is not None
    assert rt.category == '_unmapped'

def test_set_alignment_category():
    """
    Test set_alignment_category for different bundle cases.
    """
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    # Case 1: proper pair.
    aln1 = DummyAlignment("readC_UMI", True, False, iv, 0)
    aln2 = DummyAlignment("readC_UMI", True, False, iv, 0)
    bundle = [(aln1, aln2)]
    rt = set_alignment_category(bundle, bamfile="test")
    assert rt.category == ""
    assert rt.read1_almnt.read.name == "readC_UMI"
    assert rt.read2_almnt is not None

    # Case 2: one alignment is None (simulate unmapped).
    bundle_unmapped = [(aln1, None)]
    rt2 = set_alignment_category(bundle_unmapped, bamfile="test")
    assert rt2.category == "_unmapped"
    assert rt2.read2_almnt is None

    # Case 3: both alignments exist but one is not aligned.
    aln3 = DummyAlignment("readD_UMI", False, True, iv, 0)
    bundle_unmapped2 = [(aln2, aln3)]
    rt3 = set_alignment_category(bundle_unmapped2, bamfile="test")
    assert rt3.category == "_unmapped"

    # Case 4: empty bundle
    with pytest.raises(ValueError):
        set_alignment_category([], bamfile="test")

    # Case 5: multimapping but count_primary not set
    aln4 = DummyAlignment("readE_UMI", True, True, iv, 0)
    aln5 = DummyAlignment("readE_UMI", True, True, iv, 0)
    aln6 = DummyAlignment("readF_UMI", True, True, iv, 0)
    aln7 = DummyAlignment("readF_UMI", True, True, iv, 0)
    bundle_multi = [(aln4, aln5), (aln6, aln7)]
    rt5 = set_alignment_category(bundle_multi, bamfile="test", count_primary=False)
    assert rt5.category == "_multimapping"
    assert rt5.read1_almnt is aln4
    assert rt5.read2_almnt is aln5

    # Case 6: multimapping with no primary alignment
    bundle_multi = [(aln4, aln5), (aln6, aln7)]
    rt6 = set_alignment_category(bundle_multi, bamfile="test", count_primary=True)
    assert rt6.category == "_multimapping"

def test_multimapping_primary_alignments():
    """
    Test set_alignment_category for cases involving primary alignments.
    """
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")

    aln1 = DummyAlignment("readE_UMI", True, True, iv, 0)
    aln2 = DummyAlignment("readE_UMI", True, True, iv, 0)
    aln3 = DummyAlignment("readF_UMI", True, True, iv, 0)
    aln4 = DummyAlignment("readF_UMI", True, True, iv, 0)

    # Case 1: multimapping reads with a primary alignment
    aln1.not_primary_alignment = True
    aln2.not_primary_alignment = True
    aln3.not_primary_alignment = False # primary
    aln4.not_primary_alignment = False # primary
    bundle_multi = [(aln1, aln2), (aln3, aln4)]
    rt1 = set_alignment_category(bundle_multi, bamfile="test", count_primary=True)
    assert rt1.category == ""
    assert rt1.read1_almnt is aln3
    assert rt1.read2_almnt is aln4

    # set all alignments to primary
    aln1.not_primary_alignment = False
    aln2.not_primary_alignment = False
    bundle_multi = [(aln1, aln2), (aln3, aln4)] # all are primary

    # Case 3: testing multiple_primary_action raise
    with pytest.raises(ValueError):
        set_alignment_category(bundle_multi, bamfile="test", count_primary=True, 
                                multiple_primary_action='raise')

    # Case 4: testing multiple_primary_action warn
    rt4 = set_alignment_category(bundle_multi, bamfile="test", count_primary=True, 
                                  multiple_primary_action='warn')
    assert rt4.category == ""

    # Case 5: testing multiple_primary_action skip
    rt5 = set_alignment_category(bundle_multi, bamfile="test", count_primary=True, 
                                  multiple_primary_action='skip')
    assert rt5.category == "_multimapping"

    # Case 6: test nonsense multiple_primary_action
    with pytest.raises(ValueError):
        set_alignment_category(bundle_multi, bamfile="test", count_primary=True, 
                                multiple_primary_action='abc')

def test_read_at_exon_boundary():
    """
    In HTSeq, reads exactly at the boundaries of an exon interval should be
    considered overlapping the exon.
    """
    exon_iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    rt = create_readtrack("read_boundary_UMI", exon_iv)
    
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gfeatures[exon_iv] += "gene1"
    efeatures[exon_iv] += "exon1"
    dummy_eattributes = {"exon1": ["gene1", "GeneName1", "1"]}
    
    rt.find_overlap(gfeatures, efeatures, dummy_eattributes)
    rt.evaluate_overlap(dummy_eattributes)
    assert rt.gene_to_count == "gene1"
    assert rt.exon_to_count == "exon1"

def test_read_wrong_strand():
    """
    A read on the wrong strand should not count as overlapping the feature.
    """
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    gfeatures[iv] += "gene1"
    efeatures[iv] += "exon1"
    dummy_eattributes = {"exon1": ["gene1", "GeneName1", "1"]}
    
    # Create a read aligned on the minus strand.
    aln = DummyAlignment("read_strand", True, False, HTSeq.GenomicInterval("chr1", 100, 200, "-"), 0)
    rt = ReadTrack(aln, aln)
    rt.find_overlap(gfeatures, efeatures, dummy_eattributes)
    rt.evaluate_overlap(dummy_eattributes)
    assert rt.category == "_no_feature"

def test_umi_extraction_cases():
    """
    Test UMI extraction logic from read names.
    """
    test_cases = [
        ("read1_UMI123", "UMI123"),
        ("read1_lane3", "lane3"),       # underscore present, takes suffix
        ("read1", None),                # no underscore
        ("read1_", None),               # empty UMI
        ("_UMI999", "UMI999"),          # leading underscore, still extracts
    ]
    
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    for read_name, expected_umi in test_cases:
        aln = DummyAlignment(read_name, True, False, iv, 0)
        rt = ReadTrack(aln, aln)
        assert rt.umi == expected_umi

def test_readtrack_can_do_overlap():
    """
    Test the can_do_overlap method.
    """
    # Create dummy alignments with both read1 and read2 present.
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    rt = create_readtrack("readA_UMI", iv)
    assert rt.can_do_overlap() is True

    # Now set category and check that it returns False.
    rt.category = '_no_feature'
    assert rt.can_do_overlap() is False

    # Now set one of the alignments to None.
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    aln1 = DummyAlignment("readA_UMI", True, False, iv, 0)
    rt = ReadTrack(aln1, None)
    assert rt.can_do_overlap() is False

def test_readtrack_unique_overlap():
    """
    A read that overlaps a single gene with exon info should be assigned that gene.
    """
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    rt = create_readtrack("read_unique_UMI", iv)
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gfeatures[iv] += "gene1"
    efeatures[iv] += "g1exon1"
    dummy_eattributes = {"g1exon1": ["gene1", "GeneName1", "1"]}
    
    rt.find_overlap(gfeatures, efeatures, dummy_eattributes)
    rt.evaluate_overlap(dummy_eattributes)
    
    assert rt.gene_to_count == "gene1"
    assert rt.exon_to_count == "g1exon1"
    assert rt.category == ""

def test_readtrack_consensus_overlap():
    """
    If both reads from 1 bundle overlap 2 exons it should be assigned that vs another
    gene with only a single exon in that interval.
    """
    ivg1 = HTSeq.GenomicInterval("chr1", 0, 200, "+") # gene1
    ivg1e1 = HTSeq.GenomicInterval("chr1", 0, 50, "+") # gene1exon1
    ivg1e2 = HTSeq.GenomicInterval("chr1", 150, 200, "+") # gene1exon2
    ivg2 = HTSeq.GenomicInterval("chr1", 0, 200, "+") # gene2
    ivg2e1 = HTSeq.GenomicInterval("chr1", 150, 200, "+") # gene2exon1
    ivint = HTSeq.GenomicInterval("chr1", 51, 149, "+") # intergenic

    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gfeatures[ivg1] += "gene1"
    gfeatures[ivg2] += "gene2"
    efeatures[ivg1e1] += "g1exon1"
    efeatures[ivg1e2] += "g1exon2"
    efeatures[ivg1e1] += "g2exon1"
    dummy_eattributes = {
        "g1exon1": ["gene1", "GeneName1", "1"],
        "g1exon2": ["gene1", "GeneName1", "2"],
        "g2exon1": ["gene2", "GeneName2", "1"]
    }

    aln1 = DummyAlignment("read1_UMI", True, False, ivg1e1, 0)
    aln2 = DummyAlignment("read1_UMI", True, False, ivg1e2, 0)
    rt1 = ReadTrack(aln1, aln2)

    rt1.find_overlap(gfeatures, efeatures, dummy_eattributes)
    rt1.evaluate_overlap(dummy_eattributes)

    assert rt1.gene_to_count == "gene1"
    assert rt1.category == ""

    aln3 = DummyAlignment("read2_UMI", True, False, ivg1e1, 0)
    aln4 = DummyAlignment("read2_UMI", True, False, ivint, 0)
    rt2 = ReadTrack(aln3, aln4)

    rt2.find_overlap(gfeatures, efeatures, dummy_eattributes)
    rt2.evaluate_overlap(dummy_eattributes)

    assert rt2.gene_to_count == ""
    assert rt2.category == "_ambiguous"

def test_readtrack_no_feature():
    """
    A read with no overlapping features should be categorized as _no_feature.
    """
    iv = HTSeq.GenomicInterval("chr1", 300, 400, "+")
    rt = create_readtrack("read_nofeat_UMI", iv)
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    
    rt.find_overlap(gfeatures, efeatures, {})
    rt.evaluate_overlap({})
    
    assert rt.category == "_no_feature"
    assert rt.gene_to_count == ""
    assert rt.exon_to_count == ""

def test_readtrack_ambiguous():
    """
    A read overlapping two genes with exons mapping to different genes should be _ambiguous.
    """
    iv = HTSeq.GenomicInterval("chr1", 150, 250, "+")
    rt = create_readtrack("read_ambiguous", iv)
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gfeatures[iv] += "gene1"
    gfeatures[iv] += "gene2"
    efeatures[iv] += "g1exon1"
    efeatures[iv] += "g2exon1"
    dummy_eattributes = {
        "g1exon1": ["gene1", "GeneName1", "1"],
        "g2exon1": ["gene2", "GeneName2", "1"]
    }
    
    rt.find_overlap(gfeatures, efeatures, dummy_eattributes)
    rt.evaluate_overlap(dummy_eattributes)
    
    assert rt.category == "_ambiguous"
    assert rt.gene_to_count == ""
    assert rt.exon_to_count == ""

def test_readtrack_multiple_overlap_single_gene():
    """
    A read overlapping two gene annotations but with exons only from one gene
    should be assigned to that single gene.
    """
    iv = HTSeq.GenomicInterval("chr1", 150, 250, "+")
    rt = create_readtrack("read_multi_single", iv)
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gfeatures[iv] += "gene1"
    gfeatures[iv] += "gene2"
    # Both exons come from gene1:
    efeatures[iv] += "g1exon1"
    efeatures[iv] += "g1exon2"
    dummy_eattributes = {
        "g1exon1": ["gene1", "GeneName1", "1"],
        "g1exon2": ["gene1", "GeneName1", "2"]
    }
    
    rt.find_overlap(gfeatures, efeatures, dummy_eattributes)
    rt.evaluate_overlap(dummy_eattributes)
    
    assert rt.category == ""
    assert rt.gene_to_count == "gene1"
    assert rt.exon_to_count in ["g1exon1", "g1exon2"]

def test_parse_bam_and_count_simple(monkeypatch, tmp_path):
    """
    Test parse_bam_and_count with a simple single-gene scenario.
    """
    # Build a dummy GTF dump with one gene "gene1" and exon information.
    ivgene = HTSeq.GenomicInterval("chr1", 100, 400, "+")
    ivexon1 = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    ivexon2 = HTSeq.GenomicInterval("chr1", 300, 400, "+")
    genes_intervals = [
        ("gene1", "GeneName1", "exon1", ivgene, False), # overarching gene definition
        ("gene1", "GeneName1", "g1exon1", ivexon1, True), # exon1 definition for gene1
        ("gene1", "GeneName1", "g1exon2", ivexon2, True), # exon2 definition for gene1
    ]
    gtf_data = get_dummy_gtf_dump(genes_intervals)
    
    # Bundle 1: paired read with UMI "UMI1" overlapping gene1
    bundle1 = [(DummyAlignment("read1_UMI1", True, False, ivgene, 0), 
                DummyAlignment("read1_UMI1", True, False, ivgene, 0))]

    # Bundle 2: paired read without UMI overlapping exon1
    bundle2 = [(DummyAlignment("read2", True, False, ivexon1, 0), 
                DummyAlignment("read2", True, False, ivexon1, 0))]

    # Bundle 3: paired read without UMI overlapping gene1 intron
    ivintron = HTSeq.GenomicInterval("chr1", 220, 280, "+")
    bundle3 = [(DummyAlignment("read3", True, False, ivintron, 0), 
                DummyAlignment("read3", True, False, ivintron, 0))]

    # Bundle 4: unmapped (one alignment is None), with UMI
    bundle4 = [(DummyAlignment("read4_UMI3", True, False, ivgene, 0), None)]

    # Bundle 5: multimapping, without UMI, no primary alignment
    bundle5 = [(DummyAlignment("read5", True, True, ivexon1, 0), 
                DummyAlignment("read5", True, True, ivexon1, 0)), 
               (DummyAlignment("read5", True, True, ivexon2, 0), 
                DummyAlignment("read5", True, True, ivexon2, 0))]

    # Bundle 6: duplicate of Bundle 1 (same UMI "UMI1")
    bundle6 = [(DummyAlignment("read1_UMI1", True, False, ivgene, 0), 
                DummyAlignment("read1_UMI1", True, False, ivgene, 0))]

    # Bundle 7: multimapping with a primary alignment 
    bundle7 = [(DummyAlignment("read6", True, False, ivexon1, 0), # primary
                DummyAlignment("read6", True, False, ivexon1, 0)), # primary
               (DummyAlignment("read6", True, True, ivexon2, 0), 
                DummyAlignment("read6", True, True, ivexon2, 0))]

    bundles = [bundle1, bundle2, bundle3, bundle4, bundle5, bundle6, bundle7]
    monkeypatch.setattr(HTSeq, "BAM_Reader", lambda bamfile: bamfile)
    monkeypatch.setattr(HTSeq, "pair_SAM_alignments", dummy_pair_SAM_alignments_factory(bundles))
    
    counts, _ = parse_bam_and_count("dummy.bam", gtf_data, 
                                    count_primary=True,
                                    cols_to_use=['UE', 'RE', 'UI', 'RI', 'D'])
    
    gene1_counts = counts["gene1"]
    assert gene1_counts["UE"] == 1 # from bundle1
    assert gene1_counts["D"] == 1 # from bundle6
    assert gene1_counts["RE"] == 2 # from bundle2, bundle7
    assert gene1_counts["RI"] == 1 # from bundle3

    unmapped = counts["_unmapped"] # from bundle4
    assert unmapped["UE"] == 1

    multimapping = counts["_multimapping"] # from bundle5
    assert multimapping["RE"] == 1

def test_parse_bam_and_count_complex(monkeypatch, tmp_path):
    """
    Test parse_bam_and_count in a multi-gene scenario.
    """
    # Build a dummy GTF dump with two genes.
    ivgene1 = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    ivgene1exon1 = HTSeq.GenomicInterval("chr1", 100, 150, "+")
    ivgene1exon2 = HTSeq.GenomicInterval("chr1", 160, 180, "+")
    ivgene2 = HTSeq.GenomicInterval("chr1", 250, 350, "+")
    genes_intervals = [
        ("gene1", "GeneName1", "g1exon1", ivgene1, False), # gene1 
        ("gene1", "GeneName1", "g1exon1", ivgene1exon1, True), # define gene1 exon1
        ("gene1", "GeneName1", "g1exon2", ivgene1exon2, True), # define gene1 exon2, leaving some intron
        ("gene2", "GeneName2", "g2exon1", ivgene2, False), # gene2 without exon, ie assume intron
    ]
    gtf_data = get_dummy_gtf_dump(genes_intervals)
    
    # Bundle 1: proper paired read with UMI "UMI1" overlapping gene1
    bundle1 = [(DummyAlignment("read1_UMI1", True, False, ivgene1, 0), 
                DummyAlignment("read1_UMI1", True, False, ivgene1, 0))]

    # Bundle 2: proper paired read with UMI "UMI2" overlapping gene2 which is intronic
    bundle2 = [(DummyAlignment("read2_UMI2", True, False, ivgene2, 0), 
                DummyAlignment("read2_UMI2", True, False, ivgene2, 0))]

    # Bundle 3: ambiguous read overlaps gene1 and gene2, but we count gene1 because has exon
    amb_iv = HTSeq.GenomicInterval("chr1", 170, 260, "+")
    bundle3 = [(DummyAlignment("readAmb", True, False, amb_iv, 0), 
                DummyAlignment("readAmb", True, False, amb_iv, 0))]

    # Bundle 4: duplicate of Bundle 1 for gene1
    bundle4 = [(DummyAlignment("read1_UMI1", True, False, ivgene1, 0), 
                DummyAlignment("read1_UMI1", True, False, ivgene1, 0))]

    # Bundle 5: read overlaps gene1 and gene2, but ambiguous because no exons
    amb_iv = HTSeq.GenomicInterval("chr1", 190, 260, "+")
    bundle3 = [(DummyAlignment("readAmb", True, False, amb_iv, 0), 
                DummyAlignment("readAmb", True, False, amb_iv, 0))]
    
    bundles = [bundle1, bundle2, bundle3, bundle4]
    monkeypatch.setattr(HTSeq, "BAM_Reader", lambda bamfile: bamfile)
    monkeypatch.setattr(HTSeq, "pair_SAM_alignments", dummy_pair_SAM_alignments_factory(bundles))
    
    counts, _ = parse_bam_and_count("dummy.bam", gtf_data, cols_to_use=['UE', 'RE', 'UI', 'RI', 'D'])
    
    gene1_counts = counts["gene1"]
    gene2_counts = counts["gene2"]
    
    assert gene1_counts["UE"] == 1 # from bundle1
    assert gene1_counts["D"] == 1 # from bundle4

    assert gene2_counts["UI"] == 1 # from bundle2 (gene is intronic)
    assert counts["_ambiguous"]["RE"] == 1 # from bundle3

def test_parse_bam_and_count_umi_deduplication(monkeypatch, tmp_path):
    """
    Test UMI deduplication scenarios in parse_bam_and_count.
    """
    # Build a dummy GTF dump with two genes.
    ivgene1 = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    ivgene2 = HTSeq.GenomicInterval("chr1", 250, 350, "+")
    genes_intervals = [
        ("gene1", "GeneName1", "g1exon1", ivgene1, True), # gene1 
        ("gene2", "GeneName2", "g2exon1", ivgene2, True), # gene2 
    ]
    gtf_data = get_dummy_gtf_dump(genes_intervals)
    
    # Bundle 1: proper paired read with UMI to gene1
    bundle1 = [(DummyAlignment("read1_UMI1", True, False, ivgene1, 0), 
                DummyAlignment("read1_UMI1", True, False, ivgene1, 0))]

    # Bundle 2: proper paired read with same UMI as Bundle1, for gene2
    bundle2 = [(DummyAlignment("read2_UMI1", True, False, ivgene2, 0), 
                DummyAlignment("read2_UMI1", True, False, ivgene2, 0))]

    # Bundle 3: proper paired read with same UMI as Bundle1, but different position
    ivgene1_shift = HTSeq.GenomicInterval("chr1", 120, 220, "+")
    bundle3 = [(DummyAlignment("read1_UMI1", True, False, ivgene1_shift, 0), 
                DummyAlignment("read1_UMI1", True, False, ivgene1_shift, 0))]

    # Bundle 4: proper paired read with same base name as Bundle1 but different UMI
    bundle4 = [(DummyAlignment("read1_UMI2", True, False, ivgene1, 0), 
                DummyAlignment("read1_UMI2", True, False, ivgene1, 0))]
    
    bundles = [bundle1, bundle2, bundle3, bundle4]
    monkeypatch.setattr(HTSeq, "BAM_Reader", lambda bamfile: bamfile)
    monkeypatch.setattr(HTSeq, "pair_SAM_alignments", dummy_pair_SAM_alignments_factory(bundles))
    
    counts, _ = parse_bam_and_count("dummy.bam", gtf_data, cols_to_use=['UE', 'RE', 'UI', 'RI', 'D'])

    gene1_counts = counts["gene1"]
    gene2_counts = counts["gene2"]
    
    assert gene1_counts["UE"] == 2 # from bundle1, bundle4
    assert gene2_counts["UE"] == 1 # from bundle2
    assert gene1_counts["D"] == 1 # from bundle3
