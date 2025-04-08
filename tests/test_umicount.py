import pickle
import pytest
import HTSeq

from umicount.umicount import (
    parse_gtf,
    ReadTrack,
    extract_first_alignment,
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
    def __init__(self, name, aligned=True, iv=None):
        DummyRead = type("DummyRead", (), {})  # create an empty dummy class
        self.read = DummyRead()
        self.read.name = name
        self.aligned = aligned
        self.iv = iv

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

def create_readtrack(name, iv, aligned=True):
    """
    Helper to create a ReadTrack from two dummy alignments.
    """
    aln1 = DummyAlignment(name, aligned, iv)
    aln2 = DummyAlignment(name, aligned, iv)
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

def test_extract_first_alignment():
    """
    Test extract_first_alignment for different bundle cases.
    """
    iv = HTSeq.GenomicInterval("chr1", 100, 200, "+")
    # Case 1: proper pair.
    aln1 = DummyAlignment("readC_UMI", aligned=True, iv=iv)
    aln2 = DummyAlignment("readC_UMI", aligned=True, iv=iv)
    bundle = [(aln1, aln2)]
    rt = extract_first_alignment(bundle)
    assert rt.category == ""
    assert rt.read1_almnt.read.name == "readC_UMI"
    assert rt.read2_almnt is not None

    # Case 2: one alignment is None (simulate unmapped).
    bundle_unmapped = [(aln1, None)]
    rt2 = extract_first_alignment(bundle_unmapped)
    assert rt2.category == "_unmapped"
    assert rt2.read2_almnt is None

    # Case 3: both alignments exist but one is not aligned.
    aln3 = DummyAlignment("readD_UMI", aligned=False, iv=iv)
    bundle_unmapped2 = [(aln3, aln2)]
    rt3 = extract_first_alignment(bundle_unmapped2)
    assert rt3.category == "_unmapped"
    assert rt3.read2_almnt is None

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
    
    rt.find_overlap(gfeatures, efeatures)
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
    aln = DummyAlignment("read_strand", True, HTSeq.GenomicInterval("chr1", 100, 200, "-"))
    rt = ReadTrack(aln, aln)
    rt.find_overlap(gfeatures, efeatures)
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
        aln = DummyAlignment(read_name, True, iv)
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
    aln1 = DummyAlignment("readA_UMI", aligned=True, iv=iv)
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
    efeatures[iv] += "exon1"
    dummy_eattributes = {"exon1": ["gene1", "GeneName1", "1"]}
    
    rt.find_overlap(gfeatures, efeatures)
    rt.evaluate_overlap(dummy_eattributes)
    
    assert rt.gene_to_count == "gene1"
    assert rt.exon_to_count == "exon1"
    assert rt.category == ""

def test_readtrack_no_feature():
    """
    A read with no overlapping features should be categorized as _no_feature.
    """
    iv = HTSeq.GenomicInterval("chr1", 300, 400, "+")
    rt = create_readtrack("read_nofeat_UMI", iv)
    gfeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    efeatures = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    
    rt.find_overlap(gfeatures, efeatures)
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
    efeatures[iv] += "exon1"
    efeatures[iv] += "exon2"
    dummy_eattributes = {
        "exon1": ["gene1", "GeneName1", "1"],
        "exon2": ["gene2", "GeneName2", "1"]
    }
    
    rt.find_overlap(gfeatures, efeatures)
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
    efeatures[iv] += "exon1"
    efeatures[iv] += "exon2"
    dummy_eattributes = {
        "exon1": ["gene1", "GeneName1", "1"],
        "exon2": ["gene1", "GeneName1", "2"]
    }
    
    rt.find_overlap(gfeatures, efeatures)
    rt.evaluate_overlap(dummy_eattributes)
    
    assert rt.category == ""
    assert rt.gene_to_count == "gene1"
    assert rt.exon_to_count == "gene1"

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
        ("gene1", "GeneName1", "exon1", ivexon1, True), # exon1 definition for gene1
        ("gene1", "GeneName1", "exon2", ivexon2, True), # exon2 definition for gene1
    ]
    gtf_data = get_dummy_gtf_dump(genes_intervals)
    
    # Bundle 1: paired read with UMI "UMI1" overlapping gene1
    bundle1 = [(DummyAlignment("read1_UMI1", True, ivgene), 
                DummyAlignment("read1_UMI1", True, ivgene))]

    # Bundle 2: paired read without UMI overlapping exon1
    bundle2 = [(DummyAlignment("read2", True, ivexon1), 
                DummyAlignment("read2", True, ivexon1))]

    # Bundle 3: paired read without UMI overlapping gene1 intron
    ivintron = HTSeq.GenomicInterval("chr1", 220, 280, "+")
    bundle3 = [(DummyAlignment("read3", True, ivintron), 
                DummyAlignment("read3", True, ivintron))]

    # Bundle 4: unmapped (one alignment is None), with UMI
    bundle4 = [(DummyAlignment("read4_UMI3", True, ivgene), None)]

    # Bundle 5: multimapping, without UMI
    bundle5 = [(DummyAlignment("read5", True, ivexon1), 
                DummyAlignment("read5", True, ivexon1)), 
               (DummyAlignment("read5", True, ivexon2), 
                DummyAlignment("read5", True, ivexon2))]

    # Bundle 6: duplicate of Bundle 1 (same UMI "UMI1")
    bundle6 = [(DummyAlignment("read1_UMI1", True, ivgene), 
                DummyAlignment("read1_UMI1", True, ivgene))]
    
    bundles = [bundle1, bundle2, bundle3, bundle4, bundle5, bundle6]
    monkeypatch.setattr(HTSeq, "BAM_Reader", lambda bamfile: bamfile)
    monkeypatch.setattr(HTSeq, "pair_SAM_alignments", dummy_pair_SAM_alignments_factory(bundles))
    
    counts, _ = parse_bam_and_count("dummy.bam", gtf_data, cols_to_use=['UE', 'RE', 'UI', 'RI', 'D'])
    
    gene1_counts = counts["gene1"]
    assert gene1_counts["UE"] == 1 # from bundle1
    assert gene1_counts["D"] == 1 # from bundle6
    assert gene1_counts["RE"] == 1 # from bundle2
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
        ("gene1", "GeneName1", "exon1", ivgene1, False), # gene1 
        ("gene1", "GeneName1", "exon1", ivgene1exon1, True), # define gene1 exon1
        ("gene1", "GeneName1", "exon2", ivgene1exon2, True), # define gene1 exon2, leaving some intron
        ("gene2", "GeneName2", "exon1", ivgene2, False), # gene2 without exon, ie assume intron
    ]
    gtf_data = get_dummy_gtf_dump(genes_intervals)
    
    # Bundle 1: proper paired read with UMI "UMI1" overlapping gene1
    bundle1 = [(DummyAlignment("read1_UMI1", True, ivgene1), 
                DummyAlignment("read1_UMI1", True, ivgene1))]

    # Bundle 2: proper paired read with UMI "UMI2" overlapping gene2 which is intronic
    bundle2 = [(DummyAlignment("read2_UMI2", True, ivgene2), 
                DummyAlignment("read2_UMI2", True, ivgene2))]

    # Bundle 3: ambiguous read overlaps gene1 and gene2, but we count gene1 because has exon
    amb_iv = HTSeq.GenomicInterval("chr1", 170, 260, "+")
    bundle3 = [(DummyAlignment("readAmb", True, amb_iv), 
                DummyAlignment("readAmb", True, amb_iv))]

    # Bundle 4: duplicate of Bundle 1 for gene1
    bundle4 = [(DummyAlignment("read1_UMI1", True, ivgene1), 
                DummyAlignment("read1_UMI1", True, ivgene1))]

    # Bundle 5: read overlaps gene1 and gene2, but ambiguous because no exons
    amb_iv = HTSeq.GenomicInterval("chr1", 190, 260, "+")
    bundle3 = [(DummyAlignment("readAmb", True, amb_iv), 
                DummyAlignment("readAmb", True, amb_iv))]
    
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
        ("gene1", "GeneName1", "exon1", ivgene1, True), # gene1 
        ("gene2", "GeneName2", "exon1", ivgene2, True), # gene2 
    ]
    gtf_data = get_dummy_gtf_dump(genes_intervals)
    
    # Bundle 1: proper paired read with UMI to gene1
    bundle1 = [(DummyAlignment("read1_UMI1", True, ivgene1), 
                DummyAlignment("read1_UMI1", True, ivgene1))]

    # Bundle 2: proper paired read with same UMI as Bundle1, for gene2
    bundle2 = [(DummyAlignment("read2_UMI1", True, ivgene2), 
                DummyAlignment("read2_UMI1", True, ivgene2))]

    # Bundle 3: proper paired read with same UMI as Bundle1, but different position
    ivgene1_shift = HTSeq.GenomicInterval("chr1", 120, 220, "+")
    bundle3 = [(DummyAlignment("read1_UMI1", True, ivgene1_shift), 
                DummyAlignment("read1_UMI1", True, ivgene1_shift))]

    # Bundle 4: proper paired read with same base name as Bundle1 but different UMI
    bundle4 = [(DummyAlignment("read1_UMI2", True, ivgene1), 
                DummyAlignment("read1_UMI2", True, ivgene1))]
    
    bundles = [bundle1, bundle2, bundle3, bundle4]
    monkeypatch.setattr(HTSeq, "BAM_Reader", lambda bamfile: bamfile)
    monkeypatch.setattr(HTSeq, "pair_SAM_alignments", dummy_pair_SAM_alignments_factory(bundles))
    
    counts, _ = parse_bam_and_count("dummy.bam", gtf_data, cols_to_use=['UE', 'RE', 'UI', 'RI', 'D'])

    gene1_counts = counts["gene1"]
    gene2_counts = counts["gene2"]
    
    assert gene1_counts["UE"] == 2 # from bundle1, bundle4
    assert gene2_counts["UE"] == 1 # from bundle2
    assert gene1_counts["D"] == 1 # from bundle3
