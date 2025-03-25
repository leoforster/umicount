import pytest
import gzip
import tempfile
import os
import HTSeq as htseq
import re
import regex

from umicount.umiextract import process_entry

############################
# Tests for process_entry
############################

def test_process_entry_default_matching():
    """Test that exact UMI hits are extracted properly"""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGACGTC", "IIIIIIIIIIIIIIIIIIIIIII"),
        ("seq2", "ATTGCGCAATGTTTTGGGACGTC", "IIIIIIIIIIIIIIIIIIIIIII"),
        ("seq3", "ATTGCGCTATGTTTTGGGACGTC", "IIIIIIIIIIIIIIIIIIIIIII"), # 1 mismatch
    ]

    umi_len = 4
    anchor_seq = 'ATTGCGCAATG'
    trailing_seq = 'GGG'

    pattern = re.compile(f"({anchor_seq})[NGCAT]{{{umi_len}}}({trailing_seq})")

    entry_names = []
    for r in r1_reads:
        entry = htseq.SequenceWithQualities(r[1].encode(), r[0], r[2].encode())
        entry_processed, umi = process_entry(entry, pattern, umi_len, only_umi=False,
                                             search_region=-1, min_remaining_seqlen=-1,
                                             anchor_seq=anchor_seq, trailing_seq=trailing_seq,
                                             fuzzy_umi_params=None)
        entry_names.append(entry_processed)

    assert entry_names[0].name == "seq1_ACGT"
    assert entry_names[1].name == "seq2_TTTT"
    assert entry_names[2].name == "seq3"

def test_process_entry_fuzzy_matching_basic():
    """Test that fuzzy 1-mismatch UMI hits are extracted properly"""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGACGTC", "IIIIIIIIIIIIIIIIIIIIIII"), # 0 mismatch
        ("seq2", "ATTGCGCTATGACACGGGACGTC", "IIIIIIIIIIIIIIIIIIIIIII"), # 1 mismatch
        ("seq3", "ATTGCGCTTTGTATAGGGACGTC", "IIIIIIIIIIIIIIIIIIIIIII"), # 2 mismatches
    ]

    umi_len = 4
    anchor_seq = 'ATTGCGCAATG'
    trailing_seq = 'GGG'

    anchor_max_mismatch = 1
    anchor_max_indel = 0
    min_trailing_G = 3

    anchor_fuzzy = rf"({anchor_seq}){{e<={anchor_max_mismatch + anchor_max_indel}}}"
    umi_capture = rf"([ACGTN]{{{umi_len}}})"
    pattern = regex.compile(anchor_fuzzy + umi_capture, flags=regex.BESTMATCH)

    entry_names = []
    for r in r1_reads:
        entry = htseq.SequenceWithQualities(r[1].encode(), r[0], r[2].encode())
        entry_processed, umi = process_entry(entry, pattern, umi_len, only_umi=False,
                                             search_region=-1, min_remaining_seqlen=-1,
                                             anchor_seq=anchor_seq, trailing_seq=trailing_seq,
                                             fuzzy_umi_params={'anchor_max_mismatch':anchor_max_mismatch,
                                                               'anchor_max_indel':anchor_max_indel,
                                                               'min_trailing_G':min_trailing_G})
        entry_names.append(entry_processed)

    assert entry_names[0].name == "seq1_ACGT"
    assert entry_names[1].name == "seq2_ACAC"
    assert entry_names[2].name == "seq3"

def test_process_entry_fuzzy_matching():
    """Test that different fuzzy UMI hits are extracted properly"""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGACGTC"), # 0 mismatch
        ("seq2", "ATTGCGCTATGACACGGGACGTC"), # 1 mismatch
        ("seq3", "ATTGCGCTTTGTATAGGGACGTC"), # 2 mismatches
        ("seq4", "ATTGCGCTCTTGTATAGGGACGTC"), # 2 mismatches, 1 indel
        ("seq5", "ATTGCGCTCCTTGTATAGGGACGTC"), # 2 mismatches, 2 indel
        ("seq6", "ATTGCGCATTGTATAGGGACGTC"), # 1 mismatch, 1 indel
        ("seq7", "ATTGCGCTTGTCTCGGGACGTC"), # 1 indel
        ("seq8", "ATTGCGCATTGTATACCGACGTC"), # 1 mismatch, 1 Gs
        ("seq9", "ATTGCGCATTGTGTGCGGACGTC"), # 1 mismatch, 2 Gs
    ]

    umi_len = 4
    anchor_seq = 'ATTGCGCAATG'
    trailing_seq = 'GGG'

    anchor_max_mismatch = 1
    anchor_max_indel = 1
    min_trailing_G = 2

    anchor_fuzzy = rf"({anchor_seq}){{e<={anchor_max_mismatch + anchor_max_indel}}}"
    umi_capture = rf"([ACGTN]{{{umi_len}}})"
    pattern = regex.compile(anchor_fuzzy + umi_capture, flags=regex.BESTMATCH)

    entry_names = []
    for r in r1_reads:
        entry = htseq.SequenceWithQualities(r[1].encode(), r[0], ('I'*len(r[1])).encode())
        entry_processed, umi = process_entry(entry, pattern, umi_len, only_umi=False,
                                             search_region=-1, min_remaining_seqlen=-1,
                                             anchor_seq=anchor_seq, trailing_seq=trailing_seq,
                                             fuzzy_umi_params={'anchor_max_mismatch':anchor_max_mismatch,
                                                               'anchor_max_indel':anchor_max_indel,
                                                               'min_trailing_G':min_trailing_G})
        entry_names.append(entry_processed)

    assert entry_names[0].name == "seq1_ACGT"
    assert entry_names[1].name == "seq2_ACAC"
    assert entry_names[2].name == "seq3"
    assert entry_names[3].name == "seq4"
    assert entry_names[4].name == "seq5"
    assert entry_names[5].name == "seq6_TATA"
    assert entry_names[6].name == "seq7_TCTC"
    assert entry_names[7].name == "seq8"
    assert entry_names[8].name == "seq9_TGTG"

def test_process_entry_fuzzy_matching_onlyumi():
    """Test that different fuzzy UMI hits are extracted properly"""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGACGTC"), # 0 mismatch
        ("seq2", "ATTGCGCTATGACACGGGACGTC"), # 1 mismatch
        ("seq3", "ATTGCGCTTTGTATAGGGACGTC"), # 2 mismatches
        ("seq4", "ATTGCGCTCTTGTATAGGGACGTC"), # 2 mismatches, 1 indel
        ("seq5", "ATTGCGCTCCTTGTATAGGGACGTC"), # 2 mismatches, 2 indel
        ("seq6", "ATTGCGCATTGTATAGGGACGTC"), # 1 mismatch, 1 indel
        ("seq7", "ATTGCGCTTGTCTCGGGACGTC"), # 1 indel
        ("seq8", "ATTGCGCATTGTATACCGACGTC"), # 1 mismatch, 1 Gs
        ("seq9", "ATTGCGCATTGTGTGCGGACGTC"), # 1 mismatch, 2 Gs
    ]

    umi_len = 4
    anchor_seq = 'ATTGCGCAATG'
    trailing_seq = 'GGG'

    anchor_max_mismatch = 1
    anchor_max_indel = 1
    min_trailing_G = 2

    anchor_fuzzy = rf"({anchor_seq}){{e<={anchor_max_mismatch + anchor_max_indel}}}"
    umi_capture = rf"([ACGTN]{{{umi_len}}})"
    pattern = regex.compile(anchor_fuzzy + umi_capture, flags=regex.BESTMATCH)

    entry_names = []
    for r in r1_reads:
        entry = htseq.SequenceWithQualities(r[1].encode(), r[0], ('I'*len(r[1])).encode())
        entry_processed, umi = process_entry(entry, pattern, umi_len, only_umi=True,
                                             search_region=-1, min_remaining_seqlen=-1,
                                             anchor_seq=anchor_seq, trailing_seq=trailing_seq,
                                             fuzzy_umi_params={'anchor_max_mismatch':anchor_max_mismatch,
                                                               'anchor_max_indel':anchor_max_indel,
                                                               'min_trailing_G':min_trailing_G})
        entry_names.append(entry_processed)

    assert entry_names[0].name == "seq1_ACGT"
    assert entry_names[1].name == "seq2_ACAC"
    assert entry_names[2] is None
    assert entry_names[3] is None
    assert entry_names[4] is None
    assert entry_names[5].name == "seq6_TATA"
    assert entry_names[6].name == "seq7_TCTC"
    assert entry_names[7] is None
    assert entry_names[8].name == "seq9_TGTG"

def test_process_entry_searchregion_length():
    """Test that long search_region with short reads errors"""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGACGTC"),
        ("seq2", "ATTGCGCAATGACGTGGGACGTC"),
    ]

    search_region = 10

    umi_len = 4
    anchor_seq = 'ATTGCGCAATG'
    trailing_seq = 'GGG'

    pattern = re.compile(f"({anchor_seq})[NGCAT]{{{umi_len}}}({trailing_seq})")

    with pytest.raises(SystemExit):
        for r in r1_reads:
            entry = htseq.SequenceWithQualities(r[1].encode(), r[0], ('I'*len(r[1])).encode())
            entry_processed, umi = process_entry(entry, pattern, umi_len, only_umi=False,
                                                 search_region=search_region, min_remaining_seqlen=-1,
                                                 anchor_seq=anchor_seq, trailing_seq=trailing_seq,
                                                 fuzzy_umi_params=None)

def test_process_entry_minremainingseqlen_length():
    """Test that min_remaining_seqlen excludes shorter reads"""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGACGTC"), # 5 remaining seqlen
        ("seq2", "ATTGCGCAATGACGTGGGACGTCACGTC"), # 10 remaining seqlen
    ]

    min_remaining_seqlen = 5

    umi_len = 4
    anchor_seq = 'ATTGCGCAATG'
    trailing_seq = 'GGG'

    pattern = re.compile(f"({anchor_seq})[NGCAT]{{{umi_len}}}({trailing_seq})")

    entry_names = []
    for r in r1_reads:
        entry = htseq.SequenceWithQualities(r[1].encode(), r[0], ('I'*len(r[1])).encode())
        entry_processed, umi = process_entry(entry, pattern, umi_len, only_umi=False,
                                             search_region=-1, min_remaining_seqlen=min_remaining_seqlen,
                                             anchor_seq=anchor_seq, trailing_seq=trailing_seq,
                                             fuzzy_umi_params=None)
        entry_names.append(entry_processed)

    assert entry_names[0] is None
    assert entry_names[1].name == "seq2_ACGT"
