import pytest
import gzip
import tempfile
import os
import HTSeq as htseq
import re
import regex

from umicount.umiextract import (
    process_fastq,
    process_entry,
    slice_SequenceWithQualities,
    get_fastq_str
)

# Helper function to write a FASTQ file
def write_fastq(filename, reads):
    with gzip.open(filename, "wt") as f:
        for read in reads:
            f.write(f"@{read[0]}\n{read[1]}\n+\n{read[2]}\n")

############################
# Tests for process_fastq
############################

def test_process_fastq_basic():
    """Test processing with a simple valid case."""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("seq2", "ATTGCGCAATGTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("seq3", "ATGGCCGAACGAGATACGTAGTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
    ]
    r2_reads = [
        ("seq1", "GGGACGTAAAGGG", "IIIIIIIIIIIII"),
        ("seq2", "GGGTTTTAAA", "IIIIIIIIII"),
        ("seq3", "ACGATACGATCGA", "IIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r2_path = os.path.join(tmpdir, "r2.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")
        r2_out_path = os.path.join(tmpdir, "r2_out.fastq.gz")

        write_fastq(r1_path, r1_reads)
        write_fastq(r2_path, r2_reads)

        process_fastq((r1_path, r2_path), (r1_out_path, r2_out_path), 
                      umi_len=4, only_umi=False, fuzzy_umi_params=None)

        # Read the output
        with htseq.FastqReader(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 3
        assert out_reads[0].name.endswith("_ACGT")
        assert out_reads[1].name.endswith("_TTTT")
        assert out_reads[2].name == "seq3"

def test_process_fastq_umi_len():
    """Test processing for variable UMI length."""
    r1_reads = [
        ("seq1", "ATTGCGCAATGATTCGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIIIII"),
        ("seq2", "ATTGCGCAATGTCTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        process_fastq((r1_path, None), (r1_out_path, None), 
                      umi_len=6, only_umi=False, fuzzy_umi_params=None)

        with htseq.FastqReader(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 2
        assert out_reads[0].name.endswith("_ATTCGT")
        assert out_reads[1].name.endswith("_TCTTTT")

def test_process_fastq_empty_input():
    """Test processing when input is empty."""
    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r2_path = os.path.join(tmpdir, "r2.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")
        r2_out_path = os.path.join(tmpdir, "r2_out.fastq.gz")

        write_fastq(r1_path, [])
        write_fastq(r2_path, [])

        process_fastq((r1_path, r2_path), (r1_out_path, r2_out_path), 
                      umi_len=4, only_umi=False, fuzzy_umi_params=None)

        with htseq.FastqReader(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 0

def test_process_fastq_malformed_fastq():
    """Test that mismatched sequence and quality string lengths raise error"""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIII"),
        ("seq2", "ATTGCGCAATGTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        with pytest.raises(ValueError):  # parser raises error
            process_fastq((r1_path, None), (r1_out_path, None), 
                          umi_len=4, only_umi=False, fuzzy_umi_params=None)

def test_process_fastq_no_matching_umi():
    """Test processing when no UMI pattern is found."""
    r1_reads = [
        ("seq1", "GGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIII"),
        ("seq2", "CCCCCCCCCCCCCCCCCC", "IIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        process_fastq((r1_path, None), (r1_out_path, None), 
                      umi_len=4, only_umi=True, fuzzy_umi_params=None)

        with htseq.FastqReader(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 0  # No UMI, so no reads should be written

def test_process_fastq_duplicate_read_name():
    """Test that a duplicate read name raises an error."""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("seq1", "ATTGCGCAATGTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        with pytest.raises(SystemExit):  # Should exit due to duplicate read name
            process_fastq((r1_path, None), (r1_out_path, None), 
                          umi_len=4, only_umi=False, fuzzy_umi_params=None)

def test_process_fastq_mismatched_read_names():
    """Test that mismatched R1 and R2 readnames raise an error."""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGGACGTC", "IIIIIIIIIIIIIIIIIIIIIII"),
        ("seq2", "ATTGCGCAATGTTTTGGGACGTG", "IIIIIIIIIIIIIIIIIIIIIII"),
    ]
    r2_reads = [
        ("seq1", "GGGACGTAAAGGG", "IIIIIIIIIIIII"),
        ("seq3", "GGGTTTTAAA", "IIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r2_path = os.path.join(tmpdir, "r2.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")
        r2_out_path = os.path.join(tmpdir, "r2_out.fastq.gz")

        write_fastq(r1_path, r1_reads)
        write_fastq(r2_path, r2_reads)

        with pytest.raises(SystemExit):  # Should exit due to mismatched read names
            process_fastq((r1_path, r2_path), (r1_out_path, r2_out_path), 
                          umi_len=4, only_umi=False, fuzzy_umi_params=None)

def test_process_fastq_empty_sequences():
    """Test that empty sequences are properly removed."""
    r1_reads = [
        ("seq1", "ATTGCGCAATGACGTGGG", "IIIIIIIIIIIIIIIIII"),
        ("seq2", "ATTGCGCAATGTTTTGGG", "IIIIIIIIIIIIIIIIII"),
    ]
    r2_reads = [
        ("seq1", "GGGACGTAAAGGG", "IIIIIIIIIIIII"),
        ("seq2", "GGGTTTTAAA", "IIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r2_path = os.path.join(tmpdir, "r2.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")
        r2_out_path = os.path.join(tmpdir, "r2_out.fastq.gz")

        write_fastq(r1_path, r1_reads)
        write_fastq(r2_path, r2_reads)

        process_fastq((r1_path, r2_path), (r1_out_path, r2_out_path), 
                      umi_len=4, only_umi=False, fuzzy_umi_params=None)

        with htseq.FastqReader(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 0  # Both sequences should be removed

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

############################
# Tests for utility functions
############################

def test_slice_sequence_with_qualities():
    """Test that slicing a SequenceWithQualities works as expected."""
    read = htseq.SequenceWithQualities(b"ACGTACGT", "read1", b"IIIIIIII")
    sliced = slice_SequenceWithQualities(read, 4)
    assert sliced.seq == b"ACGT"
    assert sliced.qualstr == b"IIII"
    assert sliced.name == "read1"

def test_get_fastq_str():
    """Test that the FASTQ string is correctly generated."""
    read = htseq.SequenceWithQualities(b"ACGT", "read1", b"IIII")
    fastq = get_fastq_str(read)
    expected = b"@read1\nACGT\n+\nIIII\n"
    assert fastq == expected

if __name__ == "__main__":
    pytest.main()

