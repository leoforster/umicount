import pytest
import gzip
import tempfile
import os
import HTSeq as htseq
import re
import regex

from umicount.umiextract import (slice_SequenceWithQualities, get_fastq_str)

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
