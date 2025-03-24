import pytest
import gzip
import tempfile
import os
import HTSeq as htseq

from umicount.umiextract import process_fastq

# Helper function to write a FASTQ file
def write_fastq(filename, reads):
    with gzip.open(filename, "wt") as f:
        for read in reads:
            f.write(f"@{read[0]}\n{read[1]}\n+\n{read[2]}\n")

# Test cases
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
                      umi_len=4, only_umi=False)

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
                      umi_len=6, only_umi=False)

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
                      umi_len=4, only_umi=False)

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
                          umi_len=4, only_umi=False)

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
                      umi_len=4, only_umi=True)

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
                          umi_len=4, only_umi=False)

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
                          umi_len=4, only_umi=False)

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
                      umi_len=4, only_umi=False)

        with htseq.FastqReader(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 0  # Both sequences should be removed

if __name__ == "__main__":
    pytest.main()

