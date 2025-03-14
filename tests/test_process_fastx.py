import pytest
import gzip
import tempfile
import os
import pysam

from umicount.umiextract import process_fastx

# Helper function to write a FASTQ file
def write_fastq(filename, reads):
    with gzip.open(filename, "wt") as f:
        for read in reads:
            f.write(f"@{read[0]}\n{read[1]}\n+\n{read[2]}\n")

# Test cases
def test_process_fastx_basic():
    """Test processing with a simple valid case."""
    r1_reads = [
        ("read1", "ATTGCGCAATGACGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read2", "ATTGCGCAATGTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read3", "ATGGCCGAACGAGATACGTAGTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
    ]
    r2_reads = [
        ("read1", "GGGACGTAAAGGG", "IIIIIIIIIIIII"),
        ("read2", "GGGTTTTAAA", "IIIIIIIIII"),
        ("read3", "ACGATACGATCGA", "IIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r2_path = os.path.join(tmpdir, "r2.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")
        r2_out_path = os.path.join(tmpdir, "r2_out.fastq.gz")

        write_fastq(r1_path, r1_reads)
        write_fastq(r2_path, r2_reads)

        process_fastx((r1_path, r2_path), (r1_out_path, r2_out_path), umilen=4, only_umi=False)

        # Read the output
        with pysam.FastxFile(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 3
        assert out_reads[0].name.endswith("_ACGT")
        assert out_reads[1].name.endswith("_TTTT")
        assert out_reads[2].name == "read3"

def test_process_fastx_umilen():
    """Test processing for variable UMI length."""
    r1_reads = [
        ("read1", "ATTGCGCAATGATTCGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIIIII"),
        ("read2", "ATTGCGCAATGTCTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        process_fastx((r1_path, None), (r1_out_path, None), umilen=6, only_umi=False)

        # Read the output
        with pysam.FastxFile(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 2
        assert out_reads[0].name.endswith("_ATTCGT")
        assert out_reads[1].name.endswith("_TCTTTT")

def test_process_fastx_empty_input():
    """Test processing when input is empty."""
    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r2_path = os.path.join(tmpdir, "r2.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")
        r2_out_path = os.path.join(tmpdir, "r2_out.fastq.gz")

        write_fastq(r1_path, [])
        write_fastq(r2_path, [])

        process_fastx((r1_path, r2_path), (r1_out_path, r2_out_path), umilen=4, only_umi=False)

        with pysam.FastxFile(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 0

def test_process_fastx_malformed_fastq():
    """Test that a duplicate read name raises an error."""
    r1_reads = [
        ("read1", "ATTGCGCAATGACGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIII"),
        ("read1", "ATTGCGCAATGTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        with pytest.raises(ValueError):  # pysam raises error
            process_fastx((r1_path, None), (r1_out_path, None), umilen=4, only_umi=False)

def test_process_fastx_no_matching_umi():
    """Test processing when no UMI pattern is found."""
    r1_reads = [
        ("read1", "GGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIII"),
        ("read2", "CCCCCCCCCCCCCCCCCC", "IIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        process_fastx((r1_path, None), (r1_out_path, None), umilen=4, only_umi=True)

        with pysam.FastxFile(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 0  # No UMI, so no reads should be written

def test_process_fastx_duplicate_read_name():
    """Test that a duplicate read name raises an error."""
    r1_reads = [
        ("read1", "ATTGCGCAATGACGTGGGTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read1", "ATTGCGCAATGTTTTGGGCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")

        write_fastq(r1_path, r1_reads)

        with pytest.raises(SystemExit):  # Should exit due to duplicate read name
            process_fastx((r1_path, None), (r1_out_path, None), umilen=4, only_umi=False)

def test_process_fastx_empty_sequences():
    """Test that empty sequences are properly removed."""
    r1_reads = [
        ("read1", "ATTGCGCAATGACGTGGG", "IIIIIIIIIIIIIIIIII"),
        ("read2", "ATTGCGCAATGTTTTGGG", "IIIIIIIIIIIIIIIIII"),
    ]
    r2_reads = [
        ("read1", "GGGACGTAAAGGG", "IIIIIIIIIIIII"),
        ("read2", "GGGTTTTAAA", "IIIIIIIIII"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        r1_path = os.path.join(tmpdir, "r1.fastq.gz")
        r2_path = os.path.join(tmpdir, "r2.fastq.gz")
        r1_out_path = os.path.join(tmpdir, "r1_out.fastq.gz")
        r2_out_path = os.path.join(tmpdir, "r2_out.fastq.gz")

        write_fastq(r1_path, r1_reads)
        write_fastq(r2_path, r2_reads)

        process_fastx((r1_path, r2_path), (r1_out_path, r2_out_path), umilen=4, only_umi=False)

        with pysam.FastxFile(r1_out_path) as f:
            out_reads = [entry for entry in f]

        assert len(out_reads) == 0  # Both sequences should be removed

if __name__ == "__main__":
    pytest.main()

