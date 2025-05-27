import pytest
import tempfile
import os

from umicount.umicount import (
    write_counts_for_col
)

def read_lines(path):
    with open(path, 'r') as f:
        return f.read().splitlines()

sample_filecounts = {
        'file1': {
            'gene1': {'A': 1, 'B': 2},
            'gene2': {'A': 3}
        },
        'file2': {
            'gene1': {'A': 4},
            'gene3': {'B': 5}
        }
    }

def test_basic_tsv():
    geneorder = ['gene1', 'gene2', 'gene3']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            filecounts=sample_filecounts,
            col='A',
            outdir=tmpdir,
            geneorder=geneorder,
            sep='\t'
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.A.tsv'))
        expected = [
            'feature\tfile1\tfile2',
            'gene1\t1\t4',
            'gene2\t3\t0',
            'gene3\t0\t0'
        ]
        assert lines == expected

def test_basic_csv():
    geneorder = ['gene1', 'gene2']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            filecounts=sample_filecounts,
            col='B',
            outdir=tmpdir,
            geneorder=geneorder,
            sep=','
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.B.csv'))
        expected = [
            'feature,file1,file2',
            'gene1,2,0',
            'gene2,0,0'
        ]
        assert lines == expected

def test_txt_extension_for_other_sep():
    geneorder = ['gene1']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            filecounts=sample_filecounts,
            col='B',
            outdir=tmpdir,
            geneorder=geneorder,
            sep='|'
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.B.txt'))
        expected = [
            'feature|file1|file2',
            'gene1|2|0'
        ]
        assert lines == expected

def test_empty_inputs():
    # No files, no genes
    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            filecounts={}, 
            col='A', 
            outdir=tmpdir, 
            geneorder=[], 
            sep='\t')

        lines = read_lines(os.path.join(tmpdir, 'umicounts.A.tsv'))
        assert lines == ['feature']

def test_nonexistent_gene_defaults_to_zero():
    # geneorder includes a gene not in any filecounts
    geneorder = ['geneX']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            filecounts=sample_filecounts,
            col='A',
            outdir=tmpdir,
            geneorder=geneorder,
            sep='\t'
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.A.tsv'))
        expected = [
            'feature\tfile1\tfile2',
            'geneX\t0\t0',
        ]
        assert lines == expected
