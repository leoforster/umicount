import pytest
import pickle
import tempfile
import shutil
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

@pytest.fixture
def temp_pkl_files():
    """Creates a temporary directory with pickle files and yields their paths."""
    tmpdir = tempfile.mkdtemp()
    pkl_files = []
    try:
        for fname, data in sample_filecounts.items():
            path = os.path.join(tmpdir, fname + '.pkl')
            with open(path, 'wb') as pf:
                pickle.dump(data, pf)
            pkl_files.append(path)
        yield pkl_files
    finally:
        shutil.rmtree(tmpdir)

def test_basic_tsv(temp_pkl_files):
    geneorder = ['gene1', 'gene2', 'gene3']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            pkl_files=temp_pkl_files,
            col='A',
            outdir=tmpdir,
            geneorder=geneorder,
            sep='\t'
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.A.tsv'))
        expected = [
            'samples\t_unmapped\t_multimapping\t_no_feature\t_ambiguous\tgene1\tgene2\tgene3',
            'file1\t0\t0\t0\t0\t1\t3\t0',
            'file2\t0\t0\t0\t0\t4\t0\t0'
        ]
        assert lines == expected

def test_basic_csv(temp_pkl_files):
    geneorder = ['gene1', 'gene2']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            pkl_files=temp_pkl_files,
            col='B',
            outdir=tmpdir,
            geneorder=geneorder,
            sep=','
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.B.csv'))
        expected = [
            'samples,_unmapped,_multimapping,_no_feature,_ambiguous,gene1,gene2',
            'file1,0,0,0,0,2,0',
            'file2,0,0,0,0,0,0'
        ]
        assert lines == expected

def test_txt_extension_for_other_sep(temp_pkl_files):
    geneorder = ['gene1']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            pkl_files=temp_pkl_files,
            col='B',
            outdir=tmpdir,
            geneorder=geneorder,
            sep='|'
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.B.txt'))
        expected = [
            'samples|_unmapped|_multimapping|_no_feature|_ambiguous|gene1',
            'file1|0|0|0|0|2',
            'file2|0|0|0|0|0'
        ]
        assert lines == expected

def test_empty_inputs():
    # No files, no genes
    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            pkl_files=[], 
            col='A', 
            outdir=tmpdir, 
            geneorder=[], 
            sep='\t')

        lines = read_lines(os.path.join(tmpdir, 'umicounts.A.tsv'))
        assert lines == ['samples\t_unmapped\t_multimapping\t_no_feature\t_ambiguous']

def test_nonexistent_gene_defaults_to_zero(temp_pkl_files):
    # geneorder includes a gene not in any filecounts
    geneorder = ['geneX']

    with tempfile.TemporaryDirectory() as tmpdir:
        write_counts_for_col(
            pkl_files=temp_pkl_files,
            col='A',
            outdir=tmpdir,
            geneorder=geneorder,
            sep='\t'
        )
        lines = read_lines(os.path.join(tmpdir, 'umicounts.A.tsv'))
        expected = [
            'samples\t_unmapped\t_multimapping\t_no_feature\t_ambiguous\tgeneX',
            'file1\t0\t0\t0\t0\t0',
            'file2\t0\t0\t0\t0\t0'
        ]
        assert lines == expected
