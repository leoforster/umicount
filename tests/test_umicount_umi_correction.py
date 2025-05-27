import pytest

from umicount.umicount import umi_correction

############################
# Tests for UMI correction
############################

def test_empty_dict():
    """
    Test empty dictionary input
    """
    result = umi_correction({}, 
                            countratio=2, 
                            hamming_threshold=1)
    assert result == {}

def test_umi_length_mismatch():
    """
    Test mixed length UMIs
    """
    with pytest.raises(ValueError): # from RapidFuzz
        umi_correction({'AAAA': 10, 'AAA': 1}, 
                        countratio=2,
                        hamming_threshold=1)

def test_zero_counts():
    """
    Test input with 0-count UMI
    """
    with pytest.raises(ValueError):
        umi_correction({'AAAA': 5, 'AAAT': 2, 'ATAA':0}, 
                        countratio=2,
                        hamming_threshold=1)

def test_single_umi():
    """
    Test single UMI input
    """
    input_dict = {'AAAA': 10}
    result = umi_correction(input_dict, 
                            countratio=2, 
                            hamming_threshold=1)
    assert result == {'AAAA': 10}

def test_no_merge():
    """
    Test distant UMIs merging
    """
    input_dict = {'AAAA': 10, 'CCCC': 5}
    result = umi_correction(input_dict, 
                            countratio=2, 
                            hamming_threshold=1)
    assert result == input_dict # no merge made

def test_basic_merge():
    """
    Test basic UMI merge
    """
    input_dict = {'AAAA': 20, 'AAAT': 5}
    result = umi_correction(input_dict, 
                            countratio=2, 
                            hamming_threshold=1)
    assert result == {'AAAA': 25}

def test_multiple_merges():
    """
    Test merging multiple UMIs
    """
    input_dict = {
        'AAAA': 20,   # seed candidate
        'AAAT': 5,    # 1-Hamming from 'AAAA'
        'AACA': 7,    # 1-Hamming from 'AAAA'
        'AATA': 3,    # 1-Hamming from 'AAAA'
        'CCCC': 15    # Not related (Hamming distance 4)
    }
    result = umi_correction(input_dict, 
                            countratio=2, 
                            hamming_threshold=1)
    assert result == {'AAAA': 35, 'CCCC': 15}

def test_countsratio_threshold():
    """
    Test merging multiple high-count UMIs
    """
    input_dict = {
        'AAAA': 20, # seed candidate
        'ATAA': 10, # passes seed>(c*count)-1 at c=2, not 3
        'AATA': 5,  # always merged
    }
    result = umi_correction(input_dict, 
                            countratio=2, # lower merging threshold 
                            hamming_threshold=1)
    assert result == {'AAAA': 35}

    result = umi_correction(input_dict, 
                            countratio=3, # higher merging threshold
                            hamming_threshold=1)
    assert result == {'AAAA': 25, 'ATAA':10}

def test_chain_merging():
    """
    Test that UMIs are not chain-merged
    """
    input_dict = {
        'AAAA': 20,   # seed candidate
        'AAAT': 5,    # 1-Hamming from 'AAAA'
        'AATT': 2,    # 1-Hamming from 'AAAT' (and from 'AAAA' if chain merging)
        'TTTT': 10    # Unrelated
    }
    result = umi_correction(input_dict, 
                            countratio=2,
                            hamming_threshold=1)
    assert result == {'TTTT': 10, 'AAAA': 25, 'AATT': 2}
