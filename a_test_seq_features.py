import pytest
import seq_features

def test_number_negatives_single_E_or_D():
    """Perform unit tests on number_negative for single AA"""
    assert seq_features.number_negatives('E') == 1
    assert seq_features.number_negatives('D') == 1

def test_number_negatives_for_empty():
    """Perform unit tests on number_negative for empty entry"""
    assert seq_features.number_negatives('') == 0

def test_number_negatives_for_short_sequences():
    """Perform unit tests on number_negative for short sequence"""
    assert seq_features.number_negatives('ACKLWTTAE') == 1
    assert seq_features.number_negatives('DDDDEEEE') == 8

    # TESTS NEGATIVES FOR PH DEPENDENCE
def test_low_pH():
    assert seq_features.number_negatives('ACDDEK', pH = 3) == 0
    assert seq_features.number_negatives('ACDDEK', pH = 5) == 3

def test_number_negatives_for_lowercase():
    """Perform unit tests on number_negative for lowercase"""
    assert seq_features.number_negatives('acklwttae') == 1

def test_number_negatives_for_invalid_amino_acid():
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.number_negatives('Z')
    excinfo.match("Z is not a valid amino acid")

def test_number_positives_single_R_K_or_H():
    """Perform unit tests on number_positives for single AA"""
    assert seq_features.number_positives('R') == 1
    assert seq_features.number_positives('K') == 1
    assert seq_features.number_positives('H') == 1

def test_number_positives_for_empty():
    """Perform unit tests on number_positives for empty entry"""
    assert seq_features.number_positives('') == 0

def test_number_positives_for_short_sequences():
    """Perform unit tests on number_positives for short sequence"""
    assert seq_features.number_positives('RCKLWTTRE') == 3
    assert seq_features.number_positives('DDDDEEEE') == 0
    
    # TESTS POSITIVES FOR PH DEPENDENCE
def test_low_pH():
    assert seq_features.number_positives('ACKLHGR', pH = 7.7) == 3
    assert seq_features.number_positives('ACKLHGA', pH = 9) == 1
    assert seq_features.number_positives('ACKLHGR', pH = 11) == 0

def test_number_positives_for_lowercase():
    """Perform unit tests on number_positives for lowercase"""
    assert seq_features.number_positives('rcklwttre') == 3
    
def test_number_positives_for_invalid_amino_acid():
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.number_positives('Z')
    excinfo.match("Z is not a valid amino acid")

def test_number_negatives_for_invalid_amino_acid_anywhere():
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.number_negatives('AZK')
    excinfo.match("Z is not a valid amino acid")
    
def test_number_positives_for_invalid_amino_acid_anywhere():
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.number_positives('AZK')
    excinfo.match("Z is not a valid amino acid")

def test_is_valid_sequence_for_invalid_amino_acid():
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.is_valid_sequence('Z')
    excinfo.match("Z is not a valid amino acid")    
    
def test_is_valid_sequence_for_invalid_amino_acid_anywhere():
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.is_valid_sequence('AZK')
    excinfo.match("Z is not a valid amino acid")

def test_is_valid_sequence_for_other_invalid_amino_acid_anywhere():
    assert seq_features.is_valid_sequence('ALKSAYGS') is None
    
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.is_valid_sequence('AZLL')
    excinfo.match("Z is not a valid amino acid")
    
    with pytest.raises(RuntimeError) as excinfo:
        seq_features.is_valid_sequence('ALLBJ')
    excinfo.match("B is not a valid amino acid")

    with pytest.raises(RuntimeError) as excinfo:
        seq_features.is_valid_sequence('AL%J')
    excinfo.match("% is not a valid amino acid")

    #tests net_charge
def test_net_charge():
    assert seq_features.net_charge('ACCAVDDKL') == -1
    assert seq_features.net_charge('ACAHHDKLH', pH = 9) == 0
    assert seq_features.net_charge('ACAHHDKLH', pH = 11) == -1
    assert seq_features.net_charge('ACAHHDKLH', pH = 2) == 4