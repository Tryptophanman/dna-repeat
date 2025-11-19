from dna_repeat.core import clean_and_check
from dna_repeat.error import EmptySequenceError, InvalidKmerError, InvalidSequenceError
import pytest

test_dict = {
    'seq_ok_caps': 'AAGTTGATGTTAAGATGGACAAGAATGTAACTTGAACAAAAGCTGAATCATCTCTTCAGCCACTAGTATG',
    'seq_ok_lowercase': 'agccgtagcagGATCTgcgcgcgatattttttcgagacGGG',
    'seq_invalid_chars': 'AAGTTGATGOTAAGATGGACAAGAATGT?ACTTGAASTEVEAAAGCTGAATCAECTCTTCAGCCACTA6TATG',
    'seq_spaces': 'AAGATGG   AGAATGTAACTTGAA   CAAAAGCTGAATC    ATCTCTT',
    'seq_empty': ' ',
    'seq_short': 'GT',
    'seq_protein': 'MHHHHHHTIPPEW**'
    }
kmer_length = 20

# SUCCESS
def test_all_caps():
    cleaned_id, cleaned_seq = clean_and_check('ok_caps', test_dict['seq_ok_caps'], kmer_length)
    assert cleaned_id == 'ok_caps'
    assert cleaned_seq == 'AAGTTGATGTTAAGATGGACAAGAATGTAACTTGAACAAAAGCTGAATCATCTCTTCAGCCACTAGTATG'

def test_all_caps_but_no_name():
    cleaned_id, cleaned_seq = clean_and_check('', test_dict['seq_ok_caps'], kmer_length)
    assert cleaned_id == f'No-name sequence ({len(cleaned_seq)} bp)'
    assert cleaned_seq == 'AAGTTGATGTTAAGATGGACAAGAATGTAACTTGAACAAAAGCTGAATCATCTCTTCAGCCACTAGTATG'

def test_some_lowercase():
    cleaned_id, cleaned_seq = clean_and_check('ok_lowercase', test_dict['seq_ok_lowercase'], kmer_length)
    assert cleaned_seq == 'AGCCGTAGCAGGATCTGCGCGCGATATTTTTTCGAGACGGG'

def test_some_spaces():
    cleaned_id, cleaned_seq = clean_and_check('ok_spaces', test_dict['seq_spaces'], kmer_length)
    assert cleaned_seq == 'AAGATGGAGAATGTAACTTGAACAAAAGCTGAATCATCTCTT'

# RAISED
def test_empty_Sequence_error():
# Error -> empty sequence
    with pytest.raises(EmptySequenceError):
        cleaned_id, cleaned_seq = clean_and_check('empty', test_dict['seq_empty'], kmer_length)

def test_invalid_chars_error():
    with pytest.raises(InvalidSequenceError):
        cleaned_id, cleaned_seq = clean_and_check('invalid_chars', test_dict['seq_invalid_chars'], kmer_length)
        cleaned_id, cleaned_seq = clean_and_check('protein', test_dict['seq_protein'], kmer_length)

def test_invalid_kmer():
    with pytest.raises(InvalidKmerError):
        cleaned_id, cleaned_seq = clean_and_check('short', test_dict['seq_short'], kmer_length)