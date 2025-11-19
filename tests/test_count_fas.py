from dna_repeat.core import count_fasta_seqs
from pathlib import Path
import pytest
from dna_repeat.error import InvalidFASTAError

def test_count_fast_seqs():

    number_of_seqs = count_fasta_seqs(Path('tests/test.fasta'))     # fas file with 4 seqs
    print(f'\n{number_of_seqs}')
    assert number_of_seqs == 4

def test_broken_fas_file():
    with pytest.raises(InvalidFASTAError):
        number_of_seqs = count_fasta_seqs(Path('tests/bad.fasta'))  # fas file with some text, no seqs
        print(number_of_seqs)   # dummy line to give use to variable to avoid "problem"