# src/dna_repeat/core.py
    # Core functions

from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord         # import the class that defines the attributes
from pathlib import Path
from typing import Iterator

def read_fasta(input_filepath: Path) -> Iterator[tuple[str, str]]:
    '''Iterator for sequences in fasta file'''
    for rec in SeqIO.parse(handle=input_filepath, format='fasta'):
        yield rec.id, str(rec.seq)

def hamming_distance(a: str, b: str) -> int:
    '''Returns the hamming distance for two strings'''
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def find_repeats(seq: str, rep_length: int, allowed_mismatches: int) -> list[tuple[str, str, int]]:
    '''Brute force nested loop to search for repeats (strings)'''
    hits = []
    for i in range(len(seq) - rep_length + 1):
        a = seq[i:rep_length + i]
        for j in range(i + 1, len(seq) - rep_length + 1):
            b = seq[j:rep_length + j]
            mismatches = hamming_distance(a, b)
            if mismatches <= allowed_mismatches:
                hits.append((a, b, mismatches))
    return hits
