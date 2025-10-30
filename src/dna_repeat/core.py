# src/dna_repeat/core.py
    # Core functions

from Bio import SeqIO
from pathlib import Path
from typing import Iterator
from dataclasses import dataclass
from dna_repeat.constants import COMPLEMENT

@dataclass
class RepeatHit:
    '''Stores information for a found repeat'''
    record_id: str
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    query_seq: str
    subject_seq: str
    mismatches: int
    kmer_length: int

def iter_fasta(input_filepath: Path) -> Iterator[tuple[str, str]]:
    '''Generator that yields (record ID, sequence) tuples from a FASTA file'''
    for rec in SeqIO.parse(handle=input_filepath, format='fasta'):
        yield rec.id, str(rec.seq)

def hamming_distance(a: str, b: str) -> int:
    '''Returns the hamming distance for two strings'''
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def find_repeats(rec_id: str, seq: str, kmer_length: int, allowed_mismatches: int) -> list[RepeatHit]:
    hits: list[RepeatHit] = []
    for i in range(len(seq) - kmer_length + 1):
        a = seq[i:kmer_length + i]
        for j in range(i + 1, len(seq) - kmer_length + 1):
            b = seq[j:kmer_length + j]
            mismatches = hamming_distance(a, b)
            if mismatches <= allowed_mismatches:
                hits.append(RepeatHit(
                    record_id = rec_id,
                    query_start = i + 1,        # '+ 1' = convert to 1-based index for positions (DNA)
                    query_end = i + kmer_length,
                    subject_start = j + 1,
                    subject_end = j + kmer_length,
                    query_seq = a,
                    subject_seq = b,
                    mismatches = mismatches,
                    kmer_length = kmer_length
                    ))
    return hits

def reverse_complement(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]