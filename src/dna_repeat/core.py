# src/dna_repeat/core.py
# Core functions

from Bio import SeqIO
from pathlib import Path
from typing import Iterator
from dataclasses import dataclass
from dna_repeat.constants import COMPLEMENT
from dna_repeat.error import InvalidFASTAError, EmptySequence, InvalidSequence, InvalidKmer
import re

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
    orientation: str

def iter_fasta(input_filepath: Path) -> Iterator[tuple[str, str]]:
    '''Generator that yields (record ID, sequence) tuples from a FASTA file'''
    records = SeqIO.parse(handle=input_filepath, format='fasta')

    try:
        first_record = next(records)
    except StopIteration:
        raise InvalidFASTAError(
            message=f"FASTA file {input_filepath} contains no records and/or is not formatted correctly.",
            details="StopIteration exception in core.py: iter_fasta()",
        )
    yield first_record.id, str(first_record.seq)
    for rec in records:
        yield rec.id, str(rec.seq)

def hamming_distance(a: str, b: str) -> int:
    '''Returns the hamming distance for two strings'''
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def find_repeats(rec_id: str, seq: str, kmer_length: int, allowed_mismatches: int) -> list[RepeatHit]:
    '''Find direct repeats by string slicing and nested loop'''
    if len(seq) < kmer_length:
        return []
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
                    kmer_length = kmer_length,
                    orientation = 'direct'
                    ))
    return hits

def find_invert_repeats(rec_id: str, seq: str, kmer_length: int, allowed_mismatches: int) -> list[RepeatHit]:
    '''Find indirect repeats by string slicing and nested loop'''
    if len(seq) < kmer_length:
        return []
    seq_rc = reverse_complement(seq)
    hits: list[RepeatHit] = []
    for i in range(len(seq) - kmer_length + 1):
        a = seq[i:kmer_length + i]
        for j in range(i, len(seq) - kmer_length + 1):
            b = seq_rc[j:kmer_length + j]
            mismatches = hamming_distance(a, b)
            if mismatches <= allowed_mismatches:
                hits.append(RepeatHit(
                    record_id = rec_id,
                    query_start = i + 1,
                    query_end = i + kmer_length,
                    subject_start = len(seq) - j - kmer_length + 1,
                    subject_end = len(seq) - j,
                    query_seq = a,
                    subject_seq = reverse_complement(b),
                    mismatches = mismatches,
                    kmer_length = kmer_length,
                    orientation = 'inverted'
                    ))
    return hits

def reverse_complement(s: str) -> str:
    '''Get the reverse complement of a dna sequence (no wobbles)'''
    return s.translate(COMPLEMENT)[::-1]

def clean_and_check(rec_id: str, seq: str, kmer_length: int) -> tuple[str, str]:
    clean_seq = re.sub(r'[^A-Z]', '', seq).upper()          # Remove non-alphas and spaces. Make uppercase.
    if kmer_length > len(clean_seq):                        # Check that k-mer lengh isn't longer than the query sequence
        raise InvalidKmer(rec_id)
    if rec_id == '':                                        # Handle no-name seqs, not a deal breaker
        clean_rec_id = f'No-name sequence ({len(clean_seq)} bp)'
    else:
        clean_rec_id = rec_id
    if clean_seq == '':                                     # Check that there's a sequence
        raise EmptySequence(rec_id)
    invalid_chars = []
    for match in re.finditer(r'[^ACGT]', clean_seq):        # Check for invalid characters e.g. wobble bases, etc.
        char = match.group(0)
        pos = match.start() + 1                             # The position of the base is the start position plus one character
        invalid_chars.append(f"'{char}' at pos. {pos}")
    if invalid_chars:
        raise InvalidSequence(rec_id, invalid_chars)
    return clean_rec_id, clean_seq 