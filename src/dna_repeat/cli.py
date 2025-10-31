# src/dna_repeat/cli.py
# Top-level functions for dna-repeat CLI

import sys
import argparse
# import csv
from dna_repeat import __version__
from dna_repeat.core import iter_fasta, find_repeats, find_invert_repeats, RepeatHit
from pathlib import Path
import pandas as pd

def init_argparser() -> argparse.ArgumentParser:
    '''Initializes argument parser for CLI'''
    parser = argparse.ArgumentParser(prog='dna-repeat', description=f'dna-repeat {__version__} - Finds repeats in a DNA sequence')
    parser.add_argument("input_filepath", help='filepath to input file (FASTA)')
    parser.add_argument("-o", "--output", help='desired directory for output file (output.csv). Default: cwd', nargs='?', const='.', default=None, dest='output_directory')
    parser.add_argument("-k", "--length", help='repeat length in bp (min 4, max 30)', type=int, default=20, dest='kmer_length')
    parser.add_argument("-m", "--mismatches", help='number of bp mismatches allowed', type=int, default=0, dest='allowed_mismatches')
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inverted-only", help='only look for inverted repeats', action='store_true', dest='inverted_only')
    group.add_argument("-d", "--direct-only", help='only look for direct repeats', action='store_true', dest='direct_only')
    return parser

def main(argv: str | None = None) -> int:
    parser = init_argparser()
    args = parser.parse_args(argv)

    input_filepath = Path(args.input_filepath).resolve()
    if args.output_directory:
        output_directory = Path(args.output_directory).resolve()
        output_directory.mkdir(parents=True, exist_ok=True)
        output_filepath: Path | None = output_directory / 'output.csv'
    else:
        output_filepath = None
    kmer_length: int = args.kmer_length
    allowed_mismatches: int = args.allowed_mismatches

    do_direct = not args.inverted_only
    do_inverted = not args.direct_only

    if not input_filepath.is_file():
        print(f'Input file not found: {input_filepath}', file=sys.stderr)
        return 1
    if not 4 <= kmer_length <= 30:
        print('repeat length (-k, --length) must be in the range 4-30', file=sys.stderr)
        return 1
    if allowed_mismatches > kmer_length / 2:
        print('m must be <= kmer_length / 2', file=sys.stderr)
        return 1       

    results = []

    for rec_id, seq in iter_fasta(input_filepath):
        if do_direct:
            repeats: list[RepeatHit] = find_repeats(
                rec_id=rec_id,
                seq=seq,
                kmer_length=kmer_length,
                allowed_mismatches=allowed_mismatches,
            )
            if repeats:
                results.extend(repeats)
        if do_inverted:
            repeats: list[RepeatHit] = find_invert_repeats(
                rec_id=rec_id,
                seq=seq,
                kmer_length=kmer_length,
                allowed_mismatches=allowed_mismatches,
            )
            if repeats:
                results.extend(repeats)    
    if results:
        df = pd.DataFrame(results)
    else:
        df = pd.DataFrame(columns=list(RepeatHit.__dataclass_fields__.keys()))  # empty results
    df.to_csv(sys.stdout if output_filepath is None else output_filepath, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
